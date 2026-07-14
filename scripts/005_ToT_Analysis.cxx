#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TTree.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSamplePayloadLineCount = 38;
constexpr std::size_t kSampleElinkCount = 4;
constexpr std::size_t kValuesPerSample = kSamplePayloadLineCount * kSampleElinkCount;

struct PlotContext {
    std::string input_path;
    std::string run_arguments;
    std::string run_time;
};

std::string file_name_only(const std::string& path)
{
    return std::filesystem::path(path).filename().string();
}

std::string format_count(std::size_t value)
{
    std::string digits = std::to_string(value);
    for (int insert_position = static_cast<int>(digits.size()) - 3; insert_position > 0; insert_position -= 3) {
        digits.insert(static_cast<std::size_t>(insert_position), ",");
    }
    return digits;
}

std::string format_run_time()
{
    const auto now = std::chrono::system_clock::now();
    const auto now_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time{};
    localtime_r(&now_time, &local_time);
    std::ostringstream stream;
    stream << std::put_time(&local_time, "%Y-%m-%d %H:%M:%S");
    return stream.str();
}

std::string command_line(int argc, char** argv)
{
    std::ostringstream stream;
    stream << "command";
    for (int index = 0; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--input" || argument == "-i") {
            ++index;
            continue;
        }
        if (argument.rfind("--input=", 0) == 0 || argument.rfind("-i=", 0) == 0) {
            continue;
        }
        stream << ' ' << argv[index];
    }
    return stream.str();
}

void draw_page_header(const PlotContext& context, const std::string& plot_title)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(13);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.14, 0.88, "#bf{FoCal TB2026 July}");
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.14, 0.835, plot_title.c_str());
    latex.DrawLatex(0.14, 0.795, file_name_only(context.input_path).c_str());
    latex.DrawLatex(0.14, 0.755, context.run_arguments.c_str());
    latex.DrawLatex(0.14, 0.715, context.run_time.c_str());
}

std::filesystem::path default_output_path(const std::string& input_path)
{
    const auto input_file = std::filesystem::path(input_path);
    auto parent_name = input_file.parent_path().filename();
    if (parent_name.empty()) {
        parent_name = input_file.stem();
    }
    return std::filesystem::path("dump") / "005" / parent_name / "tot_val1_accumulated_waveform.pdf";
}

std::filesystem::path tot_sum_histogram_root_path(const std::string& output_path)
{
    return std::filesystem::path(output_path).parent_path() / "event_tot_sum_histograms.root";
}

constexpr const char* kSampleEventTreeName = "FOCAL";
constexpr const char* kLegacySampleEventTreeName = "sample_events";
constexpr const char* kSampleEventBranchPrefix = "FOCALHCAL.";

std::string g_sample_event_branch_prefix = kSampleEventBranchPrefix;

std::string sample_event_branch_name(const char* member_name)
{
    return g_sample_event_branch_prefix + member_name;
}

void use_only_branches(TTree& tree, const std::vector<std::string>& branch_names)
{
    tree.SetBranchStatus("*", 0);
    for (const auto& branch_name : branch_names) {
        tree.SetBranchStatus(branch_name.c_str(), 1);
    }
}

void restore_tree_read_state(TTree& tree)
{
    tree.ResetBranchAddresses();
    tree.SetBranchStatus("*", 1);
}

struct AnalysisLayout {
    std::size_t sample_bins = 1;
    std::vector<std::uint8_t> active_lpgbts;
};

struct LpgbtAnalysis {
    std::vector<std::unique_ptr<TH2D>> channel_histograms;
    std::vector<double> tot_sums;
};

struct AnalysisProducts {
    std::unique_ptr<TH2D> waveform_histogram;
    std::map<std::uint8_t, LpgbtAnalysis> lpgbts;
};

AnalysisLayout scan_analysis_layout(TTree& tree)
{
    use_only_branches(tree, {sample_event_branch_name("n_samples"), sample_event_branch_name("lpgbt_id")});
    std::vector<std::uint32_t>* n_samples = nullptr;
    std::vector<std::uint8_t>* lpgbt_id = nullptr;
    tree.SetBranchAddress(sample_event_branch_name("n_samples").c_str(), &n_samples);
    tree.SetBranchAddress(sample_event_branch_name("lpgbt_id").c_str(), &lpgbt_id);

    AnalysisLayout layout;
    std::size_t max_count = 0;
    const auto entries = tree.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree.GetEntry(entry);
        if (n_samples != nullptr) {
            for (const auto count : *n_samples) {
                max_count = std::max(max_count, static_cast<std::size_t>(count));
            }
        }
        if (lpgbt_id != nullptr) {
            for (const auto id : *lpgbt_id) {
                if (std::find(layout.active_lpgbts.begin(), layout.active_lpgbts.end(), id) == layout.active_lpgbts.end()) {
                    layout.active_lpgbts.push_back(id);
                }
            }
        }
    }
    layout.sample_bins = max_count == 0 ? 1 : max_count;
    std::sort(layout.active_lpgbts.begin(), layout.active_lpgbts.end());
    restore_tree_read_state(tree);
    return layout;
}

std::vector<std::unique_ptr<TH2D>> make_lpgbt_channel_histograms(std::uint8_t id, std::size_t sample_bins)
{
    std::vector<std::unique_ptr<TH2D>> histograms;
    histograms.reserve(kValuesPerSample);
    for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
        const auto row_index = channel_index / kSampleElinkCount;
        const auto elink_index = channel_index % kSampleElinkCount;
        auto hist = std::make_unique<TH2D>(("h_tot_val1_lpgbt_" + std::to_string(id) + "_ch_" + std::to_string(channel_index)).c_str(),
            ("lpGBT " + std::to_string(id) + " row " + std::to_string(row_index) + " elink " + std::to_string(elink_index) + ";sample index;ToT;count").c_str(),
            static_cast<int>(sample_bins), -0.5, static_cast<double>(sample_bins) - 0.5,
            4096, -0.5, 4095.5);
        histograms.push_back(std::move(hist));
    }
    return histograms;
}

std::uint16_t decode_tot(std::uint16_t raw_value)
{
    return raw_value < 512u ? raw_value : static_cast<std::uint16_t>((raw_value - 512u) * 8u);
}

AnalysisProducts collect_analysis_products(TTree& tree, const AnalysisLayout& layout)
{
    AnalysisProducts products;
    products.waveform_histogram = std::make_unique<TH2D>("h_tot_val1_accumulated_waveform",
        "Accumulated ToT waveform;sample index;ToT;count",
        static_cast<int>(layout.sample_bins), -0.5, static_cast<double>(layout.sample_bins) - 0.5,
        4096, -0.5, 4095.5);
    for (const auto id : layout.active_lpgbts) {
        products.lpgbts[id].channel_histograms = make_lpgbt_channel_histograms(id, layout.sample_bins);
    }

    use_only_branches(tree, {sample_event_branch_name("val1"), sample_event_branch_name("lpgbt_id"), sample_event_branch_name("n_samples"), sample_event_branch_name("payload_offset"), sample_event_branch_name("payload_count")});
    std::vector<std::uint16_t>* val1 = nullptr;
    std::vector<std::uint8_t>* lpgbt_id = nullptr;
    std::vector<std::uint32_t>* n_samples = nullptr;
    std::vector<std::uint32_t>* payload_offset = nullptr;
    std::vector<std::uint32_t>* payload_count = nullptr;
    tree.SetBranchAddress(sample_event_branch_name("val1").c_str(), &val1);
    tree.SetBranchAddress(sample_event_branch_name("lpgbt_id").c_str(), &lpgbt_id);
    tree.SetBranchAddress(sample_event_branch_name("n_samples").c_str(), &n_samples);
    tree.SetBranchAddress(sample_event_branch_name("payload_offset").c_str(), &payload_offset);
    tree.SetBranchAddress(sample_event_branch_name("payload_count").c_str(), &payload_count);

    const auto entries = tree.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree.GetEntry(entry);
        if (val1 == nullptr || lpgbt_id == nullptr || n_samples == nullptr || payload_offset == nullptr || payload_count == nullptr) {
            continue;
        }
        for (std::size_t lpgbt_index = 0; lpgbt_index < lpgbt_id->size() && lpgbt_index < n_samples->size() && lpgbt_index < payload_offset->size() && lpgbt_index < payload_count->size(); ++lpgbt_index) {
            const auto id = (*lpgbt_id)[lpgbt_index];
            const auto begin = static_cast<std::size_t>((*payload_offset)[lpgbt_index]);
            const auto count = static_cast<std::size_t>((*payload_count)[lpgbt_index]);
            const auto sample_count = static_cast<std::size_t>((*n_samples)[lpgbt_index]);
            const auto end = std::min(begin + count, val1->size());
            auto lpgbt_iter = products.lpgbts.find(id);
            for (std::size_t index = begin; index < end; ++index) {
                const auto local_index = index - begin;
                const auto sample_index = static_cast<double>(local_index / kValuesPerSample);
                const auto value = static_cast<double>(decode_tot((*val1)[index]));
                products.waveform_histogram->Fill(sample_index, value);
                if (lpgbt_iter != products.lpgbts.end()) {
                    const auto channel_index = local_index % kValuesPerSample;
                    if (channel_index < lpgbt_iter->second.channel_histograms.size()) {
                        lpgbt_iter->second.channel_histograms[channel_index]->Fill(sample_index, value);
                    }
                }
            }

            if (sample_count == 0 || count < kValuesPerSample || begin + count > val1->size() || lpgbt_iter == products.lpgbts.end()) {
                continue;
            }
            double tot_sum = 0.0;
            for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
                for (std::size_t sample_index = 0; sample_index < sample_count; ++sample_index) {
                    const auto value_index = begin + sample_index * kValuesPerSample + channel_index;
                    if (value_index >= begin + count || value_index >= val1->size()) {
                        break;
                    }
                    const auto raw_value = (*val1)[value_index];
                    if (raw_value != 0) {
                        tot_sum += static_cast<double>(decode_tot(raw_value));
                        break;
                    }
                }
            }
            lpgbt_iter->second.tot_sums.push_back(tot_sum);
        }
    }
    restore_tree_read_state(tree);
    return products;
}

void draw_lpgbt_channel_waveform_pages(const std::map<std::uint8_t, LpgbtAnalysis>& lpgbts, const std::string& output_path)
{
    constexpr int columns = 4;
    constexpr int rows = 4;
    constexpr std::size_t channels_per_page = columns * rows;
    for (const auto& [id, analysis] : lpgbts) {
        const auto& histograms = analysis.channel_histograms;
        for (std::size_t page_start = 0; page_start < histograms.size(); page_start += channels_per_page) {
            const auto page_index = page_start / channels_per_page;
            TCanvas canvas(("tot_lpgbt_channel_waveform_canvas_" + std::to_string(id) + "_page_" + std::to_string(page_index)).c_str(),
                ("ToT val1 channel waveforms for lpGBT " + std::to_string(id)).c_str(), 1600, 1200);
            canvas.Divide(columns, rows, 0.001, 0.001);
            const auto page_end = std::min(page_start + channels_per_page, histograms.size());
            for (std::size_t channel_index = page_start; channel_index < page_end; ++channel_index) {
                canvas.cd(static_cast<int>(channel_index - page_start + 1));
                gPad->SetRightMargin(0.18);
                gPad->SetLeftMargin(0.12);
                gPad->SetBottomMargin(0.13);
                gPad->SetLogz();
                histograms[channel_index]->GetXaxis()->SetLabelSize(0.045);
                histograms[channel_index]->GetYaxis()->SetLabelSize(0.045);
                histograms[channel_index]->GetZaxis()->SetLabelSize(0.045);
                histograms[channel_index]->GetXaxis()->SetTitleSize(0.05);
                histograms[channel_index]->GetYaxis()->SetTitleSize(0.05);
                histograms[channel_index]->Draw("COLZ");
            }
            canvas.Print(output_path.c_str());
        }
    }
}

void draw_tot_sum_histogram_pages(const std::vector<double>& sums, const std::string& output_path, const std::string& name, const std::string& title, const PlotContext& context)
{
    if (sums.empty()) {
        return;
    }
    const auto max_sum = *std::max_element(sums.begin(), sums.end());
    const auto x_max = max_sum <= 0.0 ? 1.0 : max_sum * 1.05;
    TH1D hist(("h_" + name).c_str(), (title + ";#Sigma_{channels} first non-zero ToT;events").c_str(), 120, 0.0, x_max);
    for (const auto sum : sums) {
        hist.Fill(sum);
    }
    hist.SetLineColor(kAzure + 1);
    hist.SetLineWidth(2);
    hist.SetFillStyle(0);

    const auto draw_page = [&](bool log_y) {
        const auto page_title = title + (log_y ? " log-y" : " linear-y");
        TCanvas canvas((name + (log_y ? "_log" : "_linear")).c_str(), page_title.c_str(), 1100, 800);
        canvas.SetLogy(log_y ? 1 : 0);
        canvas.SetLeftMargin(0.12);
        canvas.SetRightMargin(0.08);
        canvas.SetBottomMargin(0.12);
        hist.SetMinimum(log_y ? 0.5 : 0.0);
        hist.SetMaximum(std::max(1.0, hist.GetMaximum() * (log_y ? 1.5 : 1.25)));
        hist.Draw("hist");
        draw_page_header(context, page_title);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextAlign(33);
        latex.SetTextSize(0.032);
        latex.DrawLatex(0.88, 0.88, page_title.c_str());
        latex.DrawLatex(0.88, 0.83, ("events = " + format_count(static_cast<std::size_t>(hist.GetEntries()))).c_str());
        canvas.Print(output_path.c_str());
    };

    draw_page(true);
    draw_page(false);
}

std::unique_ptr<TH1D> make_tot_sum_histogram(const std::vector<double>& sums, const std::string& name, const std::string& title)
{
    const auto max_sum = sums.empty() ? 0.0 : *std::max_element(sums.begin(), sums.end());
    const auto x_max = max_sum <= 0.0 ? 1.0 : max_sum * 1.05;
    auto hist = std::make_unique<TH1D>(name.c_str(), (title + ";#Sigma_{channels} first non-zero ToT;events").c_str(), 120, 0.0, x_max);
    for (const auto sum : sums) {
        hist->Fill(sum);
    }
    hist->SetDirectory(nullptr);
    return hist;
}

void write_all_lpgbt_tot_sum_histogram_root(const std::vector<double>& all_lpgbt_sums, const std::filesystem::path& root_path, const PlotContext& context)
{
    std::filesystem::create_directories(root_path.parent_path());
    TFile output(root_path.string().c_str(), "RECREATE");
    if (output.IsZombie()) {
        throw std::runtime_error("failed to create ToT sum histogram ROOT file: " + root_path.string());
    }
    TNamed("ToT_sum_input_file", context.input_path.c_str()).Write();
    TNamed("ToT_sum_histogram_note", "All lpGBT event ToT sum from 005_ToT_Analysis; per channel uses first non-zero val1 sample; val1>=512 maps to (val1-512)*8").Write();
    auto hist = make_tot_sum_histogram(all_lpgbt_sums, "h_event_tot_sum_all_lpgbt", "All lpGBT event ToT sum");
    hist->Write();
    output.Write();
    output.Close();
}

void draw_lpgbt_tot_sum_distribution_pages(const std::map<std::uint8_t, LpgbtAnalysis>& lpgbts, const std::string& output_path, const std::filesystem::path& tot_sum_root_path, const PlotContext& context)
{
    if (lpgbts.empty()) {
        write_all_lpgbt_tot_sum_histogram_root({}, tot_sum_root_path, context);
        return;
    }

    std::vector<double> all_lpgbt_sums;
    const auto event_count = std::max_element(lpgbts.begin(), lpgbts.end(), [](const auto& left, const auto& right) {
        return left.second.tot_sums.size() < right.second.tot_sums.size();
    })->second.tot_sums.size();
    all_lpgbt_sums.assign(event_count, 0.0);
    for (const auto& [id, analysis] : lpgbts) {
        (void)id;
        const auto& sums = analysis.tot_sums;
        for (std::size_t index = 0; index < sums.size() && index < all_lpgbt_sums.size(); ++index) {
            all_lpgbt_sums[index] += sums[index];
        }
    }
    for (const auto& [id, analysis] : lpgbts) {
        draw_tot_sum_histogram_pages(analysis.tot_sums, output_path, "event_tot_sum_lpgbt_" + std::to_string(id), "lpGBT " + std::to_string(id) + " event ToT sum", context);
    }
    write_all_lpgbt_tot_sum_histogram_root(all_lpgbt_sums, tot_sum_root_path, context);
    draw_tot_sum_histogram_pages(all_lpgbt_sums, output_path, "event_tot_sum_all_lpgbt", "All lpGBT event ToT sum", context);
}

void draw_waveform(const std::string& input_path, const std::string& output_path, const PlotContext& context)
{
    TFile input(input_path.c_str(), "READ");
    if (input.IsZombie()) {
        throw std::runtime_error("failed to open input ROOT file: " + input_path);
    }

    auto* tree = dynamic_cast<TTree*>(input.Get(kSampleEventTreeName));
    if (tree == nullptr) {
        tree = dynamic_cast<TTree*>(input.Get(kLegacySampleEventTreeName));
        g_sample_event_branch_prefix.clear();
    } else {
        g_sample_event_branch_prefix = kSampleEventBranchPrefix;
    }
    if (tree == nullptr) {
        throw std::runtime_error("input ROOT file does not contain TTree 'FOCAL' or legacy TTree 'sample_events': " + input_path);
    }
    if (tree->GetBranch(sample_event_branch_name("val1").c_str()) == nullptr) {
        throw std::runtime_error("sample event tree does not contain branch '" + sample_event_branch_name("val1") + "': " + input_path);
    }
    if (tree->GetBranch(sample_event_branch_name("n_samples").c_str()) == nullptr || tree->GetBranch(sample_event_branch_name("payload_offset").c_str()) == nullptr || tree->GetBranch(sample_event_branch_name("payload_count").c_str()) == nullptr) {
        throw std::runtime_error("sample event tree is missing sample layout branches under '" + g_sample_event_branch_prefix + "': " + input_path);
    }
    if (tree->GetBranch(sample_event_branch_name("lpgbt_id").c_str()) == nullptr) {
        throw std::runtime_error("sample event tree does not contain branch '" + sample_event_branch_name("lpgbt_id") + "': " + input_path);
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const auto layout = scan_analysis_layout(*tree);
    auto products = collect_analysis_products(*tree, layout);

    std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
    const auto tot_sum_root_path = tot_sum_histogram_root_path(output_path);
    TCanvas canvas("tot_waveform_canvas", "ToT val1 accumulated waveform", 1100, 800);
    canvas.SetRightMargin(0.16);
    canvas.SetLeftMargin(0.11);
    canvas.SetBottomMargin(0.11);
    canvas.SetLogz();
    products.waveform_histogram->GetZaxis()->SetTitle("count");
    products.waveform_histogram->Draw("COLZ");
    draw_page_header(context, "Accumulated ToT waveform");
    canvas.Print((output_path + "[").c_str());
    canvas.Print(output_path.c_str());
    draw_lpgbt_tot_sum_distribution_pages(products.lpgbts, output_path, tot_sum_root_path, context);
    draw_lpgbt_channel_waveform_pages(products.lpgbts, output_path);
    canvas.Print((output_path + "]").c_str());

    spdlog::info("Read {} entries from {}", tree->GetEntries(), input_path);
    spdlog::info("Wrote ToT waveform PDF to {}", output_path);
    spdlog::info("Wrote ToT sum histogram ROOT file to {}", tot_sum_root_path.string());
}

} // namespace

int main(int argc, char** argv)
{
    const std::string script_full_name = argv[0];
    const auto script_name = script_full_name.substr(script_full_name.find_last_of("/\\") + 1);
    spdlog::info("Running script: {}", script_name);
    spdlog::info("----------------------------------------");

    cxxopts::Options options(script_name, "Draw accumulated val1 ToT waveform and event ToT sums from 001_Rootifier sample_events.root FOCAL tree");
    options.add_options()
        ("i,input", "Input sample_events.root file containing FOCAL/FOCALHCAL branches", cxxopts::value<std::string>())
        ("o,output", "Output PDF path", cxxopts::value<std::string>())
        ("h,help", "Print help");

    cxxopts::ParseResult parsed_opts;
    try {
        parsed_opts = options.parse(argc, argv);
    } catch (const std::exception& error) {
        spdlog::error("Argument error: {}", error.what());
        std::cout << options.help() << std::endl;
        return 2;
    }

    if (parsed_opts.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (!parsed_opts.count("input")) {
        spdlog::error("Input file is required. Use --input <sample_events.root>.");
        std::cout << options.help() << std::endl;
        return 2;
    }

    try {
        const auto input_path = parsed_opts["input"].as<std::string>();
        const auto output_path = parsed_opts.count("output")
            ? parsed_opts["output"].as<std::string>()
            : default_output_path(input_path).string();
        const PlotContext context{input_path, command_line(argc, argv), format_run_time()};
        draw_waveform(input_path, output_path, context);
    } catch (const std::exception& error) {
        spdlog::error("{}", error.what());
        return 1;
    }

    return 0;
}
