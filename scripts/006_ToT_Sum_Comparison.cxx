#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct InputHistogram {
    std::string input_path;
    std::string label;
    double parameter_value = 0.0;
    bool has_parameter_value = false;
    std::unique_ptr<TH1D> hist;
};

struct ConfigEntry {
    std::string raw_path;
    std::string histogram_path;
    std::string label;
    double parameter_value = 0.0;
};

struct ComparisonConfig {
    std::string title = "All lpGBT ToT sum comparison";
    std::string parameter_name;
    std::string parameter_unit;
    std::string output_path;
    std::vector<ConfigEntry> entries;
};

std::string current_time_minute()
{
    const auto now = std::time(nullptr);
    std::tm tm{};
    localtime_r(&now, &tm);
    std::ostringstream stream;
    stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return stream.str();
}

std::string file_label(const std::string& path)
{
    const auto input = std::filesystem::path(path);
    const auto parent = input.parent_path().filename().string();
    return parent.empty() ? input.stem().string() : parent;
}

std::string read_text_file(const std::string& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open config file: " + path);
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::optional<std::string> extract_json_string(const std::string& text, const std::string& key)
{
    const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        return std::nullopt;
    }
    return match[1].str();
}

std::optional<double> extract_json_number(const std::string& text, const std::string& key)
{
    const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*(-?[0-9]+(?:\\.[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        return std::nullopt;
    }
    return std::stod(match[1].str());
}

std::filesystem::path histogram_path_from_raw_path(const std::string& raw_path)
{
    const auto raw_file = std::filesystem::path(raw_path);
    return std::filesystem::path("dump") / "005" / raw_file.filename() / "event_tot_sum_histograms.root";
}

ComparisonConfig load_config(const std::string& config_path)
{
    const auto text = read_text_file(config_path);
    ComparisonConfig config;
    if (const auto title = extract_json_string(text, "title")) {
        config.title = *title;
    }
    if (const auto parameter_name = extract_json_string(text, "parameter_name")) {
        config.parameter_name = *parameter_name;
    }
    if (const auto parameter_unit = extract_json_string(text, "parameter_unit")) {
        config.parameter_unit = *parameter_unit;
    }
    if (const auto output_path = extract_json_string(text, "output")) {
        config.output_path = *output_path;
    }

    const auto entries_pos = text.find("\"entries\"");
    if (entries_pos == std::string::npos) {
        throw std::runtime_error("config file is missing entries array: " + config_path);
    }
    const auto array_begin = text.find('[', entries_pos);
    const auto array_end = text.find(']', array_begin);
    if (array_begin == std::string::npos || array_end == std::string::npos || array_end <= array_begin) {
        throw std::runtime_error("config file has malformed entries array: " + config_path);
    }
    const auto entries_text = text.substr(array_begin + 1, array_end - array_begin - 1);
    const std::regex object_pattern("\\{([^{}]*)\\}");
    for (auto iter = std::sregex_iterator(entries_text.begin(), entries_text.end(), object_pattern); iter != std::sregex_iterator(); ++iter) {
        const auto object_text = (*iter)[1].str();
        ConfigEntry entry;
        if (const auto raw_path = extract_json_string(object_text, "raw")) {
            entry.raw_path = *raw_path;
        }
        if (const auto histogram_path = extract_json_string(object_text, "histogram")) {
            entry.histogram_path = *histogram_path;
        }
        if (const auto value = extract_json_number(object_text, "value")) {
            entry.parameter_value = *value;
        } else {
            throw std::runtime_error("config entry is missing numeric value: " + config_path);
        }
        if (entry.histogram_path.empty()) {
            if (entry.raw_path.empty()) {
                throw std::runtime_error("config entry needs raw or histogram path: " + config_path);
            }
            entry.histogram_path = histogram_path_from_raw_path(entry.raw_path).string();
        }
        std::ostringstream label;
        if (!config.parameter_name.empty()) {
            label << config.parameter_name << " = ";
        }
        label << std::fixed << std::setprecision(std::abs(entry.parameter_value - std::round(entry.parameter_value)) < 1e-9 ? 0 : 1) << entry.parameter_value;
        if (!config.parameter_unit.empty()) {
            label << " " << config.parameter_unit;
        }
        if (!entry.raw_path.empty()) {
            label << " (" << std::filesystem::path(entry.raw_path).stem().string() << ")";
        }
        entry.label = label.str();
        config.entries.push_back(entry);
    }
    if (config.entries.empty()) {
        throw std::runtime_error("config file has no entries: " + config_path);
    }
    return config;
}

std::vector<std::string> discover_default_inputs()
{
    std::vector<std::string> inputs;
    const std::filesystem::path dump_005("dump/005");
    if (!std::filesystem::exists(dump_005)) {
        return inputs;
    }
    for (const auto& entry : std::filesystem::directory_iterator(dump_005)) {
        if (!entry.is_directory()) {
            continue;
        }
        const auto candidate = entry.path() / "event_tot_sum_histograms.root";
        if (std::filesystem::exists(candidate)) {
            inputs.push_back(candidate.string());
        }
    }
    std::sort(inputs.begin(), inputs.end());
    return inputs;
}

std::string join_arguments(int argc, char** argv)
{
    std::ostringstream stream;
    bool first_argument = true;
    for (int index = 0; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--config" || argument == "-c") {
            ++index;
            continue;
        }
        if (argument.rfind("--config=", 0) == 0) {
            continue;
        }
        if (!first_argument) {
            stream << ' ';
        }
        stream << argument;
        first_argument = false;
    }
    return stream.str();
}

std::vector<std::string> collect_input_paths_from_argv(int argc, char** argv)
{
    std::vector<std::string> inputs;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--input=", 0) == 0) {
            inputs.push_back(argument.substr(std::string("--input=").size()));
            continue;
        }
        if (argument != "--input" && argument != "-i") {
            continue;
        }
        for (++index; index < argc; ++index) {
            const std::string value = argv[index];
            if (!value.empty() && value.front() == '-') {
                --index;
                break;
            }
            inputs.push_back(value);
        }
    }
    return inputs;
}

std::string parameter_display_name(const std::string& parameter_name)
{
    return parameter_name == "temperature" ? "temp" : parameter_name;
}

std::vector<InputHistogram> load_histograms(const std::vector<std::string>& input_paths)
{
    std::vector<InputHistogram> inputs;
    for (const auto& input_path : input_paths) {
        TFile input(input_path.c_str(), "READ");
        if (input.IsZombie()) {
            throw std::runtime_error("failed to open input ROOT file: " + input_path);
        }
        auto* source_hist = dynamic_cast<TH1D*>(input.Get("h_event_tot_sum_all_lpgbt"));
        if (source_hist == nullptr) {
            throw std::runtime_error("input ROOT file does not contain TH1D 'h_event_tot_sum_all_lpgbt': " + input_path);
        }
        auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(source_hist->Clone((std::string("h_compare_") + std::to_string(inputs.size())).c_str())));
        if (!hist) {
            throw std::runtime_error("failed to clone TH1D from: " + input_path);
        }
        hist->SetDirectory(nullptr);
        inputs.push_back(InputHistogram{input_path, file_label(input_path), 0.0, false, std::move(hist)});
    }
    return inputs;
}

std::vector<InputHistogram> load_config_histograms(const ComparisonConfig& config)
{
    std::vector<InputHistogram> inputs;
    for (const auto& entry : config.entries) {
        TFile input(entry.histogram_path.c_str(), "READ");
        if (input.IsZombie()) {
            throw std::runtime_error("failed to open input ROOT file: " + entry.histogram_path);
        }
        auto* source_hist = dynamic_cast<TH1D*>(input.Get("h_event_tot_sum_all_lpgbt"));
        if (source_hist == nullptr) {
            throw std::runtime_error("input ROOT file does not contain TH1D 'h_event_tot_sum_all_lpgbt': " + entry.histogram_path);
        }
        auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(source_hist->Clone((std::string("h_compare_") + std::to_string(inputs.size())).c_str())));
        if (!hist) {
            throw std::runtime_error("failed to clone TH1D from: " + entry.histogram_path);
        }
        hist->SetDirectory(nullptr);
        inputs.push_back(InputHistogram{entry.histogram_path, entry.label, entry.parameter_value, true, std::move(hist)});
    }
    return inputs;
}

void draw_header(const std::string& title, const std::string& run_arguments, const std::string& run_time, double x = 0.12)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(13);
    latex.SetTextSize(0.042);
    latex.DrawLatex(x, 0.88, "#bf{FoCal TB2026 July}");
    latex.SetTextSize(0.032);
    latex.DrawLatex(x, 0.835, title.c_str());
    latex.DrawLatex(x, 0.795, run_arguments.c_str());
    latex.DrawLatex(x, 0.755, run_time.c_str());
}

void draw_comparison(std::vector<InputHistogram>& inputs, const std::string& output_path, const std::string& title,
    const std::string& parameter_name, const std::string& parameter_unit, const std::string& run_arguments, const std::string& run_time)
{
    if (inputs.empty()) {
        throw std::runtime_error("no input histograms available");
    }
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const std::vector<int> colors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2, kPink + 7};
    std::vector<std::unique_ptr<TH1D>> normalized_hists;
    double x_max = 1.0;
    double y_max = 0.0;
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        auto& hist = *inputs[index].hist;
        const auto color = colors[index % colors.size()];
        hist.SetLineColor(color);
        hist.SetLineWidth(2);
        hist.SetFillStyle(0);
        x_max = std::max(x_max, hist.GetXaxis()->GetXmax());

        auto normalized_hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(hist.Clone((std::string(hist.GetName()) + "_normalized").c_str())));
        if (!normalized_hist) {
            throw std::runtime_error("failed to clone histogram for normalized display: " + inputs[index].input_path);
        }
        normalized_hist->SetDirectory(nullptr);
        const auto integral = normalized_hist->Integral();
        const auto scale = integral > 0.0 ? 1.0 / integral : 1.0;
        normalized_hist->Scale(scale);
        normalized_hist->SetLineColor(color);
        normalized_hist->SetLineWidth(2);
        normalized_hist->SetFillStyle(0);
        y_max = std::max(y_max, normalized_hist->GetMaximum());
        normalized_hists.push_back(std::move(normalized_hist));
    }

    std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
    TCanvas canvas("tot_sum_comparison_canvas", "ToT sum comparison", 1300, 900);
    canvas.Print((output_path + "[").c_str());

    const auto draw_page = [&](bool log_y) {
        canvas.Clear();
        canvas.SetLogy(log_y ? 1 : 0);
        canvas.SetLeftMargin(0.10);
        canvas.SetRightMargin(0.08);
        canvas.SetBottomMargin(0.12);
        TH1D frame((log_y ? "h_tot_frame_log" : "h_tot_frame_linear"), ";#Sigma_{channels} first non-zero ToT;events / event", 120, 0.0, x_max);
        frame.SetMinimum(log_y ? std::max(1.0e-6, y_max * 1.0e-3) : 0.0);
        frame.SetMaximum(std::max(1.0e-6, y_max * (log_y ? 10.0 : 1.4)));
        frame.Draw("axis");
        TLegend legend(0.68, 0.70, 0.88, 0.86);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.SetTextSize(0.030);
        for (std::size_t index = 0; index < inputs.size(); ++index) {
            normalized_hists[index]->Draw("hist same");
            std::ostringstream label;
            if (inputs[index].has_parameter_value && !parameter_name.empty()) {
                label << parameter_display_name(parameter_name) << " = " << std::fixed
                    << std::setprecision(std::abs(inputs[index].parameter_value - std::round(inputs[index].parameter_value)) < 1e-9 ? 0 : 1)
                    << inputs[index].parameter_value;
                if (!parameter_unit.empty()) {
                    label << " " << parameter_unit;
                }
            } else {
                label << inputs[index].label;
            }
            legend.AddEntry(normalized_hists[index].get(), label.str().c_str(), "l");
        }
        legend.Draw();
        draw_header(title + (log_y ? " log-y" : " linear-y"), run_arguments, run_time);
        canvas.Print(output_path.c_str());
    };

    draw_page(false);
    draw_page(true);
    canvas.Print((output_path + "]").c_str());
}

} // namespace

int main(int argc, char** argv)
{
    const std::string script_name = "006_ToT_Sum_Comparison";
    try {
        spdlog::info("Running script: {}", script_name);
        spdlog::info("----------------------------------------");

        cxxopts::Options options(script_name, "Compare all-lpGBT ToT sum histograms from 005_ToT_Analysis ROOT outputs");
        options.add_options()
            ("c,config", "JSON comparison config with the same schema as 003", cxxopts::value<std::string>())
            ("i,input", "Input event_tot_sum_histograms.root file(s); repeat --input for multiple files", cxxopts::value<std::vector<std::string>>())
            ("o,output", "Output comparison PDF", cxxopts::value<std::string>()->default_value("dump/006/tot_sum_comparison.pdf"))
            ("h,help", "Print usage");

        const auto parsed = options.parse(argc, argv);
        if (parsed.count("help")) {
            std::cout << options.help() << '\n';
            return 0;
        }

        std::string output_path = parsed["output"].as<std::string>();
        std::string title = "All lpGBT ToT sum comparison";
        std::string parameter_name;
        std::string parameter_unit;
        std::vector<InputHistogram> inputs;
        if (parsed.count("config")) {
            const auto config = load_config(parsed["config"].as<std::string>());
            if (!config.output_path.empty() && !parsed.count("output")) {
                output_path = config.output_path;
            }
            title = config.title;
            parameter_name = config.parameter_name;
            parameter_unit = config.parameter_unit;
            inputs = load_config_histograms(config);
        }

        auto input_paths = collect_input_paths_from_argv(argc, argv);
        if (input_paths.empty() && parsed.count("input")) {
            input_paths = parsed["input"].as<std::vector<std::string>>();
        }
        if (inputs.empty() && input_paths.empty()) {
            input_paths = discover_default_inputs();
        }
        if (inputs.empty() && input_paths.empty()) {
            throw std::runtime_error("no input ROOT files provided and none found under dump/005/*/event_tot_sum_histograms.root");
        }
        if (inputs.empty()) {
            inputs = load_histograms(input_paths);
        }
        draw_comparison(inputs, output_path, title, parameter_name, parameter_unit, join_arguments(argc, argv), current_time_minute());

        spdlog::info("Read {} ToT sum histogram ROOT files", inputs.size());
        spdlog::info("Wrote ToT sum comparison PDF to {}", output_path);
        return 0;
    } catch (const std::exception& error) {
        spdlog::error("{}", error.what());
        return 1;
    }
}
