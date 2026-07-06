#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
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
	return std::filesystem::path("dump") / "002" / parent_name / "adc_val0_accumulated_waveform.pdf";
}

std::filesystem::path adc_sum_histogram_root_path(const std::string& output_path)
{
	return std::filesystem::path(output_path).parent_path() / "event_adc_sum_histograms.root";
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
	std::vector<double> adc_sums;
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
		auto hist = std::make_unique<TH2D>(("h_adc_val0_lpgbt_" + std::to_string(id) + "_ch_" + std::to_string(channel_index)).c_str(),
			("lpGBT " + std::to_string(id) + " row " + std::to_string(row_index) + " elink " + std::to_string(elink_index) + ";sample index;val0;count").c_str(),
			static_cast<int>(sample_bins), -0.5, static_cast<double>(sample_bins) - 0.5,
			1024, -0.5, 1023.5);
		histograms.push_back(std::move(hist));
	}
	return histograms;
}

AnalysisProducts collect_analysis_products(TTree& tree, const AnalysisLayout& layout)
{
	AnalysisProducts products;
	products.waveform_histogram = std::make_unique<TH2D>("h_adc_val0_accumulated_waveform",
		"Accumulated ADC waveform;sample index;val0;count",
		static_cast<int>(layout.sample_bins), -0.5, static_cast<double>(layout.sample_bins) - 0.5,
		1024, -0.5, 1023.5);
	for (const auto id : layout.active_lpgbts) {
		products.lpgbts[id].channel_histograms = make_lpgbt_channel_histograms(id, layout.sample_bins);
	}

	use_only_branches(tree, {sample_event_branch_name("val0"), sample_event_branch_name("lpgbt_id"), sample_event_branch_name("n_samples"), sample_event_branch_name("payload_offset"), sample_event_branch_name("payload_count")});
	std::vector<std::uint16_t>* val0 = nullptr;
	std::vector<std::uint8_t>* lpgbt_id = nullptr;
	std::vector<std::uint32_t>* n_samples = nullptr;
	std::vector<std::uint32_t>* payload_offset = nullptr;
	std::vector<std::uint32_t>* payload_count = nullptr;
	tree.SetBranchAddress(sample_event_branch_name("val0").c_str(), &val0);
	tree.SetBranchAddress(sample_event_branch_name("lpgbt_id").c_str(), &lpgbt_id);
	tree.SetBranchAddress(sample_event_branch_name("n_samples").c_str(), &n_samples);
	tree.SetBranchAddress(sample_event_branch_name("payload_offset").c_str(), &payload_offset);
	tree.SetBranchAddress(sample_event_branch_name("payload_count").c_str(), &payload_count);

	const auto entries = tree.GetEntries();
	for (Long64_t entry = 0; entry < entries; ++entry) {
		tree.GetEntry(entry);
		if (val0 == nullptr || lpgbt_id == nullptr || n_samples == nullptr || payload_offset == nullptr || payload_count == nullptr) {
			continue;
		}
		for (std::size_t lpgbt_index = 0; lpgbt_index < lpgbt_id->size() && lpgbt_index < n_samples->size() && lpgbt_index < payload_offset->size() && lpgbt_index < payload_count->size(); ++lpgbt_index) {
			const auto id = (*lpgbt_id)[lpgbt_index];
			const auto begin = static_cast<std::size_t>((*payload_offset)[lpgbt_index]);
			const auto count = static_cast<std::size_t>((*payload_count)[lpgbt_index]);
			const auto end = std::min(begin + count, val0->size());
			auto lpgbt_iter = products.lpgbts.find(id);
			for (std::size_t index = begin; index < end; ++index) {
				const auto local_index = index - begin;
				const auto sample_index = static_cast<double>(local_index / kValuesPerSample);
				const auto value = static_cast<double>((*val0)[index]);
				products.waveform_histogram->Fill(sample_index, value);
				if (lpgbt_iter != products.lpgbts.end()) {
					const auto channel_index = local_index % kValuesPerSample;
					if (channel_index < lpgbt_iter->second.channel_histograms.size()) {
						lpgbt_iter->second.channel_histograms[channel_index]->Fill(sample_index, value);
					}
				}
			}

			const auto sample_count = static_cast<std::size_t>((*n_samples)[lpgbt_index]);
			if (sample_count == 0 || count < kValuesPerSample || begin + count > val0->size() || lpgbt_iter == products.lpgbts.end()) {
				continue;
			}
			double adc_sum = 0.0;
			for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
				const auto pedestal = (*val0)[begin + channel_index];
				std::uint16_t peak = pedestal;
				for (std::size_t sample_index = 1; sample_index < sample_count; ++sample_index) {
					const auto value_index = begin + sample_index * kValuesPerSample + channel_index;
					if (value_index >= begin + count || value_index >= val0->size()) {
						break;
					}
					peak = std::max(peak, (*val0)[value_index]);
				}
				adc_sum += static_cast<double>(peak) - static_cast<double>(pedestal);
			}
			lpgbt_iter->second.adc_sums.push_back(adc_sum);
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
			TCanvas canvas(("adc_lpgbt_channel_waveform_canvas_" + std::to_string(id) + "_page_" + std::to_string(page_index)).c_str(),
				("ADC val0 channel waveforms for lpGBT " + std::to_string(id)).c_str(), 1600, 1200);
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

struct FitSummary {
	double mean = 0.0;
	double mean_stat = 0.0;
	double mean_sys = 0.0;
	double sigma = 0.0;
	double sigma_stat = 0.0;
	double sigma_sys = 0.0;
	double ratio = 0.0;
	double ratio_stat = 0.0;
	double ratio_sys = 0.0;
	std::unique_ptr<TF1> central_fit;
	std::vector<std::unique_ptr<TF1>> scan_fits;
};

double rms_spread(const std::vector<double>& values, double center)
{
	if (values.empty()) {
		return 0.0;
	}
	double squared_sum = 0.0;
	for (const auto value : values) {
		const auto delta = value - center;
		squared_sum += delta * delta;
	}
	return std::sqrt(squared_sum / static_cast<double>(values.size()));
}

double ratio_uncertainty(double mean, double sigma, double mean_error, double sigma_error)
{
	if (mean == 0.0) {
		return 0.0;
	}
	const auto d_ratio_d_sigma = 1.0 / mean;
	const auto d_ratio_d_mean = -sigma / (mean * mean);
	const auto variance = d_ratio_d_sigma * d_ratio_d_sigma * sigma_error * sigma_error
		+ d_ratio_d_mean * d_ratio_d_mean * mean_error * mean_error;
	return std::sqrt(std::max(0.0, variance));
}

std::optional<FitSummary> fit_adc_sum_histogram(TH1D& hist, double x_max)
{
	if (hist.GetEntries() < 3 || hist.GetRMS() <= 0.0) {
		return std::nullopt;
	}
	const auto stats_mean = hist.GetMean();
	const auto stats_sigma = hist.GetRMS();
	const auto prefit_min = std::max(0.0, stats_mean - 2.5 * stats_sigma);
	const auto prefit_max = std::min(x_max, stats_mean + 2.5 * stats_sigma);
	TF1 prefit((std::string("prefit_") + hist.GetName()).c_str(), "gaus", prefit_min, prefit_max);
	prefit.SetParameters(hist.GetMaximum(), stats_mean, stats_sigma);
	hist.Fit(&prefit, "QR0");

	const auto prefit_mean = prefit.GetParameter(1);
	const auto prefit_sigma = std::abs(prefit.GetParameter(2));
	if (prefit_sigma <= 0.0) {
		return std::nullopt;
	}
	const auto central_fit_min = std::max(0.0, prefit_mean - 2.5 * prefit_sigma);
	const auto central_fit_max = std::min(x_max, prefit_mean + 2.5 * prefit_sigma);
	auto central_fit = std::make_unique<TF1>((std::string("fit_") + hist.GetName()).c_str(), "gaus", central_fit_min, central_fit_max);
	central_fit->SetParameters(prefit.GetParameter(0), prefit_mean, prefit_sigma);
	central_fit->SetLineColor(kRed + 1);
	central_fit->SetLineStyle(2);
	central_fit->SetLineWidth(3);
	hist.Fit(central_fit.get(), "QR0");

	const auto mean = central_fit->GetParameter(1);
	const auto mean_stat = central_fit->GetParError(1);
	const auto sigma = std::abs(central_fit->GetParameter(2));
	const auto sigma_stat = central_fit->GetParError(2);
	const auto ratio = mean == 0.0 ? 0.0 : sigma / mean;
	const auto ratio_stat = ratio_uncertainty(mean, sigma, mean_stat, sigma_stat);

	std::vector<std::unique_ptr<TF1>> scan_fits;
	std::vector<double> scan_means;
	std::vector<double> scan_sigmas;
	std::vector<double> scan_ratios;
	for (double width_sigma = 2.0; width_sigma <= 3.5 + 1e-9; width_sigma += 0.5) {
		for (double center_shift_sigma = -0.5; center_shift_sigma <= 0.5 + 1e-9; center_shift_sigma += 0.1) {
			const auto scan_center = mean + center_shift_sigma * sigma;
			const auto fit_min = std::max(0.0, scan_center - width_sigma * sigma);
			const auto fit_max = std::min(x_max, scan_center + width_sigma * sigma);
			if (fit_max <= fit_min) {
				continue;
			}
			auto scan_fit = std::make_unique<TF1>((std::string("fit_scan_") + hist.GetName() + "_" + std::to_string(scan_fits.size())).c_str(), "gaus", fit_min, fit_max);
			scan_fit->SetParameters(central_fit->GetParameter(0), mean, sigma);
			scan_fit->SetLineColorAlpha(kRed + 1, 0.10);
			scan_fit->SetLineStyle(1);
			scan_fit->SetLineWidth(1);
			hist.Fit(scan_fit.get(), "QR0");
			const auto scan_mean = scan_fit->GetParameter(1);
			const auto scan_sigma = std::abs(scan_fit->GetParameter(2));
			if (scan_sigma <= 0.0 || scan_mean == 0.0) {
				continue;
			}
			scan_means.push_back(scan_mean);
			scan_sigmas.push_back(scan_sigma);
			scan_ratios.push_back(scan_sigma / scan_mean);
			scan_fits.push_back(std::move(scan_fit));
		}
	}

	const auto mean_sys = rms_spread(scan_means, mean);
	const auto sigma_sys = rms_spread(scan_sigmas, sigma);
	const auto ratio_sys = rms_spread(scan_ratios, ratio);
	return FitSummary{mean, mean_stat, mean_sys, sigma, sigma_stat, sigma_sys, ratio, ratio_stat, ratio_sys, std::move(central_fit), std::move(scan_fits)};
}

void draw_adc_sum_histogram_pages(const std::vector<double>& sums, const std::string& output_path, const std::string& name, const std::string& title, const PlotContext& context)
{
	if (sums.empty()) {
		return;
	}
	const auto max_sum = *std::max_element(sums.begin(), sums.end());
	const auto x_max = max_sum <= 0.0 ? 1.0 : max_sum * 1.05;
	TH1D hist(("h_" + name).c_str(), (title + ";#Sigma_{channels}(max(val0)-pedestal);events").c_str(), 120, 0.0, x_max);
	for (const auto sum : sums) {
		hist.Fill(sum);
	}
	hist.SetLineColor(kAzure + 1);
	hist.SetLineWidth(2);
	hist.SetFillStyle(0);
	const auto fit = fit_adc_sum_histogram(hist, x_max);

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
		if (fit.has_value()) {
			for (const auto& scan_fit : fit->scan_fits) {
				scan_fit->Draw("same");
			}
			fit->central_fit->Draw("same");
		}
		draw_page_header(context, page_title);

		TLatex latex;
		latex.SetNDC();
		latex.SetTextFont(42);
		latex.SetTextAlign(33);
		latex.SetTextSize(0.032);
		latex.DrawLatex(0.88, 0.88, page_title.c_str());
		latex.DrawLatex(0.88, 0.83, ("events = " + format_count(static_cast<std::size_t>(hist.GetEntries()))).c_str());
		if (fit.has_value()) {
			std::ostringstream mean_text;
			mean_text << "#mu = " << std::fixed << std::setprecision(2) << fit->mean
				<< " #pm " << fit->mean_stat << " (stat) #pm " << fit->mean_sys << " (sys)";
			latex.DrawLatex(0.88, 0.78, mean_text.str().c_str());
			std::ostringstream sigma_text;
			sigma_text << "#sigma = " << std::fixed << std::setprecision(2) << fit->sigma
				<< " #pm " << fit->sigma_stat << " (stat) #pm " << fit->sigma_sys << " (sys)";
			latex.DrawLatex(0.88, 0.73, sigma_text.str().c_str());
			std::ostringstream ratio_text;
			ratio_text << "#sigma/#mu = " << std::fixed << std::setprecision(4) << fit->ratio
				<< " #pm " << fit->ratio_stat << " (stat) #pm " << fit->ratio_sys << " (sys)";
			latex.DrawLatex(0.88, 0.68, ratio_text.str().c_str());
			std::ostringstream ratio_percent_text;
			ratio_percent_text << "resolution = " << std::fixed << std::setprecision(2) << fit->ratio * 100.0
				<< "% #pm " << fit->ratio_stat * 100.0 << "% (stat) #pm " << fit->ratio_sys * 100.0 << "% (sys)";
			latex.DrawLatex(0.88, 0.63, ratio_percent_text.str().c_str());
		}
		canvas.Print(output_path.c_str());
	};

	draw_page(true);
	draw_page(false);
}

std::unique_ptr<TH1D> make_adc_sum_histogram(const std::vector<double>& sums, const std::string& name, const std::string& title)
{
	const auto max_sum = sums.empty() ? 0.0 : *std::max_element(sums.begin(), sums.end());
	const auto x_max = max_sum <= 0.0 ? 1.0 : max_sum * 1.05;
	auto hist = std::make_unique<TH1D>(name.c_str(), (title + ";#Sigma_{channels}(max(val0)-pedestal);events").c_str(), 120, 0.0, x_max);
	for (const auto sum : sums) {
		hist->Fill(sum);
	}
	hist->SetDirectory(nullptr);
	return hist;
}

void write_all_lpgbt_adc_sum_histogram_root(const std::vector<double>& all_lpgbt_sums, const std::filesystem::path& root_path, const PlotContext& context)
{
	std::filesystem::create_directories(root_path.parent_path());
	TFile output(root_path.string().c_str(), "RECREATE");
	if (output.IsZombie()) {
		throw std::runtime_error("failed to create ADC sum histogram ROOT file: " + root_path.string());
	}
	TNamed("ADC_sum_input_file", context.input_path.c_str()).Write();
	TNamed("ADC_sum_histogram_note", "All lpGBT event ADC sum from 002_ADC_Analysis").Write();
	auto hist = make_adc_sum_histogram(all_lpgbt_sums, "h_event_adc_sum_all_lpgbt", "All lpGBT event ADC sum");
	hist->Write();
	output.Write();
	output.Close();
}

void draw_lpgbt_adc_sum_distribution_pages(const std::map<std::uint8_t, LpgbtAnalysis>& lpgbts, const std::string& output_path, const std::filesystem::path& adc_sum_root_path, const PlotContext& context)
{
	if (lpgbts.empty()) {
		write_all_lpgbt_adc_sum_histogram_root({}, adc_sum_root_path, context);
		return;
	}

	std::vector<double> all_lpgbt_sums;
	const auto event_count = std::max_element(lpgbts.begin(), lpgbts.end(), [](const auto& left, const auto& right) {
		return left.second.adc_sums.size() < right.second.adc_sums.size();
	})->second.adc_sums.size();
	all_lpgbt_sums.assign(event_count, 0.0);
	for (const auto& [id, analysis] : lpgbts) {
		(void)id;
		const auto& sums = analysis.adc_sums;
		for (std::size_t index = 0; index < sums.size() && index < all_lpgbt_sums.size(); ++index) {
			all_lpgbt_sums[index] += sums[index];
		}
	}
	for (const auto& [id, analysis] : lpgbts) {
		draw_adc_sum_histogram_pages(analysis.adc_sums, output_path, "event_adc_sum_lpgbt_" + std::to_string(id), "lpGBT " + std::to_string(id) + " event ADC sum", context);
	}
	write_all_lpgbt_adc_sum_histogram_root(all_lpgbt_sums, adc_sum_root_path, context);
	draw_adc_sum_histogram_pages(all_lpgbt_sums, output_path, "event_adc_sum_all_lpgbt", "All lpGBT event ADC sum", context);
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
	if (tree->GetBranch(sample_event_branch_name("val0").c_str()) == nullptr) {
		throw std::runtime_error("sample event tree does not contain branch '" + sample_event_branch_name("val0") + "': " + input_path);
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
	const auto adc_sum_root_path = adc_sum_histogram_root_path(output_path);
	TCanvas canvas("adc_waveform_canvas", "ADC val0 accumulated waveform", 1100, 800);
	canvas.SetRightMargin(0.16);
	canvas.SetLeftMargin(0.11);
	canvas.SetBottomMargin(0.11);
	canvas.SetLogz();
	products.waveform_histogram->GetZaxis()->SetTitle("count");
	products.waveform_histogram->Draw("COLZ");
	draw_page_header(context, "Accumulated ADC waveform");
	canvas.Print((output_path + "[").c_str());
	canvas.Print(output_path.c_str());
	draw_lpgbt_adc_sum_distribution_pages(products.lpgbts, output_path, adc_sum_root_path, context);
	draw_lpgbt_channel_waveform_pages(products.lpgbts, output_path);
	canvas.Print((output_path + "]").c_str());

	spdlog::info("Read {} entries from {}", tree->GetEntries(), input_path);
	spdlog::info("Wrote ADC waveform PDF to {}", output_path);
	spdlog::info("Wrote ADC sum histogram ROOT file to {}", adc_sum_root_path.string());
}

} // namespace

int main(int argc, char** argv)
{
	const std::string script_full_name = argv[0];
	const auto script_name = script_full_name.substr(script_full_name.find_last_of("/\\") + 1);
	spdlog::info("Running script: {}", script_name);
	spdlog::info("----------------------------------------");

	cxxopts::Options options(script_name, "Draw accumulated val0 waveform from 001_Rootifier sample_events.root FOCAL tree");
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
