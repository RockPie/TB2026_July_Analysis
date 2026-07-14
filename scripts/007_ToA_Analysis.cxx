#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TTree.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
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
constexpr double kSamplePeriodNs = 25.0;
constexpr double kToaLsbBinWidthNs = 0.025;
constexpr double kToaMsbBinWidthNs = 6.25;
constexpr double kDetailedWaveformBinWidthNs = 1.5;
constexpr std::array<std::uint16_t, 4> kToaShiftThresholdCandidates = {0u, 270u, 540u, 810u};

constexpr const char* kSampleEventTreeName = "FOCAL";
constexpr const char* kLegacySampleEventTreeName = "sample_events";
constexpr const char* kSampleEventBranchPrefix = "FOCALHCAL.";

std::string g_sample_event_branch_prefix = kSampleEventBranchPrefix;

struct PlotContext {
	std::string input_path;
	std::string run_arguments;
	std::string run_time;
};

struct AnalysisLayout {
	std::size_t sample_bins = 1;
	std::vector<std::uint8_t> active_lpgbts;
};

struct LpgbtAnalysis {
	std::vector<std::unique_ptr<TH2D>> channel_correlations;
	std::vector<std::unique_ptr<TH2D>> detailed_waveforms;
	std::vector<std::uint16_t> toa_shift_thresholds;
};

struct AnalysisProducts {
	std::map<std::uint8_t, LpgbtAnalysis> lpgbts;
};

struct ToaHit {
	std::size_t sample_index = 0;
	std::uint16_t raw_value = 0;
};

std::string file_name_only(const std::string& path)
{
	return std::filesystem::path(path).filename().string();
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

std::filesystem::path default_output_path(const std::string& input_path)
{
	const auto input_file = std::filesystem::path(input_path);
	auto parent_name = input_file.parent_path().filename();
	if (parent_name.empty()) {
		parent_name = input_file.stem();
	}
	return std::filesystem::path("dump") / "007" / parent_name / "toa_val2_adc_peak_correlation.pdf";
}

std::filesystem::path default_detailed_output_path(const std::string& output_path)
{
	const auto path = std::filesystem::path(output_path);
	return path.parent_path() / (path.stem().string() + "_detailed_waveform" + path.extension().string());
}

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

void draw_page_header(const PlotContext& context, const std::string& plot_title)
{
	TLatex latex;
	latex.SetNDC();
	latex.SetTextFont(42);
	latex.SetTextAlign(13);
	latex.SetTextSize(0.04);
	latex.DrawLatex(0.12, 0.94, "#bf{FoCal TB2026 July}");
	latex.SetTextSize(0.027);
	latex.DrawLatex(0.12, 0.905, plot_title.c_str());
	latex.DrawLatex(0.12, 0.875, file_name_only(context.input_path).c_str());
	latex.DrawLatex(0.12, 0.845, context.run_arguments.c_str());
	latex.DrawLatex(0.12, 0.815, context.run_time.c_str());
}

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

std::vector<std::unique_ptr<TH2D>> make_lpgbt_channel_correlations(std::uint8_t id, std::size_t sample_bins)
{
	std::vector<std::unique_ptr<TH2D>> histograms;
	histograms.reserve(kValuesPerSample);
	// const auto time_max = static_cast<double>(sample_bins) * kSamplePeriodNs + 16.0;
    const auto time_max = 200.0;
	for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
		const auto row_index = channel_index / kSampleElinkCount;
		const auto elink_index = channel_index % kSampleElinkCount;
		auto hist = std::make_unique<TH2D>(("h_toa_adc_peak_lpgbt_" + std::to_string(id) + "_ch_" + std::to_string(channel_index)).c_str(),
			("lpGBT " + std::to_string(id) + " row " + std::to_string(row_index) + " elink " + std::to_string(elink_index) + ";ToA time [ns];ADC peak = max(val0)-first(val0);events").c_str(),
			std::max(1, static_cast<int>(time_max / 0.5)), 0.0, time_max,
			512, -0.5, 1023.5);
		hist->SetDirectory(nullptr);
		histograms.push_back(std::move(hist));
	}
	return histograms;
}

std::vector<std::unique_ptr<TH2D>> make_lpgbt_detailed_waveforms(std::uint8_t id, std::size_t sample_bins)
{
	std::vector<std::unique_ptr<TH2D>> histograms;
	histograms.reserve(kValuesPerSample);
	const auto time_span = std::max(kSamplePeriodNs, static_cast<double>(sample_bins) * kSamplePeriodNs);
	const auto x_min = -100.0;
	const auto x_max = time_span;
	const auto x_bins = std::max(1, static_cast<int>((x_max - x_min) / kDetailedWaveformBinWidthNs));
	for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
		const auto row_index = channel_index / kSampleElinkCount;
		const auto elink_index = channel_index % kSampleElinkCount;
		auto hist = std::make_unique<TH2D>(("h_toa_aligned_adc_waveform_lpgbt_" + std::to_string(id) + "_ch_" + std::to_string(channel_index)).c_str(),
			("lpGBT " + std::to_string(id) + " row " + std::to_string(row_index) + " elink " + std::to_string(elink_index) + ";sample time - ToA time [ns];ADC - first sample ADC;events").c_str(),
			x_bins, x_min, x_max,
			512, -0.5, 1023.5);
		hist->SetDirectory(nullptr);
		histograms.push_back(std::move(hist));
	}
	return histograms;
}

double decoded_toa_time_ns(std::uint16_t raw_value)
{
	const auto lsb = raw_value & 0x00ffu;
	const auto msb = (raw_value >> 8u) & 0x0003u;
	return static_cast<double>(msb) * kToaMsbBinWidthNs + static_cast<double>(lsb) * kToaLsbBinWidthNs;
}

double corrected_toa_time_ns(std::size_t sample_index, std::uint16_t raw_value, std::uint16_t shift_threshold)
{
	auto time_ns = static_cast<double>(sample_index) * kSamplePeriodNs + decoded_toa_time_ns(raw_value);
	if (raw_value > shift_threshold) {
		time_ns -= kSamplePeriodNs;
	}
	return time_ns;
}

double toa_distribution_roughness(const std::vector<ToaHit>& hits, std::uint16_t shift_threshold, std::size_t sample_bins)
{
	if (hits.size() < 3) {
		return 0.0;
	}
	const auto x_min = -kSamplePeriodNs;
	const auto x_max = static_cast<double>(sample_bins) * kSamplePeriodNs + kSamplePeriodNs;
	constexpr double bin_width = 0.5;
	const auto bin_count = std::max<std::size_t>(3, static_cast<std::size_t>((x_max - x_min) / bin_width));
	std::vector<double> counts(bin_count, 0.0);
	for (const auto& hit : hits) {
		const auto time_ns = corrected_toa_time_ns(hit.sample_index, hit.raw_value, shift_threshold);
		if (time_ns < x_min || time_ns >= x_max) {
			continue;
		}
		const auto bin = static_cast<std::size_t>((time_ns - x_min) / bin_width);
		if (bin < counts.size()) {
			counts[bin] += 1.0;
		}
	}
	double roughness = 0.0;
	for (std::size_t index = 1; index + 1 < counts.size(); ++index) {
		const auto second_difference = counts[index - 1] - 2.0 * counts[index] + counts[index + 1];
		roughness += second_difference * second_difference;
	}
	return roughness / static_cast<double>(hits.size());
}

std::map<std::uint8_t, std::vector<std::uint16_t>> select_toa_shift_thresholds(TTree& tree, const AnalysisLayout& layout)
{
	std::map<std::uint8_t, std::vector<std::vector<ToaHit>>> hits_by_lpgbt;
	for (const auto id : layout.active_lpgbts) {
		hits_by_lpgbt[id].resize(kValuesPerSample);
	}

	use_only_branches(tree, {sample_event_branch_name("val2"), sample_event_branch_name("lpgbt_id"), sample_event_branch_name("n_samples"), sample_event_branch_name("payload_offset"), sample_event_branch_name("payload_count")});
	std::vector<std::uint16_t>* val2 = nullptr;
	std::vector<std::uint8_t>* lpgbt_id = nullptr;
	std::vector<std::uint32_t>* n_samples = nullptr;
	std::vector<std::uint32_t>* payload_offset = nullptr;
	std::vector<std::uint32_t>* payload_count = nullptr;
	tree.SetBranchAddress(sample_event_branch_name("val2").c_str(), &val2);
	tree.SetBranchAddress(sample_event_branch_name("lpgbt_id").c_str(), &lpgbt_id);
	tree.SetBranchAddress(sample_event_branch_name("n_samples").c_str(), &n_samples);
	tree.SetBranchAddress(sample_event_branch_name("payload_offset").c_str(), &payload_offset);
	tree.SetBranchAddress(sample_event_branch_name("payload_count").c_str(), &payload_count);

	const auto entries = tree.GetEntries();
	for (Long64_t entry = 0; entry < entries; ++entry) {
		tree.GetEntry(entry);
		if (val2 == nullptr || lpgbt_id == nullptr || n_samples == nullptr || payload_offset == nullptr || payload_count == nullptr) {
			continue;
		}
		for (std::size_t lpgbt_index = 0; lpgbt_index < lpgbt_id->size() && lpgbt_index < n_samples->size() && lpgbt_index < payload_offset->size() && lpgbt_index < payload_count->size(); ++lpgbt_index) {
			const auto id = (*lpgbt_id)[lpgbt_index];
			const auto begin = static_cast<std::size_t>((*payload_offset)[lpgbt_index]);
			const auto count = static_cast<std::size_t>((*payload_count)[lpgbt_index]);
			const auto sample_count = static_cast<std::size_t>((*n_samples)[lpgbt_index]);
			const auto end = std::min(begin + count, val2->size());
			auto lpgbt_iter = hits_by_lpgbt.find(id);
			if (sample_count == 0 || count < kValuesPerSample || begin + count > val2->size() || lpgbt_iter == hits_by_lpgbt.end()) {
				continue;
			}
			for (std::size_t channel_index = 0; channel_index < kValuesPerSample; ++channel_index) {
				for (std::size_t sample_index = 0; sample_index < sample_count; ++sample_index) {
					const auto value_index = begin + sample_index * kValuesPerSample + channel_index;
					if (value_index >= end) {
						break;
					}
					const auto raw_value = (*val2)[value_index];
					if (raw_value != 0) {
						lpgbt_iter->second[channel_index].push_back(ToaHit{sample_index, raw_value});
						break;
					}
				}
			}
		}
	}
	restore_tree_read_state(tree);

	std::map<std::uint8_t, std::vector<std::uint16_t>> thresholds_by_lpgbt;
	for (const auto& [id, channel_hits] : hits_by_lpgbt) {
		auto& thresholds = thresholds_by_lpgbt[id];
		thresholds.assign(kValuesPerSample, kToaShiftThresholdCandidates.front());
		for (std::size_t channel_index = 0; channel_index < channel_hits.size(); ++channel_index) {
			double best_roughness = std::numeric_limits<double>::infinity();
			for (const auto candidate : kToaShiftThresholdCandidates) {
				const auto roughness = toa_distribution_roughness(channel_hits[channel_index], candidate, layout.sample_bins);
				if (roughness < best_roughness) {
					best_roughness = roughness;
					thresholds[channel_index] = candidate;
				}
			}
		}
	}
	return thresholds_by_lpgbt;
}

AnalysisProducts collect_analysis_products(TTree& tree, const AnalysisLayout& layout)
{
	AnalysisProducts products;
	const auto selected_thresholds = select_toa_shift_thresholds(tree, layout);
	for (const auto id : layout.active_lpgbts) {
		products.lpgbts[id].channel_correlations = make_lpgbt_channel_correlations(id, layout.sample_bins);
		products.lpgbts[id].detailed_waveforms = make_lpgbt_detailed_waveforms(id, layout.sample_bins);
		auto selected_iter = selected_thresholds.find(id);
		products.lpgbts[id].toa_shift_thresholds = selected_iter == selected_thresholds.end()
			? std::vector<std::uint16_t>(kValuesPerSample, kToaShiftThresholdCandidates.front())
			: selected_iter->second;
	}

	use_only_branches(tree, {sample_event_branch_name("val0"), sample_event_branch_name("val2"), sample_event_branch_name("lpgbt_id"), sample_event_branch_name("n_samples"), sample_event_branch_name("payload_offset"), sample_event_branch_name("payload_count")});
	std::vector<std::uint16_t>* val0 = nullptr;
	std::vector<std::uint16_t>* val2 = nullptr;
	std::vector<std::uint8_t>* lpgbt_id = nullptr;
	std::vector<std::uint32_t>* n_samples = nullptr;
	std::vector<std::uint32_t>* payload_offset = nullptr;
	std::vector<std::uint32_t>* payload_count = nullptr;
	tree.SetBranchAddress(sample_event_branch_name("val0").c_str(), &val0);
	tree.SetBranchAddress(sample_event_branch_name("val2").c_str(), &val2);
	tree.SetBranchAddress(sample_event_branch_name("lpgbt_id").c_str(), &lpgbt_id);
	tree.SetBranchAddress(sample_event_branch_name("n_samples").c_str(), &n_samples);
	tree.SetBranchAddress(sample_event_branch_name("payload_offset").c_str(), &payload_offset);
	tree.SetBranchAddress(sample_event_branch_name("payload_count").c_str(), &payload_count);

	const auto entries = tree.GetEntries();
	for (Long64_t entry = 0; entry < entries; ++entry) {
		tree.GetEntry(entry);
		if (val0 == nullptr || val2 == nullptr || lpgbt_id == nullptr || n_samples == nullptr || payload_offset == nullptr || payload_count == nullptr) {
			continue;
		}
		for (std::size_t lpgbt_index = 0; lpgbt_index < lpgbt_id->size() && lpgbt_index < n_samples->size() && lpgbt_index < payload_offset->size() && lpgbt_index < payload_count->size(); ++lpgbt_index) {
			const auto id = (*lpgbt_id)[lpgbt_index];
			const auto begin = static_cast<std::size_t>((*payload_offset)[lpgbt_index]);
			const auto count = static_cast<std::size_t>((*payload_count)[lpgbt_index]);
			const auto sample_count = static_cast<std::size_t>((*n_samples)[lpgbt_index]);
			const auto end = std::min(begin + count, std::min(val0->size(), val2->size()));
			auto lpgbt_iter = products.lpgbts.find(id);
			if (sample_count == 0 || count < kValuesPerSample || begin + count > val0->size() || begin + count > val2->size() || lpgbt_iter == products.lpgbts.end()) {
				continue;
			}
			for (std::size_t channel_index = 0; channel_index < kValuesPerSample && channel_index < lpgbt_iter->second.channel_correlations.size() && channel_index < lpgbt_iter->second.detailed_waveforms.size() && channel_index < lpgbt_iter->second.toa_shift_thresholds.size(); ++channel_index) {
				const auto first_index = begin + channel_index;
				if (first_index >= end) {
					continue;
				}
				const auto baseline = static_cast<double>((*val0)[first_index]);
				double max_adc = baseline;
				bool found_toa = false;
				std::size_t toa_sample_index = 0;
				std::uint16_t toa_raw_value = 0;
				std::vector<double> adc_delta_by_sample;
				adc_delta_by_sample.reserve(sample_count);
				for (std::size_t sample_index = 0; sample_index < sample_count; ++sample_index) {
					const auto value_index = begin + sample_index * kValuesPerSample + channel_index;
					if (value_index >= end) {
						break;
					}
					const auto adc_value = static_cast<double>((*val0)[value_index]);
					max_adc = std::max(max_adc, adc_value);
					adc_delta_by_sample.push_back(adc_value - baseline);
					const auto toa_raw = (*val2)[value_index];
					if (!found_toa && toa_raw != 0) {
						found_toa = true;
						toa_sample_index = sample_index;
						toa_raw_value = toa_raw;
					}
				}
				if (found_toa) {
					const auto toa_time_ns = corrected_toa_time_ns(toa_sample_index, toa_raw_value, lpgbt_iter->second.toa_shift_thresholds[channel_index]);
					lpgbt_iter->second.channel_correlations[channel_index]->Fill(toa_time_ns, max_adc - baseline);
					for (std::size_t sample_index = 0; sample_index < adc_delta_by_sample.size(); ++sample_index) {
						const auto relative_time_ns = static_cast<double>(sample_index) * kSamplePeriodNs - toa_time_ns;
						lpgbt_iter->second.detailed_waveforms[channel_index]->Fill(relative_time_ns, adc_delta_by_sample[sample_index]);
					}
				}
			}
		}
	}
	restore_tree_read_state(tree);
	return products;
}

void draw_lpgbt_channel_correlation_pages(const std::map<std::uint8_t, LpgbtAnalysis>& lpgbts, const std::string& output_path, const PlotContext& context)
{
	constexpr int columns = 4;
	constexpr int rows = 4;
	constexpr std::size_t channels_per_page = columns * rows;
	for (const auto& [id, analysis] : lpgbts) {
		const auto& histograms = analysis.channel_correlations;
		for (std::size_t page_start = 0; page_start < histograms.size(); page_start += channels_per_page) {
			const auto page_index = page_start / channels_per_page;
			TCanvas canvas(("toa_adc_peak_lpgbt_" + std::to_string(id) + "_page_" + std::to_string(page_index)).c_str(),
				("ToA vs ADC peak lpGBT " + std::to_string(id)).c_str(), 1800, 1300);
			canvas.Divide(columns, rows, 0.001, 0.001);
			const auto page_end = std::min(page_start + channels_per_page, histograms.size());
			for (std::size_t channel_index = page_start; channel_index < page_end; ++channel_index) {
				canvas.cd(static_cast<int>(channel_index - page_start + 1));
				gPad->SetRightMargin(0.17);
				gPad->SetLeftMargin(0.12);
				gPad->SetBottomMargin(0.12);
				gPad->SetLogz();
				histograms[channel_index]->GetXaxis()->SetLabelSize(0.045);
				histograms[channel_index]->GetYaxis()->SetLabelSize(0.045);
				histograms[channel_index]->GetZaxis()->SetLabelSize(0.045);
				histograms[channel_index]->GetXaxis()->SetTitleSize(0.05);
				histograms[channel_index]->GetYaxis()->SetTitleSize(0.05);
				histograms[channel_index]->Draw("COLZ");
			}
			(void)context;
			canvas.Print(output_path.c_str());
		}
	}
}

void draw_lpgbt_detailed_waveform_pages(const std::map<std::uint8_t, LpgbtAnalysis>& lpgbts, const std::string& output_path)
{
	constexpr int columns = 4;
	constexpr int rows = 4;
	constexpr std::size_t channels_per_page = columns * rows;
	for (const auto& [id, analysis] : lpgbts) {
		const auto& histograms = analysis.detailed_waveforms;
		for (std::size_t page_start = 0; page_start < histograms.size(); page_start += channels_per_page) {
			const auto page_index = page_start / channels_per_page;
			TCanvas canvas(("toa_aligned_adc_waveform_lpgbt_" + std::to_string(id) + "_page_" + std::to_string(page_index)).c_str(),
				("ToA-aligned ADC waveform lpGBT " + std::to_string(id)).c_str(), 1800, 1300);
			canvas.Divide(columns, rows, 0.001, 0.001);
			const auto page_end = std::min(page_start + channels_per_page, histograms.size());
			for (std::size_t channel_index = page_start; channel_index < page_end; ++channel_index) {
				canvas.cd(static_cast<int>(channel_index - page_start + 1));
				gPad->SetRightMargin(0.17);
				gPad->SetLeftMargin(0.12);
				gPad->SetBottomMargin(0.12);
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

void draw_analysis(const std::string& input_path, const std::string& output_path, const std::string& detailed_output_path, const PlotContext& context)
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
	for (const auto* branch : {"val0", "val2", "lpgbt_id", "n_samples", "payload_offset", "payload_count"}) {
		if (tree->GetBranch(sample_event_branch_name(branch).c_str()) == nullptr) {
			throw std::runtime_error("sample event tree does not contain branch '" + sample_event_branch_name(branch) + "': " + input_path);
		}
	}

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kViridis);

	const auto layout = scan_analysis_layout(*tree);
	auto products = collect_analysis_products(*tree, layout);

	std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
	TCanvas open_canvas("toa_adc_peak_open", "ToA ADC peak correlation", 1000, 800);
	open_canvas.Print((output_path + "[").c_str());
	draw_lpgbt_channel_correlation_pages(products.lpgbts, output_path, context);
	open_canvas.Print((output_path + "]").c_str());

	TCanvas detailed_open_canvas("toa_aligned_adc_waveform_open", "ToA-aligned ADC waveform", 1000, 800);
	detailed_open_canvas.Print((detailed_output_path + "[").c_str());
	draw_lpgbt_detailed_waveform_pages(products.lpgbts, detailed_output_path);
	detailed_open_canvas.Print((detailed_output_path + "]").c_str());

	spdlog::info("Read {} entries from {}", tree->GetEntries(), input_path);
	spdlog::info("Wrote ToA ADC peak correlation PDF to {}", output_path);
	spdlog::info("Wrote ToA-aligned detailed waveform PDF to {}", detailed_output_path);
}

} // namespace

int main(int argc, char** argv)
{
	const std::string script_full_name = argv[0];
	const auto script_name = script_full_name.substr(script_full_name.find_last_of("/\\") + 1);
	spdlog::info("Running script: {}", script_name);
	spdlog::info("----------------------------------------");

	cxxopts::Options options(script_name, "Draw val2 ToA time vs ADC peak correlations from 001_Rootifier sample_events.root FOCAL tree");
	options.add_options()
		("i,input", "Input sample_events.root file containing FOCAL/FOCALHCAL branches", cxxopts::value<std::string>())
		("o,output", "Output PDF path", cxxopts::value<std::string>())
		("d,detailed-output", "Output PDF path for ToA-aligned detailed waveforms", cxxopts::value<std::string>())
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
		const auto detailed_output_path = parsed_opts.count("detailed-output")
			? parsed_opts["detailed-output"].as<std::string>()
			: default_detailed_output_path(output_path).string();
		const PlotContext context{input_path, command_line(argc, argv), format_run_time()};
		draw_analysis(input_path, output_path, detailed_output_path, context);
	} catch (const std::exception& error) {
		spdlog::error("{}", error.what());
		return 1;
	}

	return 0;
}
