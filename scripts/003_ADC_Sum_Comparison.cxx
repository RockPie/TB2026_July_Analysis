#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TBox.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <limits>
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
	double parameter_error = 0.0;
	bool has_parameter_value = false;
	std::unique_ptr<TH1D> hist;
};

struct ConfigEntry {
	std::string raw_path;
	std::string histogram_path;
	std::string label;
	double parameter_value = 0.0;
	double parameter_error = 0.0;
	double parameter_error_percent = 0.0;
};

struct ComparisonConfig {
	std::string title = "All lpGBT ADC sum comparison";
	std::string parameter_name;
	std::string parameter_unit;
	double parameter_error = 0.0;
	double parameter_error_percent = 0.0;
	bool mu_linear_fit = false;
	bool compare_to_other_tb = false;
	bool preliminary = false;
	std::string output_path;
	std::vector<ConfigEntry> entries;
};

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

bool extract_json_bool_or_number(const std::string& text, const std::string& key, bool default_value = false)
{
	if (const auto value = extract_json_number(text, key)) {
		return *value != 0.0;
	}
	const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*(true|false)");
	std::smatch match;
	if (!std::regex_search(text, match, pattern)) {
		return default_value;
	}
	return match[1].str() == "true";
}

double percent_error(double value, double percent)
{
	return std::abs(value) * percent / 100.0;
}

std::filesystem::path histogram_path_from_raw_path(const std::string& raw_path)
{
	const auto raw_file = std::filesystem::path(raw_path);
	return std::filesystem::path("dump") / "002" / raw_file.filename() / "event_adc_sum_histograms.root";
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
	if (const auto parameter_error = extract_json_number(text, "parameter_error")) {
		config.parameter_error = *parameter_error;
	}
	if (const auto parameter_error_percent = extract_json_number(text, "parameter_error_percent")) {
		config.parameter_error_percent = *parameter_error_percent;
	}
	config.mu_linear_fit = extract_json_bool_or_number(text, "mu_linear_fit");
	config.compare_to_other_tb = extract_json_bool_or_number(text, "compare_to_other_tb");
	config.preliminary = extract_json_bool_or_number(text, "preliminary");
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
		entry.parameter_error_percent = config.parameter_error_percent;
		entry.parameter_error = config.parameter_error_percent != 0.0
			? percent_error(entry.parameter_value, config.parameter_error_percent)
			: config.parameter_error;
		if (const auto error_percent = extract_json_number(object_text, "error_percent")) {
			entry.parameter_error_percent = *error_percent;
			entry.parameter_error = percent_error(entry.parameter_value, *error_percent);
		}
		if (const auto error = extract_json_number(object_text, "error")) {
			entry.parameter_error = *error;
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
	const std::filesystem::path dump_002("dump/002");
	if (!std::filesystem::exists(dump_002)) {
		return inputs;
	}
	for (const auto& entry : std::filesystem::directory_iterator(dump_002)) {
		if (!entry.is_directory()) {
			continue;
		}
		const auto candidate = entry.path() / "event_adc_sum_histograms.root";
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

std::optional<FitSummary> fit_adc_sum_histogram(TH1D& hist, int color)
{
	if (hist.GetEntries() < 3 || hist.GetRMS() <= 0.0) {
		return std::nullopt;
	}
	const auto x_max = hist.GetXaxis()->GetXmax();
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
	central_fit->SetLineColor(color);
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
			scan_fit->SetLineColorAlpha(color, 0.10);
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

std::vector<InputHistogram> load_histograms(const std::vector<std::string>& input_paths)
{
	std::vector<InputHistogram> inputs;
	for (const auto& input_path : input_paths) {
		TFile input(input_path.c_str(), "READ");
		if (input.IsZombie()) {
			throw std::runtime_error("failed to open input ROOT file: " + input_path);
		}
		auto* source_hist = dynamic_cast<TH1D*>(input.Get("h_event_adc_sum_all_lpgbt"));
		if (source_hist == nullptr) {
			throw std::runtime_error("input ROOT file does not contain TH1D 'h_event_adc_sum_all_lpgbt': " + input_path);
		}
		auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(source_hist->Clone((std::string("h_compare_") + std::to_string(inputs.size())).c_str())));
		if (!hist) {
			throw std::runtime_error("failed to clone TH1D from: " + input_path);
		}
		hist->SetDirectory(nullptr);
		inputs.push_back(InputHistogram{input_path, file_label(input_path), 0.0, 0.0, false, std::move(hist)});
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
		auto* source_hist = dynamic_cast<TH1D*>(input.Get("h_event_adc_sum_all_lpgbt"));
		if (source_hist == nullptr) {
			throw std::runtime_error("input ROOT file does not contain TH1D 'h_event_adc_sum_all_lpgbt': " + entry.histogram_path);
		}
		auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(source_hist->Clone((std::string("h_compare_") + std::to_string(inputs.size())).c_str())));
		if (!hist) {
			throw std::runtime_error("failed to clone TH1D from: " + entry.histogram_path);
		}
		hist->SetDirectory(nullptr);
		inputs.push_back(InputHistogram{entry.histogram_path, entry.label, entry.parameter_value, entry.parameter_error, true, std::move(hist)});
	}
	return inputs;
}

bool all_inputs_have_parameter_values(const std::vector<InputHistogram>& inputs)
{
	return std::all_of(inputs.begin(), inputs.end(), [](const auto& input) {
		return input.has_parameter_value;
	});
}

std::string parameter_display_name(const std::string& parameter_name)
{
	return parameter_name == "temperature" ? "temp" : parameter_name;
}

void draw_preliminary(bool preliminary, double x = 0.12, double y = 0.715)
{
	if (!preliminary) {
		return;
	}
	TLatex latex;
	latex.SetNDC();
	latex.SetTextFont(42);
	latex.SetTextAlign(13);
	latex.SetTextSize(0.032);
	latex.DrawLatex(x, y, "#it{Preliminary}");
}

void draw_header(const std::string& title, const std::string& run_arguments, const std::string& run_time, bool preliminary, double x = 0.12)
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
	draw_preliminary(preliminary, x);
}

void draw_fit_trend_page(TCanvas& canvas, const std::vector<InputHistogram>& inputs, const std::vector<std::optional<FitSummary>>& fits,
	const std::string& output_path, const std::string& y_title, const std::string& page_title, const std::string& parameter_name,
	const std::string& parameter_unit, const std::string& run_arguments, const std::string& run_time, bool preliminary,
	double (*value_getter)(const FitSummary&), double (*stat_getter)(const FitSummary&), double (*sys_getter)(const FitSummary&),
	std::optional<std::pair<double, double>> y_range = std::nullopt, bool linear_fit = false)
{
	std::vector<std::size_t> indices;
	for (std::size_t index = 0; index < inputs.size() && index < fits.size(); ++index) {
		if (fits[index].has_value()) {
			indices.push_back(index);
		}
	}
	if (indices.empty()) {
		return;
	}
	std::sort(indices.begin(), indices.end(), [&](const auto left, const auto right) {
		return inputs[left].parameter_value < inputs[right].parameter_value;
	});

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> x_error;
	std::vector<double> stat;
	std::vector<double> sys;
	std::vector<double> total;
	for (const auto index : indices) {
		x.push_back(inputs[index].parameter_value);
		x_error.push_back(inputs[index].parameter_error);
		y.push_back(value_getter(*fits[index]));
		stat.push_back(stat_getter(*fits[index]));
		sys.push_back(sys_getter(*fits[index]));
		total.push_back(std::hypot(stat.back(), sys.back()));
	}
	std::vector<double> zero(x.size(), 0.0);
	if (!y_range.has_value()) {
		double y_min = std::numeric_limits<double>::max();
		double y_max = std::numeric_limits<double>::lowest();
		for (std::size_t index = 0; index < y.size(); ++index) {
			const auto error = total[index];
			y_min = std::min(y_min, y[index] - error);
			y_max = std::max(y_max, y[index] + error);
		}
		if (y_min < y_max) {
			const auto lower_padding = std::max(1.0e-9, (y_max - y_min) * 0.12);
			const auto upper_padding = std::max(1.0e-9, (y_max - y_min) * 0.35);
			y_range = std::pair<double, double>{y_min - lower_padding, y_max + upper_padding};
		}
	}

	canvas.Clear();
	canvas.SetLogy(0);
	canvas.SetLeftMargin(0.12);
	canvas.SetRightMargin(0.08);
	canvas.SetBottomMargin(0.13);

	TGraphErrors center_graph(static_cast<int>(x.size()), x.data(), y.data(), zero.data(), zero.data());
	const auto x_axis_name = parameter_display_name(parameter_name);
	center_graph.SetTitle((";" + x_axis_name + (parameter_unit.empty() ? "" : " [" + parameter_unit + "]") + ";" + y_title).c_str());
	center_graph.SetMarkerStyle(20);
	center_graph.SetMarkerSize(1.1);
	center_graph.SetMarkerColor(kBlue + 1);
	center_graph.SetLineColor(kBlue + 1);
	center_graph.SetLineWidth(2);
	center_graph.Draw("APL");
	if (y_range.has_value()) {
		center_graph.GetYaxis()->SetRangeUser(y_range->first, y_range->second);
	}
	center_graph.GetXaxis()->SetTitleSize(0.045);
	center_graph.GetYaxis()->SetTitleSize(0.045);
	center_graph.GetXaxis()->SetLabelSize(0.04);
	center_graph.GetYaxis()->SetLabelSize(0.04);

	std::vector<std::unique_ptr<TBox>> boxes;
	const auto x_minmax = std::minmax_element(x.begin(), x.end());
	const auto x_span = std::max(1.0, *x_minmax.second - *x_minmax.first);
	const auto fallback_box_half_width = x_span * 0.025;
	for (std::size_t index = 0; index < x.size(); ++index) {
		const auto box_half_width = x_error[index] > 0.0 ? x_error[index] : fallback_box_half_width;
		auto total_box = std::make_unique<TBox>(x[index] - box_half_width, y[index] - total[index], x[index] + box_half_width, y[index] + total[index]);
		total_box->SetFillColorAlpha(kRed + 1, 0.22);
		total_box->SetLineColor(kRed + 1);
		total_box->SetLineWidth(1);
		total_box->Draw("l same");
		boxes.push_back(std::move(total_box));

		auto stat_box = std::make_unique<TBox>(x[index] - box_half_width, y[index] - stat[index], x[index] + box_half_width, y[index] + stat[index]);
		stat_box->SetFillColorAlpha(kBlue + 1, 0.40);
		stat_box->SetLineColor(kBlue + 1);
		stat_box->SetLineWidth(1);
		stat_box->Draw("l same");
		boxes.push_back(std::move(stat_box));
	}
	std::unique_ptr<TF1> trend_fit;
	if (linear_fit && x.size() >= 2) {
		TGraphErrors fit_graph(static_cast<int>(x.size()), x.data(), y.data(), x_error.data(), total.data());
		trend_fit = std::make_unique<TF1>("f_mu_linear_fit", "pol1", *x_minmax.first, *x_minmax.second);
		trend_fit->SetLineColor(kGreen + 2);
		trend_fit->SetLineStyle(2);
		trend_fit->SetLineWidth(3);
		fit_graph.Fit(trend_fit.get(), "Q0");
		trend_fit->Draw("same");

		TLatex latex;
		latex.SetNDC();
		latex.SetTextFont(42);
		latex.SetTextAlign(33);
		latex.SetTextSize(0.030);
		latex.DrawLatex(0.88, 0.68, "#mu linear fit");
		std::ostringstream intercept_text;
		intercept_text << "p_{0} = " << std::fixed << std::setprecision(2) << trend_fit->GetParameter(0)
			<< " #pm " << trend_fit->GetParError(0);
		latex.DrawLatex(0.88, 0.64, intercept_text.str().c_str());
		std::ostringstream slope_text;
		slope_text << "p_{1} = " << std::scientific << std::setprecision(3) << trend_fit->GetParameter(1)
			<< " #pm " << trend_fit->GetParError(1);
		latex.DrawLatex(0.88, 0.60, slope_text.str().c_str());
		std::ostringstream chi2_text;
		chi2_text << "#chi^{2}/ndf = " << std::fixed << std::setprecision(2) << trend_fit->GetChisquare()
			<< "/" << trend_fit->GetNDF();
		latex.DrawLatex(0.88, 0.56, chi2_text.str().c_str());
	}
	center_graph.Draw("PL same");

	TLegend legend(0.70, 0.72, 0.90, 0.84);
	legend.SetBorderSize(0);
	legend.SetFillStyle(0);
	legend.SetTextSize(0.030);
	if (boxes.size() >= 2) {
		legend.AddEntry(boxes[0].get(), "stat #oplus sys", "f");
		legend.AddEntry(boxes[1].get(), "stat", "f");
	}
	legend.AddEntry(&center_graph, "central", "lp");
	if (trend_fit) {
		legend.AddEntry(trend_fit.get(), "linear fit", "l");
	}
	legend.Draw();
	draw_header(page_title, run_arguments, run_time, preliminary, 0.15);
	canvas.Print(output_path.c_str());
}

double fit_mean(const FitSummary& fit) { return fit.mean; }
double fit_mean_stat(const FitSummary& fit) { return fit.mean_stat; }
double fit_mean_sys(const FitSummary& fit) { return fit.mean_sys; }
double fit_sigma(const FitSummary& fit) { return fit.sigma; }
double fit_sigma_stat(const FitSummary& fit) { return fit.sigma_stat; }
double fit_sigma_sys(const FitSummary& fit) { return fit.sigma_sys; }
double fit_resolution_percent(const FitSummary& fit) { return fit.ratio * 100.0; }
double fit_resolution_percent_stat(const FitSummary& fit) { return fit.ratio_stat * 100.0; }
double fit_resolution_percent_sys(const FitSummary& fit) { return fit.ratio_sys * 100.0; }

void draw_resolution_reference_comparison_page(const std::vector<InputHistogram>& inputs, const std::vector<std::optional<FitSummary>>& fits,
	const std::string& output_path, const std::string& run_arguments, const std::string& run_time, bool preliminary)
{
	std::vector<std::size_t> indices;
	for (std::size_t index = 0; index < inputs.size() && index < fits.size(); ++index) {
		if (fits[index].has_value()) {
			indices.push_back(index);
		}
	}
	if (indices.empty()) {
		return;
	}
	std::sort(indices.begin(), indices.end(), [&](const auto left, const auto right) {
		return inputs[left].parameter_value < inputs[right].parameter_value;
	});

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> x_error;
	std::vector<double> y_error;
	for (const auto index : indices) {
		x.push_back(inputs[index].parameter_value);
		x_error.push_back(inputs[index].parameter_error);
		y.push_back(fit_resolution_percent(*fits[index]));
		y_error.push_back(std::hypot(fit_resolution_percent_stat(*fits[index]), fit_resolution_percent_sys(*fits[index])));
	}

	TFile reference_file("reference/scan_0_analysis.root", "READ");
	if (reference_file.IsZombie()) {
		throw std::runtime_error("failed to open reference ROOT file: reference/scan_0_analysis.root");
	}
	auto* reference_canvas = dynamic_cast<TCanvas*>(reference_file.Get("canvas_resolution"));
	if (reference_canvas == nullptr) {
		throw std::runtime_error("reference ROOT file does not contain canvas_resolution: reference/scan_0_analysis.root");
	}
	auto comparison_canvas = std::unique_ptr<TCanvas>(dynamic_cast<TCanvas*>(reference_canvas->Clone("canvas_resolution_tb2026_comparison")));
	if (!comparison_canvas) {
		throw std::runtime_error("failed to clone reference canvas_resolution");
	}
	comparison_canvas->SetTitle("Resolution comparison with other test beams");
	comparison_canvas->cd();
	std::vector<TObject*> inherited_labels;
	for (auto* object : *comparison_canvas->GetListOfPrimitives()) {
		if (dynamic_cast<TLatex*>(object) != nullptr) {
			inherited_labels.push_back(object);
		}
	}
	for (auto* label : inherited_labels) {
		comparison_canvas->GetListOfPrimitives()->Remove(label);
		delete label;
	}

	auto* tb2026_graph = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), x_error.data(), y_error.data());
	tb2026_graph->SetName("tb2026_resolution_graph");
	tb2026_graph->SetTitle("TB2026 July");
	tb2026_graph->SetMarkerStyle(20);
	tb2026_graph->SetMarkerSize(1.2);
	tb2026_graph->SetMarkerColor(kAzure + 2);
	tb2026_graph->SetLineColor(kAzure + 2);
	tb2026_graph->SetLineWidth(3);
	tb2026_graph->Draw("P same");

	if (auto* legend = dynamic_cast<TLegend*>(comparison_canvas->GetListOfPrimitives()->FindObject("TPave"))) {
		legend->AddEntry(tb2026_graph, "FoCal-H TB2026 July", "p");
		legend->Draw();
	} else {
		auto* fallback_legend = new TLegend(0.58, 0.68, 0.88, 0.86);
		fallback_legend->SetBorderSize(0);
		fallback_legend->SetFillStyle(0);
		fallback_legend->SetTextSize(0.030);
		fallback_legend->AddEntry(tb2026_graph, "FoCal-H TB2026 July", "p");
		fallback_legend->Draw();
		comparison_canvas->Update();
	}
	draw_header("Resolution comparison with other test beams", run_arguments, run_time, preliminary);
	comparison_canvas->Modified();
	comparison_canvas->Update();
	comparison_canvas->Print(output_path.c_str());
	comparison_canvas.release();
}

void draw_comparison(std::vector<InputHistogram>& inputs, const std::string& output_path, const std::string& title,
	const std::string& parameter_name, const std::string& parameter_unit, const std::string& run_arguments, const std::string& run_time,
	bool mu_linear_fit, bool compare_to_other_tb, bool preliminary)
{
	if (inputs.empty()) {
		throw std::runtime_error("no input histograms available");
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	const std::vector<int> colors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2, kPink + 7};
	std::vector<std::optional<FitSummary>> fits;
	std::vector<std::unique_ptr<TH1D>> normalized_hists;
	std::vector<double> normalization_scales;
	double x_max = 1.0;
	double y_max = 0.0;
	for (std::size_t index = 0; index < inputs.size(); ++index) {
		auto& hist = *inputs[index].hist;
		const auto color = colors[index % colors.size()];
		hist.SetLineColor(color);
		hist.SetLineWidth(2);
		hist.SetFillStyle(0);
		x_max = std::max(x_max, hist.GetXaxis()->GetXmax());
		fits.push_back(fit_adc_sum_histogram(hist, color));

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
		normalization_scales.push_back(scale);
		normalized_hists.push_back(std::move(normalized_hist));
	}

	std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
	TCanvas canvas("adc_sum_comparison_canvas", "ADC sum comparison", 1300, 900);
	canvas.Print((output_path + "[").c_str());

	const auto draw_page = [&](bool log_y) {
		canvas.Clear();
		canvas.SetLogy(log_y ? 1 : 0);
		canvas.SetLeftMargin(0.10);
		canvas.SetRightMargin(0.08);
		canvas.SetBottomMargin(0.12);
		TH1D frame((log_y ? "h_frame_log" : "h_frame_linear"), ";#Sigma_{channels}(max(val0)-pedestal);events / event", 120, 0.0, x_max);
		frame.SetMinimum(log_y ? std::max(1.0e-6, y_max * 1.0e-3) : 0.0);
		frame.SetMaximum(std::max(1.0e-6, y_max * (log_y ? 10.0 : 1.4)));
		frame.Draw("axis");
		TLegend legend(0.68, 0.70, 0.88, 0.86);
		legend.SetBorderSize(0);
		legend.SetFillStyle(0);
		legend.SetTextSize(0.030);
		std::vector<std::unique_ptr<TF1>> display_fits;
		for (std::size_t index = 0; index < inputs.size(); ++index) {
			normalized_hists[index]->Draw("hist same");
			if (fits[index].has_value()) {
				const auto scale = normalization_scales[index];
				for (const auto& scan_fit : fits[index]->scan_fits) {
					auto display_fit = std::unique_ptr<TF1>(dynamic_cast<TF1*>(scan_fit->Clone((std::string(scan_fit->GetName()) + "_display").c_str())));
					if (display_fit) {
						display_fit->SetParameter(0, display_fit->GetParameter(0) * scale);
						display_fit->Draw("same");
						display_fits.push_back(std::move(display_fit));
					}
				}
				auto display_fit = std::unique_ptr<TF1>(dynamic_cast<TF1*>(fits[index]->central_fit->Clone((std::string(fits[index]->central_fit->GetName()) + "_display").c_str())));
				if (display_fit) {
					display_fit->SetParameter(0, display_fit->GetParameter(0) * scale);
					display_fit->Draw("same");
					display_fits.push_back(std::move(display_fit));
				}
			}
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
		draw_header(title + (log_y ? " log-y" : " linear-y"), run_arguments, run_time, preliminary);
		canvas.Print(output_path.c_str());
	};

	draw_page(false);
	draw_page(true);
	if (all_inputs_have_parameter_values(inputs) && !parameter_name.empty()) {
		draw_fit_trend_page(canvas, inputs, fits, output_path, "#mu", title + " #mu trend", parameter_name, parameter_unit, run_arguments, run_time, preliminary, fit_mean, fit_mean_stat, fit_mean_sys, std::nullopt, mu_linear_fit);
		draw_fit_trend_page(canvas, inputs, fits, output_path, "#sigma", title + " #sigma trend", parameter_name, parameter_unit, run_arguments, run_time, preliminary, fit_sigma, fit_sigma_stat, fit_sigma_sys);
		draw_fit_trend_page(canvas, inputs, fits, output_path, "resolution [%]", title + " resolution trend", parameter_name, parameter_unit, run_arguments, run_time, preliminary, fit_resolution_percent, fit_resolution_percent_stat, fit_resolution_percent_sys, std::pair<double, double>{5.0, 20.0});
		if (compare_to_other_tb) {
			draw_resolution_reference_comparison_page(inputs, fits, output_path, run_arguments, run_time, preliminary);
		}
	}
	canvas.Print((output_path + "]").c_str());
}

} // namespace

int main(int argc, char** argv)
{
	const std::string script_name = "003_ADC_Sum_Comparison";
	try {
		spdlog::info("Running script: {}", script_name);
		spdlog::info("----------------------------------------");

		cxxopts::Options options(script_name, "Compare all-lpGBT ADC sum histograms from 002_ADC_Analysis ROOT outputs");
		options.add_options()
			("c,config", "JSON comparison config", cxxopts::value<std::string>())
			("i,input", "Input event_adc_sum_histograms.root file(s); repeat --input for multiple files", cxxopts::value<std::vector<std::string>>())
			("o,output", "Output comparison PDF", cxxopts::value<std::string>()->default_value("dump/003/adc_sum_comparison.pdf"))
			("h,help", "Print usage");

		const auto parsed = options.parse(argc, argv);
		if (parsed.count("help")) {
			std::cout << options.help() << '\n';
			return 0;
		}

		std::string output_path = parsed["output"].as<std::string>();
		std::string title = "All lpGBT ADC sum comparison";
		std::string parameter_name;
		std::string parameter_unit;
		bool mu_linear_fit = false;
		bool compare_to_other_tb = false;
		bool preliminary = false;
		std::vector<InputHistogram> inputs;
		if (parsed.count("config")) {
			const auto config = load_config(parsed["config"].as<std::string>());
			if (!config.output_path.empty() && !parsed.count("output")) {
				output_path = config.output_path;
			}
			title = config.title;
			parameter_name = config.parameter_name;
			parameter_unit = config.parameter_unit;
			mu_linear_fit = config.mu_linear_fit;
			compare_to_other_tb = config.compare_to_other_tb;
			preliminary = config.preliminary;
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
			throw std::runtime_error("no input ROOT files provided and none found under dump/002/*/event_adc_sum_histograms.root");
		}
		if (inputs.empty()) {
			inputs = load_histograms(input_paths);
		}
		draw_comparison(inputs, output_path, title, parameter_name, parameter_unit, join_arguments(argc, argv), current_time_minute(), mu_linear_fit, compare_to_other_tb, preliminary);

		spdlog::info("Read {} ADC sum histogram ROOT files", inputs.size());
		spdlog::info("Wrote ADC sum comparison PDF to {}", output_path);
		return 0;
	} catch (const std::exception& error) {
		spdlog::error("{}", error.what());
		return 1;
	}
}
