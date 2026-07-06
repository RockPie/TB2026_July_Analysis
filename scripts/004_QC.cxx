#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>

#include "cxxopts.hpp"

#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct QcPoint {
	std::string run_name;
	std::time_t timestamp = 0;
	double reconstructed_percent = 0.0;
	double crc_percent = 0.0;
	double all_lpgbt_peak_events = 0.0;
	std::map<std::string, double> line_percents;
	bool has_reconstructed_percent = false;
	bool has_crc_percent = false;
	bool has_all_lpgbt_peak_events = false;
};

std::string read_text_file(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("failed to open log file: " + path.string());
	}
	std::ostringstream buffer;
	buffer << input.rdbuf();
	return buffer.str();
}

std::string remove_commas(std::string text)
{
	text.erase(std::remove(text.begin(), text.end(), ','), text.end());
	return text;
}

std::time_t parse_run_timestamp(const std::string& run_name)
{
	const std::regex pattern(R"(run__(\d{4})_(\d{2})_(\d{2})__(\d{2})_(\d{2})_(\d{2})__)");
	std::smatch match;
	if (!std::regex_search(run_name, match, pattern)) {
		throw std::runtime_error("log filename does not contain run timestamp: " + run_name);
	}
	std::tm tm{};
	tm.tm_year = std::stoi(match[1].str()) - 1900;
	tm.tm_mon = std::stoi(match[2].str()) - 1;
	tm.tm_mday = std::stoi(match[3].str());
	tm.tm_hour = std::stoi(match[4].str());
	tm.tm_min = std::stoi(match[5].str());
	tm.tm_sec = std::stoi(match[6].str());
	tm.tm_isdst = -1;
	return std::mktime(&tm);
}

std::string run_name_from_log_path(const std::filesystem::path& path)
{
	return path.stem().string();
}

bool parse_log_file(const std::filesystem::path& path, QcPoint& point)
{
	point.run_name = run_name_from_log_path(path);
	point.timestamp = parse_run_timestamp(point.run_name);
	const auto text = read_text_file(path);

	const std::regex reconstructed_pattern(R"(All lpGBT sample peak events(?: \(CRC OK\))?\s*:\s*([0-9,]+)\s*\(\s*([0-9]+(?:\.[0-9]+)?)%\))");
	std::smatch match;
	if (std::regex_search(text, match, reconstructed_pattern)) {
		point.all_lpgbt_peak_events = std::stod(remove_commas(match[1].str()));
		point.reconstructed_percent = std::stod(match[2].str());
		point.has_all_lpgbt_peak_events = true;
		point.has_reconstructed_percent = true;
	}

	const std::regex crc_pattern(R"(40-line all-elink CRC OK\s*:\s*([0-9,]+)\s*\(\s*([0-9]+(?:\.[0-9]+)?)%\))");
	if (std::regex_search(text, match, crc_pattern)) {
		point.crc_percent = std::stod(match[2].str());
		point.has_crc_percent = true;
	}

	const std::vector<std::string> line_types{
		"heartbeat line",
		"phy trigger line",
		"data idle line",
		"data daq line",
		"dummy zero line",
		"other line",
	};
	for (const auto& line_type : line_types) {
		const std::regex line_pattern(line_type + R"(\s*:\s*[0-9,]+\s*\(\s*([0-9]+(?:\.[0-9]+)?)%\))");
		if (std::regex_search(text, match, line_pattern)) {
			point.line_percents[line_type] = std::stod(match[1].str());
		}
	}

	return point.has_reconstructed_percent && point.has_crc_percent && point.has_all_lpgbt_peak_events && point.line_percents.size() == line_types.size();
}

std::vector<QcPoint> load_qc_points(const std::filesystem::path& log_dir)
{
	std::vector<QcPoint> points;
	if (!std::filesystem::exists(log_dir)) {
		throw std::runtime_error("log directory does not exist: " + log_dir.string());
	}
	for (const auto& entry : std::filesystem::directory_iterator(log_dir)) {
		if (!entry.is_regular_file() || entry.path().extension() != ".log") {
			continue;
		}
		QcPoint point;
		try {
			if (parse_log_file(entry.path(), point)) {
				points.push_back(std::move(point));
			}
		} catch (const std::exception& error) {
			spdlog::warn("Skipping {}: {}", entry.path().string(), error.what());
		}
	}
	std::sort(points.begin(), points.end(), [](const auto& left, const auto& right) {
		return left.timestamp < right.timestamp;
	});
	if (points.empty()) {
		throw std::runtime_error("no complete 001 reconstruction statistics found under: " + log_dir.string());
	}
	return points;
}

std::string format_run_time(std::time_t timestamp)
{
	std::tm tm{};
	localtime_r(&timestamp, &tm);
	std::ostringstream stream;
	stream << std::put_time(&tm, "%m-%d %H:%M");
	return stream.str();
}

std::vector<double> time_hours_from_first(const std::vector<QcPoint>& points)
{
	std::vector<double> hours;
	hours.reserve(points.size());
	const auto first = points.front().timestamp;
	for (const auto& point : points) {
		hours.push_back(std::difftime(point.timestamp, first) / 3600.0);
	}
	return hours;
}

void draw_header(const std::string& title, const std::string& subtitle)
{
	TLatex latex;
	latex.SetNDC();
	latex.SetTextFont(42);
	latex.SetTextAlign(13);
	latex.SetTextSize(0.042);
	latex.DrawLatex(0.12, 0.88, "#bf{FoCal TB2026 July QC}");
	latex.SetTextSize(0.032);
	latex.DrawLatex(0.12, 0.835, title.c_str());
	latex.DrawLatex(0.12, 0.795, subtitle.c_str());
}

void label_time_axis(TAxis& axis, const std::vector<QcPoint>& points, const std::vector<double>& hours)
{
	axis.SetTitle("run time");
	axis.SetLabelSize(0.032);
	axis.SetTitleSize(0.045);
	axis.SetNdivisions(static_cast<int>(std::min<std::size_t>(points.size(), 8)), false);
	for (std::size_t index = 0; index < points.size() && index < hours.size(); ++index) {
		const auto bin = axis.FindBin(hours[index]);
		axis.SetBinLabel(bin, format_run_time(points[index].timestamp).c_str());
	}
}

void draw_single_series(TCanvas& canvas, const std::string& pdf_path, const std::vector<QcPoint>& points,
	const std::vector<double>& x, const std::vector<double>& y, const std::string& title, const std::string& y_title,
	int color, bool percent_range)
{
	canvas.Clear();
	canvas.SetLogy(0);
	canvas.SetLeftMargin(0.12);
	canvas.SetRightMargin(0.08);
	canvas.SetBottomMargin(0.16);
	TGraph graph(static_cast<int>(x.size()), x.data(), y.data());
	graph.SetTitle((";run time;" + y_title).c_str());
	graph.SetMarkerStyle(20);
	graph.SetMarkerSize(1.1);
	graph.SetMarkerColor(color);
	graph.SetLineColor(color);
	graph.SetLineWidth(2);
	graph.Draw("APL");
	if (percent_range) {
		graph.GetYaxis()->SetRangeUser(0.0, 105.0);
	}
	label_time_axis(*graph.GetXaxis(), points, x);
	graph.GetXaxis()->LabelsOption("v");
	graph.GetYaxis()->SetTitleSize(0.045);
	graph.GetYaxis()->SetLabelSize(0.04);
	draw_header(title, "source: dump/logs/001/*.log");
	canvas.Print(pdf_path.c_str());
}

void draw_line_type_page(TCanvas& canvas, const std::string& pdf_path, const std::vector<QcPoint>& points, const std::vector<double>& x)
{
	canvas.Clear();
	canvas.SetLogy(0);
	canvas.SetLeftMargin(0.12);
	canvas.SetRightMargin(0.20);
	canvas.SetBottomMargin(0.16);

	const std::vector<std::pair<std::string, int>> series{
		{"heartbeat line", kCyan + 2},
		{"phy trigger line", kOrange + 7},
		{"data idle line", kGray + 2},
		{"data daq line", kGreen + 2},
		{"dummy zero line", kBlue + 1},
		{"other line", kRed + 1},
	};
	std::vector<std::unique_ptr<TGraph>> graphs;
	TLegend legend(0.82, 0.56, 0.98, 0.86);
	legend.SetBorderSize(0);
	legend.SetFillStyle(0);
	legend.SetTextSize(0.027);
	for (std::size_t series_index = 0; series_index < series.size(); ++series_index) {
		std::vector<double> y;
		y.reserve(points.size());
		for (const auto& point : points) {
			y.push_back(point.line_percents.at(series[series_index].first));
		}
		auto graph = std::make_unique<TGraph>(static_cast<int>(x.size()), x.data(), y.data());
		graph->SetTitle(";run time;line fraction [%]");
		graph->SetMarkerStyle(20 + static_cast<int>(series_index));
		graph->SetMarkerSize(1.0);
		graph->SetMarkerColor(series[series_index].second);
		graph->SetLineColor(series[series_index].second);
		graph->SetLineWidth(2);
		graph->Draw(series_index == 0 ? "APL" : "PL same");
		if (series_index == 0) {
			graph->GetYaxis()->SetRangeUser(0.0, 105.0);
			label_time_axis(*graph->GetXaxis(), points, x);
			graph->GetXaxis()->LabelsOption("v");
			graph->GetYaxis()->SetTitleSize(0.045);
			graph->GetYaxis()->SetLabelSize(0.04);
		}
		legend.AddEntry(graph.get(), series[series_index].first.c_str(), "lp");
		graphs.push_back(std::move(graph));
	}
	legend.Draw();
	draw_header("001 line type fractions", "source: dump/logs/001/*.log");
	canvas.Print(pdf_path.c_str());
}

void write_qc_pdf(const std::vector<QcPoint>& points, const std::string& output_path)
{
	std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	const auto x = time_hours_from_first(points);
	TCanvas canvas("qc_canvas", "001 reconstruction QC", 1300, 850);
	canvas.Print((output_path + "[").c_str());

	std::vector<double> reconstructed_percent;
	std::vector<double> crc_percent;
	std::vector<double> peak_events;
	for (const auto& point : points) {
		reconstructed_percent.push_back(point.reconstructed_percent);
		crc_percent.push_back(point.crc_percent);
		peak_events.push_back(point.all_lpgbt_peak_events);
	}
	draw_single_series(canvas, output_path, points, x, reconstructed_percent,
		"All lpGBT sample peak event fraction", "reconstructed event fraction [%]", kGreen + 2, true);
	draw_single_series(canvas, output_path, points, x, crc_percent,
		"40-line all-elink CRC OK fraction", "CRC OK fraction [%]", kBlue + 1, true);
	draw_line_type_page(canvas, output_path, points, x);
	draw_single_series(canvas, output_path, points, x, peak_events,
		"All lpGBT sample peak event count", "reconstructed event count", kMagenta + 1, false);

	canvas.Print((output_path + "]").c_str());
}

} // namespace

int main(int argc, char** argv)
{
	const std::string script_name = "004_QC";
	try {
		spdlog::info("Running script: {}", script_name);
		spdlog::info("----------------------------------------");
		cxxopts::Options options(script_name, "Draw QC trends from 001_Rootifier logs");
		options.add_options()
			("l,log-dir", "Directory containing 001_Rootifier logs", cxxopts::value<std::string>()->default_value("dump/logs/001"))
			("o,output", "Output QC PDF", cxxopts::value<std::string>()->default_value("dump/004/qc_summary.pdf"))
			("h,help", "Print usage");
		const auto parsed = options.parse(argc, argv);
		if (parsed.count("help")) {
			std::cout << options.help() << '\n';
			return 0;
		}
		const auto log_dir = parsed["log-dir"].as<std::string>();
		const auto output_path = parsed["output"].as<std::string>();
		const auto points = load_qc_points(log_dir);
		write_qc_pdf(points, output_path);
		spdlog::info("Read {} complete 001 log files from {}", points.size(), log_dir);
		spdlog::info("Wrote QC PDF to {}", output_path);
		return 0;
	} catch (const std::exception& error) {
		spdlog::error("{}", error.what());
		return 1;
	}
}