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
#include <memory>
#include <optional>
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
	double run_info_events = 0.0;
	double run_info_triggers = 0.0;
	std::map<std::string, double> line_percents;
	bool has_reconstructed_percent = false;
	bool has_crc_percent = false;
	bool has_all_lpgbt_peak_events = false;
	bool has_run_info_events = false;
	bool has_run_info_triggers = false;
};

struct RunInfoRow {
	std::time_t start_timestamp = 0;
	std::time_t end_timestamp = 0;
	double events = 0.0;
	double triggers = 0.0;
	bool has_events = false;
	bool has_triggers = false;
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

std::string trim(std::string text)
{
	const auto first = text.find_first_not_of(" \t\r\n");
	if (first == std::string::npos) {
		return {};
	}
	const auto last = text.find_last_not_of(" \t\r\n");
	return text.substr(first, last - first + 1);
}

std::vector<std::string> split_csv_line(const std::string& line)
{
	std::vector<std::string> fields;
	std::string field;
	bool in_quotes = false;
	for (const auto character : line) {
		if (character == '"') {
			in_quotes = !in_quotes;
			continue;
		}
		if (character == ',' && !in_quotes) {
			fields.push_back(trim(field));
			field.clear();
			continue;
		}
		field.push_back(character);
	}
	fields.push_back(trim(field));
	return fields;
}

std::optional<double> parse_optional_count(std::string text)
{
	text = trim(text);
	if (text.empty() || text == "-") {
		return std::nullopt;
	}
	text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char character) {
		return character == ',' || std::isspace(character);
	}), text.end());
	if (text.empty() || text == "-") {
		return std::nullopt;
	}
	return std::stod(text);
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

std::optional<std::time_t> parse_run_info_timestamp(const std::string& date_text, const std::string& time_text)
{
	const auto date = trim(date_text);
	const auto time = trim(time_text);
	if (date.empty() || time.empty()) {
		return std::nullopt;
	}
	const std::regex date_pattern(R"((\d{1,2})/(\d{1,2})/(\d{2,4}))");
	const std::regex time_pattern(R"((\d{1,2}):(\d{2}))");
	std::smatch date_match;
	std::smatch time_match;
	if (!std::regex_match(date, date_match, date_pattern) || !std::regex_match(time, time_match, time_pattern)) {
		return std::nullopt;
	}
	int year = std::stoi(date_match[3].str());
	if (year < 100) {
		year += 2000;
	}
	std::tm tm{};
	tm.tm_mday = std::stoi(date_match[1].str());
	tm.tm_mon = std::stoi(date_match[2].str()) - 1;
	tm.tm_year = year - 1900;
	tm.tm_hour = std::stoi(time_match[1].str());
	tm.tm_min = std::stoi(time_match[2].str());
	tm.tm_sec = 0;
	tm.tm_isdst = -1;
	return std::mktime(&tm);
}

std::vector<RunInfoRow> load_run_info_rows(const std::filesystem::path& path)
{
	std::vector<RunInfoRow> rows;
	if (!std::filesystem::exists(path)) {
		spdlog::warn("Run info CSV does not exist: {}", path.string());
		return rows;
	}
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("failed to open run info CSV: " + path.string());
	}
	std::string line;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() < 6 || fields[0] == "Date") {
			continue;
		}
		const auto start = parse_run_info_timestamp(fields[0], fields[1]);
		const auto end = parse_run_info_timestamp(fields[0], fields[2]);
		if (!start.has_value() || !end.has_value()) {
			continue;
		}
		RunInfoRow row;
		row.start_timestamp = *start;
		row.end_timestamp = *end;
		if (row.end_timestamp < row.start_timestamp) {
			row.end_timestamp += 24 * 60 * 60;
		}
		if (const auto events = parse_optional_count(fields[4])) {
			row.events = *events;
			row.has_events = true;
		}
		if (const auto triggers = parse_optional_count(fields[5])) {
			row.triggers = *triggers;
			row.has_triggers = true;
		}
		if (row.has_events || row.has_triggers) {
			rows.push_back(row);
		}
	}
	std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
		return left.start_timestamp < right.start_timestamp;
	});
	return rows;
}

void attach_run_info_counts(std::vector<QcPoint>& points, const std::vector<RunInfoRow>& rows)
{
	for (auto& point : points) {
		const auto row = std::find_if(rows.begin(), rows.end(), [&point](const auto& candidate) {
			return candidate.start_timestamp <= point.timestamp && point.timestamp <= candidate.end_timestamp;
		});
		if (row == rows.end()) {
			continue;
		}
		if (row->has_events) {
			point.run_info_events = row->events;
			point.has_run_info_events = true;
		}
		if (row->has_triggers) {
			point.run_info_triggers = row->triggers;
			point.has_run_info_triggers = true;
		}
	}
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

void apply_qc_margins(TCanvas& canvas)
{
	canvas.SetLeftMargin(0.13);
	canvas.SetRightMargin(0.08);
	canvas.SetBottomMargin(0.16);
}

void draw_header(const std::string& title, const std::string& subtitle)
{
	TLatex latex;
	latex.SetNDC();
	latex.SetTextFont(42);
	latex.SetTextAlign(13);
	latex.SetTextSize(0.042);
	latex.DrawLatex(0.145, 0.88, "#bf{FoCal TB2026 July QC}");
	latex.SetTextSize(0.032);
	latex.DrawLatex(0.145, 0.835, title.c_str());
	latex.DrawLatex(0.145, 0.795, subtitle.c_str());
}

void label_time_axis(TAxis& axis, const std::vector<QcPoint>& points, const std::vector<double>& hours)
{
	axis.SetTitle("run time");
	axis.SetLabelSize(0.032);
	axis.SetTitleSize(0.045);
	axis.SetTitleOffset(1.25);
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
	apply_qc_margins(canvas);
	TGraph graph(static_cast<int>(x.size()), x.data(), y.data());
	graph.SetTitle((";run time;" + y_title).c_str());
	graph.SetMarkerStyle(20);
	graph.SetMarkerSize(1.1);
	graph.SetMarkerColor(color);
	graph.SetLineColor(color);
	graph.SetLineWidth(2);
	graph.Draw("APL");
	if (percent_range) {
		graph.GetYaxis()->SetRangeUser(0.0, 130.0);
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
	apply_qc_margins(canvas);

	const std::vector<std::pair<std::string, int>> series{
		{"heartbeat line", kCyan + 2},
		{"phy trigger line", kOrange + 7},
		{"data idle line", kGray + 2},
		{"data daq line", kGreen + 2},
		{"dummy zero line", kBlue + 1},
		{"other line", kRed + 1},
	};
	std::vector<std::unique_ptr<TGraph>> graphs;
	TLegend legend(0.70, 0.64, 0.90, 0.88);
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
			graph->GetYaxis()->SetRangeUser(0.0, 130.0);
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

void draw_reconstructed_event_count_page(TCanvas& canvas, const std::string& pdf_path, const std::vector<QcPoint>& points, const std::vector<double>& x)
{
	canvas.Clear();
	canvas.SetLogy(0);
	apply_qc_margins(canvas);

	std::vector<double> reconstructed_events;
	std::vector<double> run_info_events;
	std::vector<double> run_info_triggers;
	reconstructed_events.reserve(points.size());
	run_info_events.reserve(points.size());
	run_info_triggers.reserve(points.size());
	double max_y = 0.0;
	for (const auto& point : points) {
		reconstructed_events.push_back(point.all_lpgbt_peak_events);
		run_info_events.push_back(point.has_run_info_events ? point.run_info_events : 0.0);
		run_info_triggers.push_back(point.has_run_info_triggers ? point.run_info_triggers : 0.0);
		max_y = std::max(max_y, point.all_lpgbt_peak_events);
		if (point.has_run_info_events) {
			max_y = std::max(max_y, point.run_info_events);
		}
		if (point.has_run_info_triggers) {
			max_y = std::max(max_y, point.run_info_triggers);
		}
	}
	const auto y_max = max_y > 0.0 ? max_y * 1.4 : 1.0;

	TGraph reconstructed_graph(static_cast<int>(x.size()), x.data(), reconstructed_events.data());
	TGraph events_graph(static_cast<int>(x.size()), x.data(), run_info_events.data());
	TGraph triggers_graph(static_cast<int>(x.size()), x.data(), run_info_triggers.data());
	std::vector<TGraph*> graphs{&reconstructed_graph, &events_graph, &triggers_graph};
	const std::vector<int> colors{kMagenta + 1, kGreen + 2, kBlue + 1};
	const std::vector<int> markers{20, 21, 22};
	for (std::size_t index = 0; index < graphs.size(); ++index) {
		graphs[index]->SetTitle(";run time;event count");
		graphs[index]->SetMarkerStyle(markers[index]);
		graphs[index]->SetMarkerSize(1.0);
		graphs[index]->SetMarkerColor(colors[index]);
		graphs[index]->SetLineColor(colors[index]);
		graphs[index]->SetLineWidth(2);
		graphs[index]->Draw(index == 0 ? "APL" : "PL same");
	}
	reconstructed_graph.GetYaxis()->SetRangeUser(0.0, y_max);
	label_time_axis(*reconstructed_graph.GetXaxis(), points, x);
	reconstructed_graph.GetXaxis()->LabelsOption("v");
	reconstructed_graph.GetYaxis()->SetTitleSize(0.045);
	reconstructed_graph.GetYaxis()->SetLabelSize(0.04);

	TLegend legend(0.62, 0.70, 0.90, 0.88);
	legend.SetBorderSize(0);
	legend.SetFillStyle(0);
	legend.SetTextSize(0.030);
	legend.AddEntry(&reconstructed_graph, "reconstructed events", "lp");
	legend.AddEntry(&events_graph, "Run_Info number of events", "lp");
	legend.AddEntry(&triggers_graph, "Run_Info number of triggers", "lp");
	legend.Draw();
	draw_header("All lpGBT sample peak event count", "source: dump/logs/001/*.log and config/Run_Info.csv");
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
	for (const auto& point : points) {
		reconstructed_percent.push_back(point.reconstructed_percent);
		crc_percent.push_back(point.crc_percent);
	}
	draw_single_series(canvas, output_path, points, x, reconstructed_percent,
		"All lpGBT sample peak event fraction", "reconstructed event fraction [%]", kGreen + 2, true);
	draw_single_series(canvas, output_path, points, x, crc_percent,
		"40-line all-elink CRC OK fraction", "CRC OK fraction [%]", kBlue + 1, true);
	draw_line_type_page(canvas, output_path, points, x);
	draw_reconstructed_event_count_page(canvas, output_path, points, x);

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
			("r,run-info", "Run information CSV with start/end times and event/trigger counts", cxxopts::value<std::string>()->default_value("config/Run_Info.csv"))
			("o,output", "Output QC PDF", cxxopts::value<std::string>()->default_value("dump/004/qc_summary.pdf"))
			("h,help", "Print usage");
		const auto parsed = options.parse(argc, argv);
		if (parsed.count("help")) {
			std::cout << options.help() << '\n';
			return 0;
		}
		const auto log_dir = parsed["log-dir"].as<std::string>();
		const auto run_info_path = parsed["run-info"].as<std::string>();
		const auto output_path = parsed["output"].as<std::string>();
		auto points = load_qc_points(log_dir);
		const auto run_info_rows = load_run_info_rows(run_info_path);
		attach_run_info_counts(points, run_info_rows);
		write_qc_pdf(points, output_path);
		spdlog::info("Read {} complete 001 log files from {}", points.size(), log_dir);
		spdlog::info("Read {} usable run info rows from {}", run_info_rows.size(), run_info_path);
		spdlog::info("Wrote QC PDF to {}", output_path);
		return 0;
	} catch (const std::exception& error) {
		spdlog::error("{}", error.what());
		return 1;
	}
}