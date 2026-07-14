cmake_minimum_required(VERSION 3.20)

project(TB2026_July_Analysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_PREFIX_PATH
	/opt/homebrew
	/usr/local
)

find_package(ROOT REQUIRED COMPONENTS
	Core
	RIO
	Tree
	Hist
	Graf
	Gpad
	Matrix
	RHTTP
)

find_package(cxxopts CONFIG QUIET)

set(missing_dependencies)

find_path(SPDLOG_INCLUDE_DIR
	NAMES spdlog/spdlog.h
	HINTS /opt/homebrew/include /usr/local/include
)
if(SPDLOG_INCLUDE_DIR)
	add_library(spdlog_header_only INTERFACE)
	target_include_directories(spdlog_header_only INTERFACE ${SPDLOG_INCLUDE_DIR})
	target_compile_definitions(spdlog_header_only INTERFACE SPDLOG_HEADER_ONLY FMT_HEADER_ONLY)
else()
	list(APPEND missing_dependencies "spdlog (install spdlog or set SPDLOG_INCLUDE_DIR)")
endif()

if(NOT TARGET cxxopts::cxxopts)
	find_path(CXXOPTS_INCLUDE_DIR
		NAMES cxxopts.hpp
		HINTS /opt/homebrew/include /usr/local/include
	)
	if(CXXOPTS_INCLUDE_DIR)
		add_library(cxxopts_header_only INTERFACE)
		add_library(cxxopts::cxxopts ALIAS cxxopts_header_only)
		target_include_directories(cxxopts_header_only INTERFACE ${CXXOPTS_INCLUDE_DIR})
	else()
		list(APPEND missing_dependencies "cxxopts.hpp (install cxxopts or set cxxopts_DIR/CXXOPTS_INCLUDE_DIR)")
	endif()
endif()

find_path(BINPARSE_INCLUDE_DIR
	NAMES binparse/parser.hpp
	HINTS ${PROJECT_SOURCE_DIR}/include /opt/homebrew/include /usr/local/include
)
if(NOT BINPARSE_INCLUDE_DIR)
	list(APPEND missing_dependencies "binparse headers (place binparse/ under include/ or set BINPARSE_INCLUDE_DIR)")
endif()

find_path(ZSTD_INCLUDE_DIR
	NAMES zstd.h
	HINTS /opt/homebrew/include /usr/local/include /usr/include
)
find_library(ZSTD_LIBRARY
	NAMES zstd
	HINTS /opt/homebrew/lib /usr/local/lib /usr/lib /usr/lib/x86_64-linux-gnu
)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
	list(APPEND missing_dependencies "zstd (install libzstd development headers/library)")
endif()

if(missing_dependencies)
	list(JOIN missing_dependencies "\n  - " missing_dependencies_text)
	message(FATAL_ERROR "Missing dependencies for 001_Rootifier:\n  - ${missing_dependencies_text}")
endif()

add_executable(001_Rootifier
	scripts/001_Rootifier.cxx
)

add_executable(002_ADC_Analysis
	scripts/002_ADC_Analysis.cxx
)

add_executable(003_ADC_Sum_Comparison
	scripts/003_ADC_Sum_Comparison.cxx
)

add_executable(004_QC
	scripts/004_QC.cxx
)

add_executable(005_ToT_Analysis
	scripts/005_ToT_Analysis.cxx
)

add_executable(006_ToT_Sum_Comparison
	scripts/006_ToT_Sum_Comparison.cxx
)

add_executable(007_ToA_Analysis
	scripts/007_ToA_Analysis.cxx
)

set_target_properties(001_Rootifier PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(002_ADC_Analysis PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(003_ADC_Sum_Comparison PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(004_QC PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(005_ToT_Analysis PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(006_ToT_Sum_Comparison PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(007_ToA_Analysis PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

target_include_directories(001_Rootifier PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
	${ZSTD_INCLUDE_DIR}
)

target_compile_options(001_Rootifier PRIVATE
	$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,GNU>:-O3>
)

target_include_directories(002_ADC_Analysis PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_compile_options(002_ADC_Analysis PRIVATE
	$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,GNU>:-O3>
)

target_include_directories(003_ADC_Sum_Comparison PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_include_directories(004_QC PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_include_directories(005_ToT_Analysis PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_include_directories(006_ToT_Sum_Comparison PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_include_directories(007_ToA_Analysis PRIVATE
	${PROJECT_SOURCE_DIR}/include
	${BINPARSE_INCLUDE_DIR}
)

target_link_libraries(001_Rootifier PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	ROOT::Matrix
	ROOT::RHTTP
	${ZSTD_LIBRARY}
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(002_ADC_Analysis PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(003_ADC_Sum_Comparison PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(004_QC PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(005_ToT_Analysis PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(006_ToT_Sum_Comparison PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)

target_link_libraries(007_ToA_Analysis PRIVATE
	ROOT::Core
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	spdlog_header_only
	cxxopts::cxxopts
)
