#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>

Timer::Timer(const std::string& name) : name_(name) {
    reset();
}

Timer::~Timer() {
    std::cerr << name_ << " completed in " << elapsed_ms() << "ms" << std::endl;
}

void Timer::reset() {
    start_time_ = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed_ms() const {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_);
    return duration.count();
}

double Timer::elapsed_seconds() const {
    return elapsed_ms() / 1000.0;
}

size_t MemoryTracker::get_peak_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss * 1024; // Convert KB to bytes on Linux
}

size_t MemoryTracker::get_current_memory_usage() {
    std::ifstream file("/proc/self/status");
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.substr(0, 6) == "VmRSS:") {
            std::istringstream iss(line);
            std::string label, value, unit;
            iss >> label >> value >> unit;
            return std::stoul(value) * 1024; // Convert KB to bytes
        }
    }
    return 0;
}

void MemoryTracker::print_memory_stats() {
    size_t current = get_current_memory_usage();
    size_t peak = get_peak_memory_usage();
    
    std::cerr << "Estimated memory usage: " << current / (1024 * 1024) << " MB" << std::endl;
}

namespace utils {

bool file_exists(const std::string& filename) {
    std::ifstream file(filename.c_str());
    return file.good();
}

std::string error_code_to_string(ErrorCode code) {
    switch (code) {
        case ErrorCode::SUCCESS: return "Success";
        case ErrorCode::FILE_NOT_FOUND: return "File not found";
        case ErrorCode::PARSE_ERROR: return "Parse error";
        case ErrorCode::MEMORY_ERROR: return "Memory error";
        case ErrorCode::INVALID_FORMAT: return "Invalid format";
        default: return "Unknown error";
    }
}

} // namespace utils
