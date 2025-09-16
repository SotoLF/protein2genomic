#ifndef UTILS_HPP
#define UTILS_HPP

#include "common.hpp"
#include <chrono>

class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time_;
    std::string name_;
    
public:
    explicit Timer(const std::string& name = "Timer");
    ~Timer();
    
    void reset();
    double elapsed_ms() const;
    double elapsed_seconds() const;
};

class MemoryTracker {
public:
    static size_t get_peak_memory_usage();
    static size_t get_current_memory_usage();
    static void print_memory_stats();
};

namespace utils {
    // File utilities
    bool file_exists(const std::string& filename);
    
    // Error handling
    std::string error_code_to_string(ErrorCode code);
}

#endif // UTILS_HPP
