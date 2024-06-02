#pragma once
#include <chrono>

class FPSTimer {
private:
    size_t frames_per_second = 0;
    size_t curr_frames = 0;
    std::chrono::time_point<std::chrono::steady_clock> old_time = std::chrono::steady_clock::now();
public:
    void update() {
        ++curr_frames;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - old_time) >= std::chrono::milliseconds(1000)) {
            old_time = std::chrono::steady_clock::now();
            frames_per_second = curr_frames;
            curr_frames = 0;
        }
    }
    size_t get_fps() {
        return frames_per_second;
    }
};

class Timer {
private:
    std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
public:
    void reset() {
        start_time = std::chrono::steady_clock::now();
    }
    
    double elapsed_s() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    }

    double elapsed_ms() {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start_time).count();
    }
};