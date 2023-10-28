#pragma once
#include <cmath>
#include <vector>
#include <math.h> 
#include <string>
#include <unordered_map>
#include <fstream>
#include <deque>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Iir.h"
#include "tools.h"

namespace p = pybind11;

class FindPeakRealTime {
public:
    FindPeakRealTime(double fs, double filter, double sigma, int resolution, int window, int k);
    p::tuple append_filter(std::vector<float>& data, int skip);
    p::tuple append_nofilter(std::vector<float>& data, int skip);
    p::tuple append(std::vector<float>& data, int skip) { return (this->*_append)(data, skip); };
    void reset();
    void flush();

private:
    p::tuple(FindPeakRealTime::*_append)(std::vector<float>&, int);
    double data = 0;
    double datapre = 0;
    double datasuf = 0;
    double fs = 0;
    double filter = 0;
    double th = 0;
    double sigma = 0;
    size_t localminindex = 0;
    size_t localmaxindex = 0;
    double median = 0;
    double medianold = 0;
    double sd = 0;
    double stddev = 0;
    double mean = 0;
    double meanold = 0;
    double d2 = 0;
    double d = 0;
    double z = 0;
    double min = 0;
    double max = 0;
    size_t s = 0;
    size_t e = 0;
    size_t i = 0;
    int k = -1;
    int windowsizestart = 0;
    int windowsizeend = 0;
    bool flag = false;
    bool s0flag = false;
    double interval = 0;
    std::deque<std::pair<double, size_t>> quemin;
    std::deque<std::pair<double, size_t>> quemax;
    std::vector<double> peakcurrent;
    std::deque<double> windowend;
    std::deque<std::pair<double, size_t>> windowstart;
    std::vector<Peak> result;
    Iir::Butterworth::LowPass<5> f;
};

