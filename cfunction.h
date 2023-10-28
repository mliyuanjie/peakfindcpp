#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <fstream>
#include <thread>
#include <cmath>
#include <time.h>
#include "datio.h"

namespace p = pybind11;

template <typename T>
p::dict findevent(p::array_t<T> x_py, int resolution, int window, int sigma, bool direction, bool islongevent);
std::vector<float> test(const std::vector<float>&, int);
p::tuple randomWalk(int, double, double, double, double, double, double, double, double, double, double, double, double, double, int);
p::tuple randomWalkDt(double, double, double, double, double, double, double, double, double, double);
std::vector<double> randomAngleWalk(int, double, double, double);
void randomangle_thread(float* data, float* angle0, float* dipolefield, float* dangle, int m, int i, int end, int skips);
void randomAngleWalkParallel(p::array_t<float>&, p::array_t<float>&, p::array_t<float>&, p::array_t<float>&, int, int, int);
template <typename T>
class Downsampler {
public:
	Downsampler(p::array_t<T> x_py, int samplingrate);
	Downsampler(const std::string& fn, int samplingrate);
	~Downsampler();
	p::tuple downsample(double, double);

private:
	DATIO* dat;
	double interval = 0;
};


