#pragma once
#include <cmath>
#include <vector>
#include <math.h> 
#include <string>
#include <unordered_map>
#include <fstream>
#include <deque>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


namespace p = pybind11;

//void moveBaseline(double* x, double* b, int size, double mean, double sd, double thdiff, double maxmeandiff, double maxstddiff);
//gsl_vector* meanSmooth(gsl_vector* x, int w);
//std::pair<double, double> valminmax(float* pos, int start, int end);
//void cumSum(double* x, int n, int s, int e, double stepsize, double h, double stdDev, std::vector<int>& v1, std::vector<double>& v2);
//std::vector<int> extremeval(double* x, int s, int e, double stdDev, double stddiff);
//std::vector<int> extreminval(double* x, int s, int e);
//void stdDev(double* x, int size, double& stdDev);
//double meancurrent(double* x, int size, int s, int e);
//int argmin(double* x, int s, int e);
//void setmovebaseline(double* buff, double* x, int s, int e);
//void setcurrent(double* x, int s, int e);
//double stdvariance(double* x, int size, int s, int e, double baseline);
//void stdDevBaseline(double* x, double* x_d, int size, double& stdDev, double& baseline);

struct Peak {
    double currentpeak;
    double baseline;
    int s;
    int e;
    int s0;
};


std::list<Peak> findPeak_longevent(std::vector<float>&, int, int, int, int, int, int);
std::list<Peak> findPeak_shortevent(std::vector<float>&, int, int, int, int, int, int);
std::vector<float> Iir_filter(std::vector<float>&, int, int);
//p::list CUSUM(std::vector<float>&, double, double, double);


