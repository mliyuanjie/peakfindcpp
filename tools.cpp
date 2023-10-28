#include <numeric>
#include <algorithm>
#include <deque>
#include <time.h>
#include <cmath>
#include <vector>
#include <random>
#include "tools.h"
#include "Iir.h"

void calsd(std::vector<float>& data, int s, int e, float& sd) {
    if (e - s < 6000) {
        double m = std::accumulate(data.begin() + s, data.begin() + e, 0.0) / (e - s);
        double accum = 0.0;
        std::for_each(data.begin() + s, data.begin() + e, [&](const double d) {
            accum += (d - m) * (d - m);
            });
        sd = sqrt(accum / (e - s));
        return;
    }
    std::random_device seed;
    std::mt19937 gen{ seed() };
    std::uniform_int_distribution<int> dist{ s, e - 5000};
    int n = 0;
    while (n < 20) {
        auto iter = data.begin() + dist(gen);
        double sum = std::accumulate(iter, std::next(iter, 5000), 0.0);
        double m = sum / 5000;
        double accum = 0.0;
        std::for_each(iter, std::next(iter, 5000), [&](const double d) {
            accum += (d - m) * (d - m);
            });
        float stdev = sqrt(accum / (5000 - 1));
        if (n == 0) {
            sd = stdev;
        }
        else {
            sd = (stdev < sd) ? stdev : sd;
        }
        n++;
    }
    return;
}

float meancurrent(std::vector<float>& data, int s, int e) {
    float res = 0;
    for (int i = s; i <= e; i++) {
        res = (res * (i - s) + data[i]) / (i - s + 1);
    }
    return res;
}


std::vector<int> extreminval(std::vector<float>& x, int s, int e, int k) {
    std::vector<int> l;
    if (s == 0 || e == x.size())
        return l;

    for (int i = s; i <= e; i++) {
        if (k * x[i] <= k * x[i - 1] && k * x[i] <= k * x[i + 1])
            l.push_back(i);
    }
    return l;
}


std::list<Peak> findPeak_longevent(std::vector<float>& data, int s, int e, int resolution, int window, int sigma, int k) {
    //stddev calculate
    float stdev = 0;
    calsd(data, s, e, stdev);
    float th = sigma * stdev;

    std::list<Peak> res;
    int localmin = 0;
    bool flag = false;
    std::deque<std::pair<float, int>> quemin;
    std::deque<std::pair<float, int>> quemax;
    std::deque<int> medianwindow;
    float median = 0;
    float mean = 0;
    float meanold = 0;
    float d2 = 0;
    float d = 0;
    float z = 0;
    int ps = s;
    int pe = s;
    float min = 0;
    float max = 0;


    for (int i = s; i < e - 1; i++) {
        //cal moving mean for find start of peak
        if (!flag) {
            if (medianwindow.size() < resolution) {
                medianwindow.push_back(i);
                median += (data[i] - median) / medianwindow.size();

            }
            else {
                median += (data[i] - data[medianwindow.front()]) / resolution;
                medianwindow.pop_front();
                medianwindow.push_back(i);
            }
        }

        //maintain min by dequene
        while (!quemin.empty() && quemin.back().first >= k * data[i])
            quemin.pop_back();
        quemin.push_back(std::pair<float, int>(k * data[i], i));
        while (quemin.front().second <= i - window)
            quemin.pop_front();
        //time1 += double(clock() - times)/ CLOCKS_PER_SEC;
        min = quemin.front().first;
        while (!quemax.empty() && quemax.back().first <= k * data[i])
            quemax.pop_back();
        quemax.push_back(std::pair<float, int>(k * data[i], i));
        while (quemax.front().second <= i - window)
            quemax.pop_front();
        max = quemax.front().first;

        //maintain sliding std and mean, mean, d, d2
        if (i - s <= window - 1) {
            meanold = mean;
            mean += (data[i] - mean) / (i - s + 1);
            d2 += (data[i] - meanold) * (data[i] - mean);
            continue;
        }
        else {
            meanold = mean;
            mean += (data[i] - data[i - window]) / window;
            d2 += (data[i] - data[i - window]) * (data[i] + data[i - window] - mean - meanold);
            d = sqrt(d2 / window);
        }

        //float z = 0;

        if (!flag) {
            if (k * data[i] < k * median - th) {
                for (int j = i; j > pe && j == medianwindow.back(); j--) {
                    if (k * data[j] >= k * data[j - 1] && k * data[j] >= k * data[j + 1] && k * data[j] > k * median - stdev) {
                        Peak peak = { 0, 0, 0, 0, 0 };
                        res.push_back(peak);
                        res.back().s = j;
                        res.back().s0 = i;
                        res.back().baseline = median;
                        flag = true;
                        ps = i;
                        //meanpre = mean;
                        //dpre = d;
                        break;
                    }
                    medianwindow.pop_back();
                    median += (median - data[j]) / medianwindow.size();
                }
            }
        }
        else {
            z = abs((res.back().baseline - mean) / (sqrt(d * d / window + stdev * stdev / window)));
            //k * mean > k * res.back().baseline - 0.5 * d && d < 1.2 * stdev && min > k * res.back().baseline - 3 * d && max > res.back().baseline + d
            if (z < 1.96 && min > k * res.back().baseline - 0.8 * th && max > k * res.back().baseline + stdev) {
                if (localmin == 0) {
                    flag = false;
                    res.pop_back();
                    pe = i;
                    continue;
                }
                res.back().e = localmin;
                std::vector<int> l = extreminval(data, ps, localmin, k);
                if (l.size() == 1)
                    res.back().currentpeak = data[l[0]];
                else if (l.size() > 1)
                    res.back().currentpeak = meancurrent(data, ps, localmin);
                if (abs(res.back().currentpeak - res.back().baseline) < 3 * stdev)
                    res.pop_back();
                flag = false;
                pe = i;
                localmin = 0;
            }
            else if (k * data[i] <= k * data[i - 1] && k * data[i] <= k * data[i + 1] && k * data[i] < k * res.back().baseline - 0.8 * th) {
                localmin = i;
            }
        }
    }
    if (flag)
        res.pop_back();
    return res;
}

std::list<Peak> findPeak_shortevent(std::vector<float>& data, int s, int e, int resolution, int window, int sigma, int k) {
    //stddev calculate
    float stdev = 0;
    calsd(data, s, e, stdev);
    float th = sigma * stdev;

    std::list<Peak> res;
    int localmin = 0;
    bool flag = false;
    std::deque<std::pair<float, int>> quemin;
    std::deque<std::pair<float, int>> quemax;
    std::deque<int> medianwindow;
    float median = 0;
    int ps = s;
    int pe = s;
    float min = 0;
    float max = 0;


    for (int i = s; i < e - 1; i++) {
        //cal moving mean for find start of peak
        if (medianwindow.size() < resolution) {
            medianwindow.push_back(i);
            median += (data[i] - median) / medianwindow.size();
        }
        else {
            median += (data[i] - data[medianwindow.front()]) / resolution;
            medianwindow.pop_front();
            medianwindow.push_back(i);
        }

        //maintain min by dequene
        while (!quemin.empty() && quemin.back().first >= k * data[i])
            quemin.pop_back();
        quemin.push_back(std::pair<float, int>(k * data[i], i));
        while (quemin.front().second <= i - window)
            quemin.pop_front();
        //time1 += double(clock() - times)/ CLOCKS_PER_SEC;
        min = quemin.front().first;
        while (!quemax.empty() && quemax.back().first <= k * data[i])
            quemax.pop_back();
        quemax.push_back(std::pair<float, int>(k * data[i], i));
        while (quemax.front().second <= i - window)
            quemax.pop_front();
        max = quemax.front().first;

        //maintain sliding std and mean, mean, d, d2
        if (i - s <= window - 1) {
            continue;
        }
        //float z = 0;

        if (!flag) {
            if (k * data[i] < k * median - th) {
                for (int j = i; j > pe && j == medianwindow.back(); j--) {
                    if (k * data[j] >= k * data[j - 1] && k * data[j] >= k * data[j + 1] && k * data[j] > k * median - stdev) {
                        Peak peak = { 0, 0, 0, 0, 0 };
                        res.push_back(peak);
                        res.back().s = j;
                        res.back().s0 = i;
                        res.back().baseline = median;
                        flag = true;
                        ps = i;
                        //meanpre = mean;
                        //dpre = d;
                        break;
                    }
                    medianwindow.pop_back();
                    median += (median - data[j]) / medianwindow.size();
                }
            }
        }
        else {
            if (k * data[i] > k * median && min > k * median - 0.7 * th) {
                if (localmin == 0) {
                    flag = false;
                    res.pop_back();
                    pe = i;
                    continue;
                }
                res.back().e = localmin;
                std::vector<int> l = extreminval(data, ps, localmin, k);
                if (l.size() == 1)
                    res.back().currentpeak = data[l[0]];
                else if (l.size() > 1)
                    res.back().currentpeak = meancurrent(data, ps, localmin);
                if (abs(res.back().currentpeak - res.back().baseline) < 3 * stdev)
                    res.pop_back();
                flag = false;
                pe = i;
                localmin = 0;
            }
            else if (k * data[i] <= k * data[i - 1] && k * data[i] <= k * data[i + 1] && k * data[i] < k * res.back().baseline - 0.8 * th) {
                localmin = i;
            }
        }
    }
    if (flag)
        res.pop_back();
    return res;
}
std::vector<float> Iir_filter(std::vector<float>& data, int fs, int cutoff) {
    Iir::Butterworth::LowPass<5> f;
    f.setup(fs, cutoff); 
    std::vector<float> res(data.size());
    for (int i = 0; i < data.size(); i++) {
        res[i] = f.filter(data[i]);
    }
    return res;
}
/*
p::list CUSUM(std::vector<float>& data, double stepsize, double h, double stddev) {
    std::vector<float> cpos(data.size(), 0.0);
    std::vector<float> cneg(data.size(), 0.0);
    std::vector<float> gpos(data.size(), 0.0);
    std::vector<float> gneg(data.size(), 0.0);


}

void cumSum(double* x, int n, int s, int e, double stepsize, double h, double stdDev, std::vector<int>& v1, std::vector<double>& v2) {
    int size = e - s + 1;
    double logp, logn = 0;
    double* cpos = new double[size];
    double* cneg = new double[size];
    double* gpos = new double[size];
    double* gneg = new double[size];
    for (int i = s; i <= e; i++) {
        cpos[i - s] = 0;
        cneg[i - s] = 0;
        gpos[i - s] = 0;
        gneg[i - s] = 0;
    }
    double anchor = s;
    double mean = x[s];
    double variance = stdDev * stdDev;
    int nstates = 0;
    double varM = x[s];
    double varS = 0;
    int jump = s;
    v1.push_back(s);
    for (int i = s + 1; i <= e; i++) {
        double varoldM = varM; //calculate move variance and mean.from mosaic. 
        varM = varM + (x[i] - varM) / (i - anchor + 1);
        varS = varS + (x[i] - varoldM) * (x[i] - varM);
        variance = varS / (i - anchor + 1);
        mean = ((i - anchor) * mean + x[i]) / (i - anchor + 1);
        variance = (variance == 0) ? stdDev * stdDev : variance;

        logp = stepsize * stdDev / variance * (x[i] - mean - stepsize * stdDev / 2);
        logn = -1 * stepsize * stdDev / variance * (x[i] - mean + stepsize * stdDev / 2);
        cpos[i - s] = cpos[i - s - 1] + logp;
        cneg[i - s] = cneg[i - s - 1] + logn;
        gpos[i - s] = std::max(gpos[i - s - 1] + logp, 0.);
        gneg[i - s] = std::max(gneg[i - s - 1] + logn, 0.);
        if (gpos[i - s] >= h) {
            jump = s + argmin(cpos, anchor - s, i - s);
            v2.push_back(meancurrent(x, n, v1.back(), jump));
            v1.push_back(jump);
            nstates++;
            anchor = i;
            cpos[i - s] = 0;
            cneg[i - s] = 0;
            gpos[i - s] = 0;
            gneg[i - s] = 0;
            mean = x[i];
            varM = x[i];
        }
        else if (gneg[i - s] >= h) {
            jump = s + argmin(cneg, anchor - s, i - s);
            v2.push_back(meancurrent(x, n, v1.back(), jump));
            v1.push_back(jump);
            nstates++;
            anchor = i;
            cpos[i - s] = 0;
            cneg[i - s] = 0;
            gpos[i - s] = 0;
            gneg[i - s] = 0;
            mean = x[i];
            varM = x[i];
        }
    }

    v2.push_back(meancurrent(x, n, v1.back(), e));
    v1.push_back(e);

    delete[] cpos;
    delete[] cneg;
    delete[] gpos;
    delete[] gneg;
    return;
}

int argmin(double* x, int s, int e) {
    double minval = x[s];
    int k = s;
    for (int i = s + 1; i <= e; i++) {
        if (x[i] < minval) {
            minval = x[i];
            k = i;
        }
    }
    return k;
}
*/