#include "cfunction.h"
#include "realtime.h"
#include "tools.h"


std::vector<float> test(const std::vector<float>& data, int window) {

    int n = data.size();
    int l = n / window;
    int size = window * l;

    bool flag = true;
    float max = -FLT_MAX;
    float min = FLT_MAX;
    int e = 0;
    std::vector<float> res;
    res.reserve(2 * l + 2);
    for (int i = 0; i < size; i += window) {
        e = i + window;
        for (int j = i; j < e; j++) {
            if (min > data[j]) {
                min = data[j];
                flag = true;
            }
            if (max < data[j]) {
                max = data[j];
                flag = false;
            }
        }
        if (flag) {
            res.push_back(max);
            res.push_back(min);
        }
        else {
            res.push_back(min);
            res.push_back(max);
        }
        min = FLT_MAX;
        max = -FLT_MAX;
    }
    for (int i = size; i < n; i++) {
        if (min > data[i]) {
            min = data[i];
            flag = true;
        }
        if (max < data[i]) {
            max = data[i];
            flag = false;
        }
    }
    if (min == FLT_MAX)
        return res;
    if (min==max) {
        res.push_back(max);
    }
    else if (flag) {
        res.push_back(max);
        res.push_back(min);
    }
    else {
        res.push_back(min);
        res.push_back(max);
    }
    return res;
}

p::tuple randomWalkDt(double n, double dt, double dx, double prob, double xmin, double xmax, double ymax, double radius, double kon, double koff) {
    const double PI = 3.141592653589793238463;
    double datax = 0;
    double datay = 0;
    int collision = 0;
    double reacttime = 0;
    double staytime = 0;
    double rdt = 0;
    int success = 0;
    bool reactstart = false;
    double xsign = 1;
    double ysign = 1;

    std::random_device dev;
    std::mt19937_64 rng(dev());

    std::uniform_real_distribution<double> dist6(0.0, 1.0);
    std::normal_distribution<double> dist2{ 0.0, 1.0 };
    std::exponential_distribution<double> dist8(koff);
    std::exponential_distribution<double> dist9(kon);

    double runtime = 0;
    double i = 0;

    while (i < n) {
        if (datay >= ymax) {
            success = 1;
            break;
        }
        if (datay < xmin) {
            success = -1;
            break;
        }
        xsign = (dist6(rng) < 0.5) ? 1 : -1;
        ysign = (dist6(rng) <= prob) ? 1 : -1;

        if (reactstart) {
            rdt = dist8(rng);
            reacttime += rdt;
            reactstart = false;
            datax = datax + xsign * dx * abs(dist2(rng));
            datay = datay + ysign * dx * abs(dist2(rng));
            if (datax < xmin)
                datax = xmin;
            else if (datax > xmax)
                datax = xmax;
            i += rdt;
        }
        else {
            datax = datax + xsign * dx * abs(dist2(rng));
            datay = datay + ysign * dx * abs(dist2(rng));
            if (datax < xmin) {
                collision++;
                datax = xmin;
            }
            else if (datax > xmax) {
                collision++;
                datax = xmax;
            }
            else if (datax <= xmin + radius || datax >= xmax - radius) {
                if (dist9(rng) < staytime) {
                    reactstart = true;
                    staytime = 0;
                }
                else
                    staytime += dt;
            }
            i += dt;
        }
    }
    p::tuple res = p::make_tuple(i * 1e3, reacttime * 1e3, collision, success);
    return res;
}

p::tuple randomWalk(int n, double dx, double dangle, double dipolefield, double xprob, double x0, double y0, double angle0, double xmin, double xmax, double ymax, double radius, double kon, double koff, int react0) {
    const double PI = 3.141592653589793238463;
    std::vector<double> datax(n);
    std::vector<double> datay(n);
    std::vector<double> dataz(n);
    std::vector<int> reaction(n);
    int collision = 0;
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::uniform_real_distribution<double> dist6(0.0, 1.0);
    std::uniform_real_distribution<double> dist8(0.0, 1e9);
    std::normal_distribution<double> dist2{0.0, 1.0};
    datax[0] = x0;
    datay[0] = y0;
    dataz[0] = angle0;
    int end = -1;
    int success = 0;
    double xsign = 1;
    double ysign = 1;
    double zsign = 1;
    double esign = 1;
    if (dipolefield < 0)
        esign = -1;
    double prob = xprob;
    double on_prob = kon;
    double offprob = koff * 1e9;
    double prob_rorate = 0.5;
    int reactstart = react0;

    dangle = dangle * 31.6227766;
    for (int i = 1; i < n; i++) {
        if (datay[i - 1] >= ymax) {
            end = i;
            success = 1;
            break;
        }
        if (datay[i - 1] < xmin) {
            end = i;
            success = -1;
            break;
        }
        xsign = (dist6(rng) < 0.5) ? 1 : -1;
        ysign = (dist6(rng) <= prob) ? 1 : -1;
        
        if (i % 1000 == 0) {
            prob_rorate = 1.0 / (1.0 + std::exp(dipolefield * 4.05115441e-10 * (std::cos(dataz[i - 1] - dangle) - std::cos(dataz[i - 1] + dangle))));
            if (dist6(rng) > prob_rorate)
                zsign = -1 * esign;
            else
                zsign = esign;
            dataz[i] = dataz[i - 1] + dangle * zsign * abs(dist2(rng));
        }
        else {
            dataz[i] = dataz[i - 1];
        }

        
        if (reactstart > 0) {
            if (dist8(rng) < offprob) {
                datax[i] = datax[i - 1];
                datay[i] = datay[i - 1];
                reaction[i] = 1;
            }
            else {
                reactstart = -1;
                reaction[i] = 0;
                datax[i] = datax[i - 1] + xsign * dx * abs(dist2(rng));
                datay[i] = datay[i - 1] + ysign * dx * abs(dist2(rng));
                if (datax[i] < xmin + radius)
                    datax[i] = xmin + radius;
                else if (datax[i] > xmax - radius)
                    datax[i] = xmax - radius;
            }
        }
        else {
            datax[i] = datax[i - 1] + xsign * dx * abs(dist2(rng));
            datay[i] = datay[i - 1] + ysign * dx * abs(dist2(rng));
            reaction[i] = 0;
            if (datax[i] < xmin + radius) {
                collision++;
                if (dist6(rng) < on_prob) {
                    reactstart = 1;
                    reaction[i] = 1;
                    datax[i] = xmin + radius;
                }
                else {
                    datax[i] = xmin + radius;
                }
            }
            else if (datax[i] > xmax - radius) {
                collision++;
                if (dist6(rng) < on_prob) {
                    reactstart = 1;
                    reaction[i] = 1;
                    datax[i] = xmax - radius;
                }
                else {
                    datax[i] = xmax - radius;
                }
            }
            
        }
             
    }
    if (end > 0) {
        datax.resize(end);
        datay.resize(end);
        reaction.resize(end);
        dataz.resize(end);
    }
    //std::vector<double> angle = randomAngleWalk(datax.size(), dangle, angle0, dipolefield);
    p::tuple res = p::make_tuple(datax, datay, dataz, reaction, collision, success, reactstart);
    return res;
}

std::vector<double> randomAngleWalk(int n, double dangle, double angle0, double dipolefield) {
    double prob = 0.5;
    std::vector<double> data(n);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<> dist6(1, 100);
    std::normal_distribution<> dist2{ 0.0, 1.0 };
    data[0] = angle0;
    double sign = 1;
    if (dipolefield == 0) {
        for (int i = 1; i < data.size(); i++) {
            if (dist6(rng) <= 50)
                sign = 1;
            else
                sign = -1 ;
            data[i] = data[i - 1] + dangle * sign;
        }
        return data;
    }
    double esign = 1;
    if (dipolefield < 0)
        esign = -1;
    for (int i = 1; i < data.size(); i++) {
        prob = 100.0 / (1.0 + std::exp(dipolefield * 4.05115441e-10 * (std::cos(data[i - 1] - dangle) - std::cos(data[i - 1] + dangle))));
        if ((int)dist6(rng) > (int)prob)
            sign = -1 * esign;
        else
            sign = esign;
        data[i] = data[i - 1] + dangle * sign * abs(dist2(rng));
    }
    return data;
}

void randomangle_thread(float* data, float* angle0, float* dipolefield, float* dangle, int m, int i, int end, int skips) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<> dist6(1, 100);
    std::normal_distribution<> dist2{ 0.0, 1.0 };
    skips = (skips < 1) ? 1 : skips;
    for (; i < end; i++) {
        data[i * m] = angle0[i];
        double prob = 0.5;
        double sign = 1;
        float data_pre = angle0[i];
        if (dipolefield[i] == 0) {
            for (int j = 1; j < m; j++) {
                for (int skip = 0; skip < skips; skip++) {
                    if (dist6(rng) <= 50)
                        sign = 1;
                    else
                        sign = -1;
                    data_pre = data_pre + dangle[i] * sign;
                }
                data[i * m + j] = data_pre;
            }
            continue;
        }
        double esign = 1;
        if (dipolefield[i] < 0)
            esign = -1;
        for (int j = 1; j < m; j++) {
            for (int skip = 0; skip < skips; skip++) {
                prob = 100.0 / (1.0 + std::exp(dipolefield[i] * 4.05115441e-10 * (std::cos(data_pre - dangle[i]) - std::cos(data_pre + dangle[i]))));
                if ((int)dist6(rng) > (int)prob)
                    sign = -1 * esign;
                else
                    sign = esign;
                data_pre = data_pre + dangle[i] * sign * abs(dist2(rng));
            }
            data[i * m + j] = data_pre;
        }
    }
    return;
}

void randomAngleWalkParallel(p::array_t<float>& data_array, p::array_t<float>& angle0_array, p::array_t<float>& dipolefield_array, p::array_t<float>& dangle_array, int n, int m, int skips) {
    int max_thread = std::thread::hardware_concurrency();
    max_thread = (max_thread < 1) ? 1 : max_thread;
    std::vector<std::thread> threads;
    p::buffer_info buf_info = data_array.request(true);
    float* data = static_cast<float*>(buf_info.ptr);
    p::buffer_info buf_info1 = angle0_array.request(false);
    float* angle0 = static_cast<float*>(buf_info1.ptr);
    p::buffer_info buf_info2 = dipolefield_array.request(false);
    float* dipolefield = static_cast<float*>(buf_info2.ptr);
    p::buffer_info buf_info3 = dangle_array.request(false);
    float* dangle = static_cast<float*>(buf_info3.ptr);
    if (n < max_thread) {
        for (int thread_id = 0; thread_id < n; thread_id++) {
            threads.emplace_back(randomangle_thread, data, angle0, dipolefield, dangle, m, thread_id, thread_id + 1, skips);
        }
    }
    int chunk_size = n / max_thread;
    for (int thread_id = 0; thread_id < max_thread; thread_id++) {
        int start = thread_id * chunk_size;
        int end = (thread_id == max_thread - 1) ? n : (thread_id + 1) * chunk_size;
        threads.emplace_back(randomangle_thread, data, angle0, dipolefield, dangle, m, start, end, skips);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    return;
}


template <typename T>
p::dict findevent(p::array_t<T> x_py, int resolution, int window, int sigma, bool direction, bool islongevent) {
    p::buffer_info info = x_py.request();
    if (info.ndim != 1)
        throw std::runtime_error("Input array must be 1D");
    std::vector<float> data(info.size);
    auto ptr_py = x_py.unchecked<1>();
    int n = info.size;
    for (int i = 0; i < n; i++) 
        data[i] = ptr_py(i);
    int k = direction ? 1 : -1;
    
    window = (n < window) ? n / 2 : window;
    resolution = (n < resolution) ? n / 4 : resolution;

    std::list<Peak> res_c;
    if (islongevent)
        res_c = findPeak_longevent(data, 0, n, resolution, window, sigma, k);
    else
        res_c = findPeak_shortevent(data, 0, n, resolution, window, sigma, k);
    data.clear();
    //convert to python object
    p::dict y;
    p::list p1, p2, p3, p4, p5;
    std::list<Peak>::iterator it;
    for (it = res_c.begin(); it != res_c.end(); ++it) {
        p3.append(it->baseline);
        p4.append(it->currentpeak);
        p1.append(it->s);
        p2.append(it->e);
        p5.append(it->s0);
    }
    y["start"] = p1;
    y["end"] = p2;
    y["I0(pA)"] = p3;
    y["I1(pA)"] = p4;
    y["begin"] = p5;
    return y;
}

template <typename T>
Downsampler<T>::Downsampler(p::array_t<T> x_py, int samplingrate) {
    p::buffer_info info = x_py.request();
    if (info.ndim != 1)
        throw std::runtime_error("Input array must be 1D");
    std::vector<float> data(info.size);
    auto ptr_py = x_py.unchecked<1>();
    int n = info.size;
    for (int i = 0; i < n; i++)
        data[i] = ptr_py(i);
    interval = 1000 / samplingrate;
    dat = new DATIO(data);
}

template <typename T>
Downsampler<T>::Downsampler(const std::string& fn, int samplingrate) {
    interval = 1000 / samplingrate;
    dat = new DATIO(fn);
}

template <typename T>
Downsampler<T>::~Downsampler() {
    if (dat != nullptr)
        delete dat;
}

template <typename T>
p::tuple Downsampler<T>::downsample(double start, double end) {
    int s = start * 1000 / interval;
    int e = end * 1000 / interval;
    if (e >= int(dat->n) || e < s)
        e = dat->n - 1;
    s = (s<int(dat->n) && s >= 0) ? s : 0;
    std::vector<double> y = dat->datafig(s, e);
    std::vector<double> x = std::vector<double>(y.size());
    int k = ceil(double(e - s) / double(y.size()));
    for (int i = 0; i < y.size(); i++) {
        x[i] = double(s + i * k) * interval / 1000;
    }
    return p::make_tuple(x, y);
}


PYBIND11_MODULE(cfunction, m)
{
    m.def("findEvent", findevent<double>, "find the peak from 1d double array",
        p::arg("data"),
        p::arg("resolution") = 2000,
        p::arg("window") = 100,
        p::arg("threshold") = 5,
        p::arg("positive") = false,
        p::arg("islongevent") = true);

    m.def("findEvent", findevent<float>, "find the peak from 1d float array",
        p::arg("data"),
        p::arg("resolution") = 2000,
        p::arg("window") = 100,
        p::arg("threshold") = 5,
        p::arg("positive") = false,
        p::arg("islongevent") = true);
    m.def("test", test, "test array");
    m.def("Iir_filter", Iir_filter, "filter with iir low pass");
    m.def("randomWalk", randomWalk, "simulate protein translocation with rotation");
    m.def("randomWalkDt", randomWalkDt, "simulate protein translocation with dt");
    m.def("randomAngleWalk", randomAngleWalk, "simulate protein rotation");
    m.def("randomAngleWalkParallel", randomAngleWalkParallel, "simulate protein rotation with parallel algorithm");
    p::class_<Downsampler<float>>(m, "Downsampler")
        .def(p::init<const std::string&, int>())
        .def(p::init<p::array_t<float>, int>())
        .def("downsample", &Downsampler<float>::downsample);

    p::class_<FindPeakRealTime>(m, "FindPeakRealTime")
        .def(p::init<double, double, double, int, int, int>())
        .def("append", &FindPeakRealTime::append)
        .def("reset", &FindPeakRealTime::reset)
        .def("flush", &FindPeakRealTime::flush);

}
