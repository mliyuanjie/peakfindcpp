p::tuple findMultiPeak(np::ndarray& x_py, int t) {
    int size = x_py.shape(0);
    double* x = new double[size];
    double* x_g = new double[size];
    if (x_py.get_dtype().get_itemsize() == 4) {
        int strides = x_py.get_strides()[0];
        strides = (strides < 0) ? -1 : 1;
        float* tmp = reinterpret_cast<float*> (x_py.get_data());
        for (int i = 0; i < size; i++) {
            x[i] = (double)*(tmp + i * strides);
            if (i == 1)
                x_g[0] = x[1] - x[0];
            else if (i > 1)
                x_g[i - 1] = (x[i] - x[i - 2]) / 2;
        }
        x_g[size - 1] = x[size - 1] - x[size - 2];
    }
    else {
        double* tmp = reinterpret_cast<double*> (x_py.get_data());
        int strides = x_py.get_strides()[0];
        strides = (strides < 0) ? -1 : 1;
        x = new double[size];
        for (int i = 0; i < size; i++) {
            x[i] = *(tmp + i * strides);
            if (i == 1)
                x_g[0] = x[1] - x[0];
            else if (i > 1)
                x_g[i - 1] = (x[i] - x[i - 2]) / 2;
        }
        x_g[size - 1] = x[size - 1] - x[size - 2];
    }

    float sd = 0;
    float mean = 0;
    double sd_g = 0;
    stdDev(x_g, size, sd_g);
    stdDevBaseline(x, size, sd, mean);
    gsl_vector* xb = gsl_vector_alloc(size);
    
    

    double* b = xb->data;

    p::list l1;
    p::list l2;
    double th = sd * (double)t;
    bool flag = false;
    std::vector<std::pair<int, int>> res;
    for (int i = 0; i < size; i++) {
        if (x[i] < b[i] - th && !flag) {
            for (int j = i; j <= i && j < size - 1 && j >= 1; j++) {
                if (res.size() > 0 && j <= res.back().second + 200) {
                    flag = true;
                    break;
                }
                if (x[j] <= x[j - 1] && x[j] <= x[j + 1]) {
                    std::pair<int, int> peak = { 0, 0 };
                    res.push_back(peak);
                    res.back().first = j;
                    flag = true;
                    i = j;
                    break;
                }
            }
        }
        if (flag && x[i] >= b[i]) {
            for (int j = i; j >= res.back().first && j < size - 1 && j >= 1; j--) {
                if (x[j] <= x[j - 1] && x[j] <= x[j + 1] && x[j] < b[j] - 0.5 * th) {
                    if (j <= res.back().first)
                        res.pop_back();
                    else 
                        res.back().second = j;
                    flag = false;
                    break;
                }
            }
        }
    }
    p::list l_index;
    p::list l_current;
    for (int i = 0; i < res.size(); i++) {
        int ps = res[i].first;
        int pe = res[i].second;
        int dt = ((pe - ps) < 100) ? 100 : (pe - ps);
        ps = ps - dt;
        ps = (ps <= 0) ? 0 : ps;
        pe = pe + dt;
        pe = (pe >= size) ? size - 1 : pe;
        std::vector<int> tmp = extremeval(x_g, ps, pe, sd_g);
        if (tmp.size() <= 2)
            continue;
        else{
            p::list l1;
            p::list l2;
            l1.append(tmp[0]);
            for (int j = 1; j < tmp.size(); j++) {
                l2.append(meancurrent(x, size, tmp[j - 1], tmp[j]));
                l1.append(tmp[j]);
            }
        }
    }
    gsl_vector_free(xb);
    delete[] x;
    delete[] x_g;
    return p::make_tuple(l_index, l_current);
}

void stdDevBaseline(double* x, size_t size, float& stdDev, float& baseline) {
    size_t s = 0;
    size_t e;
    if (s + 100000 >= size) {
        e = size;
        s = (e >= 100000) ? e - 100000 : 0;
    }
    else
        e = s + 100000;
    gsl_histogram* h = gsl_histogram_alloc(175);
    double xmin = x[s];
    double xmax = x[s];
    for (size_t i = s + 1; i < e; i++) {
        xmin = (x[i] < xmin) ? x[i] : xmin;
        xmax = (x[i] > xmax) ? x[i] : xmax;
    }
    gsl_histogram_set_ranges_uniform(h, xmin, xmax);
    for (size_t i = s; i < e; i++)
        gsl_histogram_increment(h, x[i]);
    double maxbin = gsl_histogram_max_val(h);
    for (size_t i = 0; i < h->n; i++) {
        if (h->bin[i] >= maxbin * 0.1) {
            xmin = h->range[i];
            break;
        }
    }
    for (int i = h->n - 1; i >= 0; i--) {
        if (h->bin[i] >= maxbin * 0.1) {
            xmax = h->range[i];
            break;
        }
    }
    gsl_histogram_free(h);
    int n = 0;
    for (size_t i = s; i < e; i++) {
        if (x[i] >= xmin && x[i] <= xmax) {
            baseline += x[i];
            n++;
        }
    }
    baseline /= n;
    n = 0;
    for (size_t i = s; i < e; i++) {
        if (x[i] >= xmin && x[i] <= xmax) {
            stdDev += (x[i] - baseline) * (x[i] - baseline);
            n++;
        }
    }
    stdDev = sqrt(stdDev / n);
    return;
}