#include <iostream>
#include "realtime.h" 


double meancurrent(std::vector<double>& data) {
	double mean = 0;
	if (data.size() == 0)
		return 0;
	for (int i = 0; i < data.size(); i++) {
		mean += (data[i] - mean) / (i + 1);
	}
	return mean;
}

FindPeakRealTime::FindPeakRealTime(double samplerate, double filters, double sigmas, int resolutions, int windows, int ks) {
	fs = samplerate;
	filter = filters;
	sigma = sigmas;
	windowsizestart = resolutions;
	windowsizeend = windows;
	k = ks;
	if (filter / fs > 0.5 || filter < 0) {
		_append = &FindPeakRealTime::append_nofilter;
	}
	else {
		_append = &FindPeakRealTime::append_filter;
		f.setup(fs, filter);
	}
	interval = 1 / (samplerate * 1000);
}

void FindPeakRealTime::reset() {
	localminindex = 0;
	localmaxindex = 0;
	median = 0;
	medianold = 0;
	sd = 0;
	stddev = 0;
	mean = 0;
	meanold = 0;
	d2 = 0;
	d = 0;
	z = 0;
	min = 0;
	max = 0;
	s = 0;
	e = 0;
	i = 0;
	flag = false;
	s0flag = false;
	quemin.clear();
	quemax.clear();
	peakcurrent.clear();
	windowend.clear();
	windowstart.clear();
	result.clear();
}

void FindPeakRealTime::flush() {
	localminindex = 0;
	localmaxindex = 0;
	median = 0;
	medianold = 0;
	sd = 0;
	stddev = 0;
	flag = false;
	s0flag = false;
	peakcurrent.clear();
	windowstart.clear();
	result.clear();
}

p::tuple FindPeakRealTime::append_filter(std::vector<float>& x, int skip) {
	std::vector<std::vector<double>> resy(2);
	std::vector<std::vector<double>> rese(4);
	int number = 0;
	resy[0].reserve(ceil(x.size() / skip) * 2);
	resy[1].reserve(ceil(x.size() / skip) * 2);
	double mincurrent = FLT_MAX;
	double maxcurrent = -FLT_MAX;
	std::string eventstring;
	if (i == 0) {
		data = f.filter(x[0]);
		datasuf = f.filter(x[1]);
		i = 2;
		mincurrent = (data < datasuf) ? data : datasuf;
		maxcurrent = (data < datasuf) ? datasuf : data;
	}
	
	while (i < x.size() + s) {
		datapre = data;
		data = datasuf;	
		datasuf = f.filter(x[i - s]);
		if ((i - s) % skip == 0 && (i - s) != 0) {
			resy[0].push_back((i - 1) * interval);
			resy[0].push_back((i - 1) * interval);
			resy[1].push_back(mincurrent);
			resy[1].push_back(maxcurrent);
			mincurrent = FLT_MAX;
			maxcurrent = -FLT_MAX;
		}
		mincurrent = (mincurrent < datasuf) ? mincurrent : datasuf;
		maxcurrent = (maxcurrent > datasuf) ? maxcurrent : datasuf;

		while (!quemin.empty() && quemin.back().first >= k * data)
			quemin.pop_back();
		quemin.push_back(std::pair<double, size_t>(k * data, i - 1));
		while (quemin.front().second + windowsizeend <= i - 1)
			quemin.pop_front();
		//time1 += double(clock() - times)/ CLOCKS_PER_SEC;
		min = quemin.front().first;
		while (!quemax.empty() && quemax.back().first <= k * data)
			quemax.pop_back();
		quemax.push_back(std::pair<double, size_t>(k * data, i - 1));
		while (quemax.front().second + windowsizeend <= i - 1)
			quemax.pop_front();
		max = quemax.front().first;

		if ((int)windowend.size() < windowsizeend) {
			windowend.push_back(data);
			meanold = mean;
			mean += (data - mean) / windowend.size();
			d2 += (data - meanold) * (data - mean);
		}
		else {
			meanold = mean;
			mean += (data - windowend.front()) / windowsizeend;
			d2 += (data - windowend.front()) * (data + windowend.front() - mean - meanold);
			d = sqrt(d2 / windowsizeend);
			windowend.pop_front();
			windowend.push_back(data);
		}

		if (!flag) {
			if (int(windowstart.size()) < windowsizestart) {
				medianold = median;
				windowstart.push_back(std::pair<double, size_t>(data, i - 1));
				median += (data - median) / windowstart.size();
				sd += (data - medianold) * (data - median);
				stddev = sqrt(sd / windowstart.size());
				th = stddev * sigma;
			}
			else {
				medianold = median;
				median += (data - windowstart.front().first) / windowsizestart;
				sd += (data - windowstart.front().first) * (data + windowstart.front().first - median - medianold);
				windowstart.pop_front();
				windowstart.push_back(std::pair<double, size_t>(data, i - 1));
				stddev = sqrt(sd / windowstart.size());
				th = stddev * sigma;
			}
			if (int(windowstart.size()) < windowsizestart / 2) {
				i += 1;
				continue;
			}	
		}
		if (!flag) {
			if (k * data < k * median - th) {
				if (localmaxindex < windowstart.front().second) {
					flush();
					i += 1;
					continue;
				}
				while (windowstart.size() > 0 && windowstart.back().second != localmaxindex) {
					medianold = median;
					median += (median - windowstart.back().first) / ((int)windowstart.size() - 1);
					sd += (median - windowstart.back().first) * (windowstart.back().first - medianold);
					windowstart.pop_back();
					stddev = sqrt(sd / (int)windowstart.size());
					th = stddev * sigma;
				}
				Peak peak = { 0,0,0,0,0 };
				result.push_back(peak);
				result.back().s = localmaxindex;
				result.back().baseline = median;
				flag = true;
				s0flag = true;
				localmaxindex = 0;
			}
			else if (k * data >= k * datapre && k * data >= k * datasuf) {
				localmaxindex = i - 1;
			}
		}
		else {
			z = abs((result.back().baseline - mean) / (sqrt(d * d / windowsizeend + stddev * stddev / windowsizeend)));
			peakcurrent.push_back(data);
			if (z < 1.96 && min > k * result.back().baseline - 0.8 * th && max > k * result.back().baseline + stddev) {
				if (localminindex == 0 || s0flag) {
					flush();
					i += 1;
					continue;
				}
				for (int j = i - 1; j > localminindex; j--)
					peakcurrent.pop_back();
				result.back().e = localminindex;
				result.back().currentpeak = meancurrent(peakcurrent);
				if (k * result.back().baseline - k * result.back().currentpeak < th) {
					localminindex = 0;
					flag = false;
					peakcurrent.clear();
					i += 1;
					continue;
				}
				rese[0].push_back(result.back().s * interval);
				rese[0].push_back(result.back().s * interval);
				rese[0].push_back(localminindex * interval);
				rese[0].push_back(localminindex * interval);
				rese[1].push_back(result.back().baseline);
				rese[1].push_back(result.back().currentpeak);
				rese[1].push_back(result.back().currentpeak);
				rese[1].push_back(result.back().baseline);
				eventstring += std::to_string(result.back().s) + ',' \
					+ std::to_string(result.back().e) + ',' \
					+ std::to_string(result.back().s * interval * 1000) + ','\
					+ std::to_string(result.back().e * interval * 1000) + ','\
					+ std::to_string(result.back().baseline) + ','\
					+ std::to_string(result.back().currentpeak) + ','\
					+ std::to_string(result.back().s0) + ','\
					+ std::to_string(stddev) + '\n';
				localminindex = 0;
				flag = false;
				peakcurrent.clear();
				number++;
				
			}
			else if (k * data <= k * datapre && k * data <= k * datasuf && k * data < k * result.back().baseline - th) {
				localminindex = i - 1;
				if (s0flag) {
					result.back().s0 = localminindex;
					s0flag = false;
				}
			}
			else if (k * mean > k * result.back().baseline + 3 * stddev) {
				flush();
			}
		}
		i += 1;
	}
	resy[0].push_back((i - 1) * interval);
	resy[0].push_back((i - 1) * interval);
	resy[1].push_back(mincurrent);
	resy[1].push_back(maxcurrent);

	p::tuple res = p::make_tuple(resy, rese, median, stddev, number, eventstring, flag);
	s += x.size();
	return res;
}

p::tuple FindPeakRealTime::append_nofilter(std::vector<float>& x, int skip) {
	std::vector<std::vector<double>> resy(2);
	std::vector<std::vector<double>> rese(4);
	int number = 0;
	resy[0].reserve(ceil(x.size() / skip) * 2);
	resy[1].reserve(ceil(x.size() / skip) * 2);
	double mincurrent = FLT_MAX;
	double maxcurrent = -FLT_MAX;
	std::string eventstring;
	if (i == 0) {
		data = x[0];
		datasuf = x[1];
		i = 2;
		mincurrent = (data < datasuf) ? data : datasuf;
		maxcurrent = (data < datasuf) ? datasuf : data;
	}

	while (i < x.size() + s) {
		datapre = data;
		data = datasuf;
		datasuf = x[i - s];
		if ((i - s) % skip == 0 && (i - s) != 0) {
			resy[0].push_back((i - 1) * interval);
			resy[0].push_back((i - 1) * interval);
			resy[1].push_back(mincurrent);
			resy[1].push_back(maxcurrent);
			mincurrent = FLT_MAX;
			maxcurrent = -FLT_MAX;
		}
		mincurrent = (mincurrent < datasuf) ? mincurrent : datasuf;
		maxcurrent = (maxcurrent > datasuf) ? maxcurrent : datasuf;

		while (!quemin.empty() && quemin.back().first >= k * data)
			quemin.pop_back();
		quemin.push_back(std::pair<double, size_t>(k * data, i - 1));
		while (quemin.front().second + windowsizeend <= i - 1)
			quemin.pop_front();
		//time1 += double(clock() - times)/ CLOCKS_PER_SEC;
		min = quemin.front().first;
		while (!quemax.empty() && quemax.back().first <= k * data)
			quemax.pop_back();
		quemax.push_back(std::pair<double, size_t>(k * data, i - 1));
		while (quemax.front().second + windowsizeend <= i - 1)
			quemax.pop_front();
		max = quemax.front().first;

		if ((int)windowend.size() < windowsizeend) {
			windowend.push_back(data);
			meanold = mean;
			mean += (data - mean) / windowend.size();
			d2 += (data - meanold) * (data - mean);
		}
		else {
			meanold = mean;
			mean += (data - windowend.front()) / windowsizeend;
			d2 += (data - windowend.front()) * (data + windowend.front() - mean - meanold);
			d = sqrt(d2 / windowsizeend);
			windowend.pop_front();
			windowend.push_back(data);
		}

		if (!flag) {
			if (int(windowstart.size()) < windowsizestart) {
				medianold = median;
				windowstart.push_back(std::pair<double, size_t>(data, i - 1));
				median += (data - median) / windowstart.size();
				sd += (data - medianold) * (data - median);
				stddev = sqrt(sd / windowstart.size());
				th = stddev * sigma;
			}
			else {
				medianold = median;
				median += (data - windowstart.front().first) / windowsizestart;
				sd += (data - windowstart.front().first) * (data + windowstart.front().first - median - medianold);
				windowstart.pop_front();
				windowstart.push_back(std::pair<double, size_t>(data, i - 1));
				stddev = sqrt(sd / windowstart.size());
				th = stddev * sigma;
			}
			if (int(windowstart.size()) < windowsizestart / 2)
				i += 1;
				continue;
		}
		if (!flag) {
			if (k * data < k * median - th) {
				if (localmaxindex < windowstart.front().second) {
					flush();
					i += 1;
					continue;
				}
				while (windowstart.size() > 0 && windowstart.back().second != localmaxindex) {
					medianold = median;
					median += (median - windowstart.back().first) / ((int)windowstart.size() - 1);
					sd += (median - windowstart.back().first) * (windowstart.back().first - medianold);
					windowstart.pop_back();
					stddev = sqrt(sd / (int)windowstart.size());
					th = stddev * sigma;
				}
				Peak peak = { 0,0,0,0,0 };
				result.push_back(peak);
				result.back().s = localmaxindex;
				result.back().baseline = median;
				flag = true;
				s0flag = true;
				localmaxindex = 0;
			}
			else if (k * data >= k * datapre && k * data >= k * datasuf) {
				localmaxindex = i - 1;
			}
		}
		else {
			z = abs((result.back().baseline - mean) / (sqrt(d * d / windowsizeend + stddev * stddev / windowsizeend)));
			peakcurrent.push_back(data);
			if (z < 1.96 && min > k * result.back().baseline - 0.8 * th && max > k * result.back().baseline + stddev) {
				if (localminindex == 0 || s0flag) {
					flush();
					i += 1;
					continue;
				}
				for (int j = i - 1; j > localminindex; j--)
					peakcurrent.pop_back();
				result.back().e = localminindex;
				result.back().currentpeak = meancurrent(peakcurrent);
				if (k * result.back().baseline - k * result.back().currentpeak < th) {
					localminindex = 0;
					flag = false;
					peakcurrent.clear();
					i += 1;
					continue;
				}
				rese[0].push_back(result.back().s * interval);
				rese[0].push_back(result.back().s * interval);
				rese[0].push_back(localminindex * interval);
				rese[0].push_back(localminindex * interval);
				rese[1].push_back(result.back().baseline);
				rese[1].push_back(result.back().currentpeak);
				rese[1].push_back(result.back().currentpeak);
				rese[1].push_back(result.back().baseline);
				eventstring += std::to_string(result.back().s) + ',' \
					+ std::to_string(result.back().e) + ',' \
					+ std::to_string(result.back().s * interval * 1000) + ','\
					+ std::to_string(result.back().e * interval * 1000) + ','\
					+ std::to_string(result.back().baseline) + ','\
					+ std::to_string(result.back().currentpeak) + ','\
					+ std::to_string(result.back().s0) + ','\
					+ std::to_string(stddev) + '\n';
				localminindex = 0;
				flag = false;
				number++;
				peakcurrent.clear();

			}
			else if (k * data <= k * datapre && k * data <= k * datasuf && k * data < k * result.back().baseline - th) {
				localminindex = i - 1;
				if (s0flag) {
					result.back().s0 = localminindex;
					s0flag = false;
				}
			}
			else if (k * mean > k * result.back().baseline + 3 * stddev) {
				flush();
			}
		}
		i += 1;
	}
	resy[0].push_back((i - 1) * interval);
	resy[0].push_back((i - 1) * interval);
	resy[1].push_back(mincurrent);
	resy[1].push_back(maxcurrent);

	p::tuple res = p::make_tuple(resy, rese, median, stddev, number, eventstring, flag);
	s += x.size();
	return res;
}