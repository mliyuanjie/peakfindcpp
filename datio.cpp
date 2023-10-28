#include <algorithm>
#include <fstream>
#include <Shlobj.h>
#include "datio.h" 

std::string desktop_directory(bool path_w)
{
	if (path_w == true)
	{
		WCHAR path[MAX_PATH + 1];
		if (SHGetSpecialFolderPathW(HWND_DESKTOP, path, CSIDL_DESKTOPDIRECTORY, FALSE))
		{
			std::wstring ws(path);
			std::string str(ws.begin(), ws.end());
			return str;
		}
		else return NULL;
	}
}


DATIO::DATIO(const std::string& fn) {
	fndata = fn;
	file.open(fndata);
	size_t length = file.size();
	if (length == 0)
		return;
	n = length / sizeof(float);
	pos = (float*)file.data();

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis('a', 'z');
	for (int i = 0; i < 8; i++) {
		char c = static_cast<char>(dis(gen));
		id += c;
	}
	buildpyramid();
	//filepyramid.open("datapyramid" + std::to_string(id) + ".dat");
	//pospyramid = (float*)filepyramid.data();
	return;
}


DATIO::DATIO(std::vector<float>& x) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis('a', 'z');
	for (int i = 0; i < 8; i++) {
		char c = static_cast<char>(dis(gen));
		id += c;
	}
	fndata = desktop_directory(true) +"/_dtemp_" + id + ".dat";
	std::ofstream wf(fndata, std::ios::out | std::ios::trunc);
	std::ostream_iterator<float> output_iterator(wf);
	wf.write((char*)&x[0], sizeof(float) * x.size());
	wf.close();
	file.open(fndata);
	n = x.size();
	if (n == 0)
		return;
	pos = (float*)file.data();
	buildpyramid();
	return;
}

DATIO::~DATIO() {
	if(file.is_open())
		file.close();
	if(filepyramid.is_open())
		filepyramid.close();
	pos = nullptr;
	pospyramid = nullptr;
	layer.clear();
	if (boost::filesystem::exists(fntemp))
		remove(fntemp.c_str());
	std::string fntemp = desktop_directory(true) + "/_dtemp_" + id + ".dat";
	if (boost::filesystem::exists(fntemp))
		remove(fntemp.c_str());
}

void DATIO::buildpyramid() {
	if (n <= 6000) {
		return;
	}
	layer.clear();
	bool flag = false;
	float min = FLT_MAX;
	float max = -FLT_MAX;

	size_t e = 0;
	size_t ntmp = ceil(double(n) / 2);
	while (ntmp > 3000) {
		e += ntmp;
		ntmp = ceil(double(ntmp) / 2);
	}
	
	fntemp = desktop_directory(true) + "/_ptemp_" + id + ".dat";
	std::ofstream wf;
	wf.open(fntemp, std::ios::out | std::ios::trunc);
	wf.seekp(e * sizeof(float)-1);
	wf.write("0", 1);
	wf.close();

	filepyramid.open(fntemp);
	pospyramid = (float*)filepyramid.data();

	e = 0;
	ntmp = ceil(double(n) / 2);
	layer.push_back(0);
	while (ntmp > 3000) {
		if (e == 0) {
			for (size_t i = 0; i < ntmp; i += 2) {
				min = FLT_MAX;
				max = -FLT_MAX;
				for (int j = 0; j < 4 && 2 * i + j < n; j++) {
					if (pos[2 * i + j] < min) {
						min = pos[2 * i + j];
						flag = false;
					}
					if (pos[2 * i + j] > max) {
						max = pos[2 * i + j];
						flag = true;
					}
				}
				if (i + 1 == ntmp) {
					pospyramid[i] = min;
					break;
				}
				if (flag) {
					pospyramid[i] = min;
					pospyramid[i + 1] = max;
				}
				else {
					pospyramid[i] = max;
					pospyramid[i + 1] = min;
				}
			}
		}
		else {
			int epre = layer[layer.size() - 2];
			for (size_t i = e; i < ntmp + e; i += 2) {
				min = FLT_MAX;
				max = -FLT_MAX;
				for (size_t j = 0; j < 4 && 2 * (i - e) + j + epre < e; j++) {
					if (pospyramid[2 * i + j - 2 * e + epre] < min) {
						min = pospyramid[2 * i + j - 2 * e + epre];
						flag = false;
					}
					if (pospyramid[2 * i + j - 2 * e + epre] > max) {
						max = pospyramid[2 * i + j - 2 * e + epre];
						flag = true;
					}
				}
				if (i + 1 == ntmp + e) {
					pospyramid[i] = min;
					break;
				}
				if (flag) {
					pospyramid[i] = min;
					pospyramid[i + 1] = max;
				}
				else {
					pospyramid[i] = max;
					pospyramid[i + 1] = min;
				}
			}
		}
		e += ntmp;
		layer.push_back(e);
		ntmp = ceil(double(ntmp) / 2);
	}
	return;
}


std::vector<double> DATIO::data(size_t start, size_t end) {
	if (!file.is_open())
		return std::vector<double>();
	end = (end == 0 || end > n) ? n : end;
	std::vector<double> res(end - start);
	size_t j = 0;
	for (size_t i = start; i < end; i++) {
		res[j] = pos[i];
		j++;
	}
	return res;
}

std::vector<double> DATIO::datafig(size_t start, size_t end) {
	size_t s = start;
	size_t e = end;
	size_t ntotal = e - s;
	int numlayer = 0;
	while (ntotal > 6000) {
		s = floor(double(s) / 2);
		e = ceil(double(e) / 2);
		numlayer++;
		ntotal = e - s;
	}
	if (numlayer == 0) {
		return data(start, end);
	}
	s += layer[numlayer-1];
	e += layer[numlayer-1];
	if (e > layer[numlayer]) 
		e = layer[numlayer];	
	std::vector<double> res(e - s);
	size_t j = 0;
	for (size_t i = s; i < e; i++) {
		res[j] = pospyramid[i];
		j++;
	}
	return res;
}


