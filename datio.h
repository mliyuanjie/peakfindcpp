#ifndef DATIO_H
#define DATIO_H

#include <string>
#include <vector>
#include <cmath>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <random>

class DATIO {
public:
	DATIO(const std::string&);
	DATIO(std::vector<float>&);
	~DATIO();
	std::vector<double> data(size_t start, size_t end);
	std::vector<double> datafig(size_t start, size_t end);
	size_t n = 0;
private:
	boost::iostreams::mapped_file_source file;
	boost::iostreams::mapped_file filepyramid;
	float* pospyramid;
	std::vector<int> layer;
	std::string id;
	float* pos;
	void buildpyramid();
	std::string fndata;
	std::string fntemp;
};
#endif //DATIO_H