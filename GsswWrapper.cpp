#include <iostream>
#include <unistd.h>
#include <fstream>
#include "Aligner.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
	std::string graphFile = "";
	std::string fastqFile = "";
	std::string seedFile = "";
	std::string alignmentFile = "";
	int numThreads = 0;
	int bandwidth = 0;
	double bandwidthStddev = 0;
	int c;

	while ((c = getopt(argc, argv, "g:f:s:a:t:b:B:")) != -1)
	{
		switch(c)
		{
			case 'g':
				graphFile = std::string(optarg);
				break;
			case 'f':
				fastqFile = std::string(optarg);
				break;
			case 's':
				seedFile = std::string(optarg);
				break;
			case 'a':
				alignmentFile = std::string(optarg);
				break;
			case 't':
				numThreads = std::stoi(optarg);
				break;
			case 'b':
				bandwidth = std::stoi(optarg);
				break;
			case 'B':
				bandwidthStddev = std::stod(optarg);
				break;
		}
	}

	if (numThreads < 1)
	{
		std::cerr << "select number of threads -t" << std::endl;
		std::exit(0);
	}

	if (bandwidth < 2 && bandwidthStddev < 0.1)
	{
		std::cerr << "select either bandwidth -b or bandwidth standard deviation -B" << std::endl;
		std::exit(0);
	}

	alignReads(graphFile, fastqFile, seedFile, numThreads, bandwidth, bandwidthStddev, alignmentFile);

	return 0;
}
