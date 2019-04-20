#include <fstream>
#include <unordered_map>
#include <vector>
#include "GfaGraph.h"

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right);

struct Path
{
	std::string name;
	std::vector<NodePos> position;
	std::vector<size_t> nodeSize;
	std::vector<size_t> cumulativePrefixLength;
	std::unordered_map<NodePos, std::vector<size_t>> occurrences;
	Path Reverse() const;
	void calculateCumulativePrefixLength();
	void calculateOccurrences();
};

struct AlignmentMatch
{
	size_t leftIndex;
	size_t rightIndex;
	bool leftReverse;
	bool rightReverse;
};

struct Alignment
{
	size_t leftPath;
	size_t rightPath;
	std::vector<AlignmentMatch> alignedPairs;
	size_t leftStart;
	size_t leftEnd;
	size_t rightStart;
	size_t rightEnd;
	bool rightReverse;
	size_t alignmentLength;
	double alignmentIdentity;
};

template <typename T>
T read(std::ifstream& file)
{
	T result;
	file.read((char*)&result, sizeof(T));
	return result;
}

template <typename T>
void write(std::ofstream& file, T val)
{
	file.write((char*)&val, sizeof(T));
}

void WriteAlignment(std::ofstream& file, const Alignment& aln);

template <typename F>
void StreamAlignments(std::string filename, F f)
{
	std::ifstream file { filename, std::ios::in | std::ios::binary };
	while (file.good())
	{
		Alignment result;
		result.leftPath = read<uint32_t>(file);
		result.rightPath = read<uint32_t>(file);
		result.leftStart = read<uint32_t>(file);
		result.rightStart = read<uint32_t>(file);
		result.leftEnd = read<uint32_t>(file);
		result.rightEnd = read<uint32_t>(file);
		result.rightReverse = read<char>(file);
		result.alignmentLength = read<uint64_t>(file);
		result.alignmentIdentity = read<double>(file);
		size_t numMatches = read<uint32_t>(file);
		result.alignedPairs.resize(numMatches);
		for (size_t i = 0; i < numMatches; i++)
		{
			result.alignedPairs[i].leftIndex = read<uint32_t>(file);
			result.alignedPairs[i].rightIndex = read<uint32_t>(file);
			result.alignedPairs[i].leftReverse = read<char>(file);
			result.alignedPairs[i].rightReverse = read<char>(file);
		}
		f(result);
	}
}

bool AlignmentQualityCompareLT(const Alignment& left, const Alignment& right);

struct AlignmentComparerLT
{
public:
	AlignmentComparerLT(const std::vector<Alignment>& paths);
	AlignmentComparerLT();
	bool operator()(const Alignment& left, const Alignment& right) const;
	bool operator()(size_t left, size_t right) const;
private:
	const std::vector<Alignment>* const paths;
};

bool AlignmentQualityCompareGT(const Alignment& left, const Alignment& right);

struct AlignmentComparerGT
{
public:
	AlignmentComparerGT(const std::vector<Alignment>& paths);
	AlignmentComparerGT();
	bool operator()(const Alignment& left, const Alignment& right) const;
	bool operator()(size_t left, size_t right) const;
private:
	const std::vector<Alignment>* const paths;
};

std::vector<Path> loadAlignmentsAsPaths(std::string fileName);
std::vector<Path> addNodeLengths(const std::vector<Path>& original, const GfaGraph& graph);
std::vector<Path> filterByLength(const std::vector<Path>& paths, size_t minLen);
