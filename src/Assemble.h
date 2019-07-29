#include <fstream>
#include <unordered_map>
#include <vector>
#include "GfaGraph.h"

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right);

struct Path
{
	std::string name;
	std::vector<NodePos> position;
	Path Reverse() const;
};

struct AlignmentMatch
{
	AlignmentMatch() = default;
	AlignmentMatch(const AlignmentMatch&) = default;
	AlignmentMatch(AlignmentMatch&&) = default;
	AlignmentMatch& operator=(const AlignmentMatch&) = default;
	AlignmentMatch& operator=(AlignmentMatch&&) = default;
	size_t leftIndex;
	size_t rightIndex;
	bool leftReverse;
	bool rightReverse;
};

struct Alignment
{
	size_t alignmentID;
	size_t leftPath;
	size_t rightPath;
	std::vector<AlignmentMatch> alignedPairs;
	size_t leftStart;
	size_t leftEnd;
	size_t rightStart;
	size_t rightEnd;
	bool rightReverse;
	size_t matches;
	size_t mismatches;
	size_t alignmentLength;
};

template <typename T>
T read(std::ifstream& file)
{
	T result;
	file >> result;
	// file.read((char*)&result, sizeof(T));
	return result;
}

template <typename T>
void write(std::ofstream& file, T val)
{
	file << val << " ";
	// file.write((char*)&val, sizeof(T));
}

void WriteAlignment(std::ofstream& file, const Alignment& aln);

template <typename F>
void StreamAlignments(std::string filename, F f)
{
	std::ifstream file { filename, std::ios::in | std::ios::binary };
	while (file.good())
	{
		Alignment result;
		result.alignmentID = read<uint64_t>(file);
		result.leftPath = read<uint32_t>(file);
		result.rightPath = read<uint32_t>(file);
		result.leftStart = read<uint32_t>(file);
		result.rightStart = read<uint32_t>(file);
		result.leftEnd = read<uint32_t>(file);
		result.rightEnd = read<uint32_t>(file);
		result.rightReverse = read<char>(file);
		result.alignmentLength = read<uint64_t>(file);
		result.matches = read<uint64_t>(file);
		result.mismatches = read<uint64_t>(file);
		size_t numMatches = read<uint32_t>(file);
		result.alignedPairs.resize(numMatches);
		for (size_t i = 0; i < numMatches; i++)
		{
			result.alignedPairs[i].leftIndex = read<uint32_t>(file);
			result.alignedPairs[i].rightIndex = read<uint32_t>(file);
			result.alignedPairs[i].leftReverse = read<char>(file);
			result.alignedPairs[i].rightReverse = read<char>(file);
		}
		if (!file.good()) break;
		f(result);
	}
}

// bool AlignmentMatchCompareLT(const Alignment& left, const Alignment& right);
// bool AlignmentIdentityCompareLT(const Alignment& left, const Alignment& right);

// struct AlignmentMatchComparerLT
// {
// public:
// 	AlignmentMatchComparerLT(const std::vector<Alignment>& paths);
// 	AlignmentMatchComparerLT();
// 	bool operator()(const Alignment& left, const Alignment& right) const;
// 	bool operator()(size_t left, size_t right) const;
// private:
// 	const std::vector<Alignment>* const paths;
// };

// struct AlignmentQualityComparerLT
// {
// public:
// 	AlignmentQualityComparerLT(const std::vector<Alignment>& paths);
// 	AlignmentQualityComparerLT();
// 	bool operator()(const Alignment& left, const Alignment& right) const;
// 	bool operator()(size_t left, size_t right) const;
// private:
// 	const std::vector<Alignment>* const paths;
// };

std::vector<Path> loadAlignmentsAsPaths(std::string fileName, size_t minLen, const std::unordered_map<int, size_t>& nodeSizes);
std::unordered_map<int, size_t> getNodeSizes(const GfaGraph& graph);
