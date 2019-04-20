#include "Assemble.h"
#include "stream.hpp"
#include "vg.pb.h"

Path Path::Reverse() const
{
	Path result = *this;
	std::reverse(result.position.begin(), result.position.end());
	std::reverse(result.nodeSize.begin(), result.nodeSize.end());
	for (size_t i = 0; i < result.position.size(); i++)
	{
		result.position[i] = result.position[i].Reverse();
	}
	result.calculateOccurrences();
	result.calculateCumulativePrefixLength();
	return result;
}

void Path::calculateCumulativePrefixLength()
{
	cumulativePrefixLength.resize(nodeSize.size()+1, 0);
	for (size_t i = 1; i < cumulativePrefixLength.size(); i++)
	{
		cumulativePrefixLength[i] = cumulativePrefixLength[i-1] + nodeSize[i-1];
	}
}

void Path::calculateOccurrences()
{
	occurrences.clear();
	for (size_t i = 0; i < position.size(); i++)
	{
		occurrences[position[i]].push_back(i);
	}
}

bool AlignmentQualityCompareLT(const Alignment& left, const Alignment& right)
{
	return left.alignmentLength * left.alignmentIdentity < right.alignmentLength * right.alignmentIdentity;
}

bool AlignmentQualityCompareGT(const Alignment& left, const Alignment& right)
{
	return left.alignmentLength * left.alignmentIdentity > right.alignmentLength * right.alignmentIdentity;
}

AlignmentComparerLT::AlignmentComparerLT(const std::vector<Alignment>& paths) :
paths(&paths)
{}

AlignmentComparerLT::AlignmentComparerLT() :
paths(nullptr)
{}

bool AlignmentComparerLT::operator()(const Alignment& left, const Alignment& right) const
{
	return AlignmentQualityCompareLT(left, right);
}

bool AlignmentComparerLT::operator()(size_t left, size_t right) const
{
	assert(paths != nullptr);
	return AlignmentQualityCompareLT(paths->at(left), paths->at(right));
}

AlignmentComparerGT::AlignmentComparerGT(const std::vector<Alignment>& paths) :
paths(&paths)
{}

AlignmentComparerGT::AlignmentComparerGT() :
paths(nullptr)
{}

bool AlignmentComparerGT::operator()(const Alignment& left, const Alignment& right) const
{
	return AlignmentQualityCompareLT(left, right);
}

bool AlignmentComparerGT::operator()(size_t left, size_t right) const
{
	assert(paths != nullptr);
	return AlignmentQualityCompareGT(paths->at(left), paths->at(right));
}

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right)
{
	return left.first < right.first || (left.first == right.first && left.second < right.second);
}

void WriteAlignment(std::ofstream& file, const Alignment& aln)
{
	write(file, (uint32_t)aln.leftPath);
	write(file, (uint32_t)aln.rightPath);
	write(file, (uint32_t)aln.leftStart);
	write(file, (uint32_t)aln.rightStart);
	write(file, (uint32_t)aln.leftEnd);
	write(file, (uint32_t)aln.rightEnd);
	write(file, (char)aln.rightReverse);
	write(file, (uint64_t)aln.alignmentLength);
	write(file, (double)aln.alignmentIdentity);
	write(file, (uint32_t)aln.alignedPairs.size());
	for (size_t i = 0; i < aln.alignedPairs.size(); i++)
	{
		write(file, (uint32_t)aln.alignedPairs[i].leftIndex);
		write(file, (uint32_t)aln.alignedPairs[i].rightIndex);
		write(file, (char)aln.alignedPairs[i].leftReverse);
		write(file, (char)aln.alignedPairs[i].rightReverse);
	}
}

std::vector<Path> loadAlignmentsAsPaths(std::string fileName)
{
	std::unordered_map<std::string, std::vector<std::pair<size_t, std::vector<NodePos>>>> alnsPerRead;

	{
		std::ifstream file { fileName, std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&alnsPerRead](vg::Alignment& g) {
			std::vector<NodePos> thisAln;
			for (int i = 0; i < g.path().mapping_size(); i++)
			{
				thisAln.emplace_back(g.path().mapping(i).position().node_id(), !g.path().mapping(i).position().is_reverse());
			}
			alnsPerRead[g.name()].emplace_back(g.query_position(), thisAln);
		};
		stream::for_each(file, lambda);
	}

	std::vector<Path> result;
	for (auto pair : alnsPerRead)
	{
		assert(pair.second.size() > 0);
		std::sort(pair.second.begin(), pair.second.end(), [](const std::pair<size_t, std::vector<NodePos>>& left, const std::pair<size_t, std::vector<NodePos>>& right) { return left.first < right.first; });
		result.emplace_back();
		result.back().name = pair.first;
		for (auto p : pair.second)
		{
			for (auto pos : p.second)
			{
				result.back().position.emplace_back(pos.id, pos.end);
			}
		}
		if (result.back().position.size() == 0)
		{
			result.pop_back();
			continue;
		}
		result.back().calculateOccurrences();
	}
	return result;
}

std::vector<Path> addNodeLengths(const std::vector<Path>& original, const GfaGraph& graph)
{
	std::vector<Path> result = original;
	for (size_t i = 0; i < original.size(); i++)
	{
		result[i].nodeSize.reserve(result[i].position.size());
		for (size_t j = 0; j < original[i].position.size(); j++)
		{
			result[i].nodeSize.push_back(graph.nodes.at(original[i].position[j].id).size() - graph.edgeOverlap);
		}
		result[i].calculateCumulativePrefixLength();
	}
	return result;
}

std::vector<Path> filterByLength(const std::vector<Path>& paths, size_t minLen)
{
	std::vector<Path> result;
	for (auto path : paths)
	{
		size_t len = 0;
		for (auto nodeSize : path.nodeSize)
		{
			len += nodeSize;
		}
		if (len >= minLen) result.push_back(path);
	}
	return result;
}
