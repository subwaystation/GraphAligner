#include "Assemble.h"
#include "stream.hpp"
#include "vg.pb.h"

Path Path::Reverse() const
{
	Path result = *this;
	std::reverse(result.position.begin(), result.position.end());
	for (size_t i = 0; i < result.position.size(); i++)
	{
		result.position[i] = result.position[i].Reverse();
	}
	return result;
}

bool AlignmentMatchCompareLT(const Alignment& left, const Alignment& right)
{
	return left.alignmentLength * left.alignmentIdentity < right.alignmentLength * right.alignmentIdentity;
}

AlignmentMatchComparerLT::AlignmentMatchComparerLT(const std::vector<Alignment>& paths) :
paths(&paths)
{}

AlignmentMatchComparerLT::AlignmentMatchComparerLT() :
paths(nullptr)
{}

bool AlignmentMatchComparerLT::operator()(const Alignment& left, const Alignment& right) const
{
	return AlignmentMatchCompareLT(left, right);
}

bool AlignmentMatchComparerLT::operator()(size_t left, size_t right) const
{
	assert(paths != nullptr);
	return AlignmentMatchCompareLT(paths->at(left), paths->at(right));
}

bool AlignmentQualityCompareLT(const Alignment& left, const Alignment& right)
{
	return left.alignmentIdentity < right.alignmentIdentity;
}

AlignmentQualityComparerLT::AlignmentQualityComparerLT(const std::vector<Alignment>& paths) :
paths(&paths)
{}

AlignmentQualityComparerLT::AlignmentQualityComparerLT() :
paths(nullptr)
{}

bool AlignmentQualityComparerLT::operator()(const Alignment& left, const Alignment& right) const
{
	return AlignmentQualityCompareLT(left, right);
}

bool AlignmentQualityComparerLT::operator()(size_t left, size_t right) const
{
	assert(paths != nullptr);
	return AlignmentQualityCompareLT(paths->at(left), paths->at(right));
}

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right)
{
	return left.first < right.first || (left.first == right.first && left.second < right.second);
}

void WriteAlignment(std::ofstream& file, const Alignment& aln)
{
	write(file, (uint64_t)aln.alignmentID);
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

std::vector<Path> loadAlignmentsAsPaths(std::string fileName, size_t minLen, const std::unordered_map<int, size_t>& nodeSizes)
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
		size_t totalSize = 0;
		for (auto p : pair.second)
		{
			for (auto pos : p.second)
			{
				totalSize += nodeSizes.at(pos.id);
				result.back().position.emplace_back(pos.id, pos.end);
			}
		}
		if (result.back().position.size() == 0)
		{
			result.pop_back();
			continue;
		}
		if (totalSize < minLen) result.pop_back();
	}
	return result;
}

std::unordered_map<int, size_t> getNodeSizes(const GfaGraph& graph)
{
	std::unordered_map<int, size_t> result;
	for (auto pair : graph.nodes)
	{
		result[pair.first] = pair.second.size() - graph.edgeOverlap;
	}
	return result;
}
