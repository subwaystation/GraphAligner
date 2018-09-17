#include <iostream>
#include <fstream>
#include "vg.pb.h"
#include "CommonUtils.h"
#include "GfaGraph.h"

const int UnalignedNodeId = -1;
const int MismatchScore = -6;
const int MatchScore = 1;

struct Overlap
{
	size_t leftread;
	bool leftReverse;
	size_t rightread;
	bool rightReverse;
	size_t overlapSize;
	double overlapIdentity;
};

std::vector<std::string> getCorrectedReads(const std::vector<vg::Alignment>& alns, const GfaGraph& graph)
{
	std::vector<std::string> result;
	for (auto aln : alns)
	{
		std::string current;
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			if (aln.path().mapping(i).position().node_id() == UnalignedNodeId)
			{
				current += aln.path().mapping(i).edit(0).sequence();
			}
			else
			{
				int nodeId = aln.path().mapping(i).position().node_id();
				bool reverse = aln.path().mapping(i).position().is_reverse();
				size_t offset = aln.path().mapping(i).position().offset();
				size_t len = aln.path().mapping(i).edit(0).from_length();
				if (i == 0) offset = 0;
				if (i == aln.path().mapping_size()-1) len = graph.nodes.at(nodeId).size() - offset;
				std::string add = graph.nodes.at(nodeId);
				if (reverse) add = CommonUtils::ReverseComplement(add);
				current += add.substr(offset, len);
			}
		}
		result.push_back(current);
	}
	return result;
}

void writeGraph(std::string filename, const std::vector<std::string>& reads, const std::vector<Overlap>& overlaps)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < reads.size(); i++)
	{
		file << "S\t" << i << "\t" << reads[i] << std::endl;
	}
	for (auto overlap : overlaps)
	{
		file << "L\t" << overlap.leftread << "\t" << (overlap.leftReverse ? "-" : "+") << "\t" << overlap.rightread << "\t" << (overlap.rightReverse ? "-" : "+") << "\t" << overlap.overlapSize << "M" << std::endl;
	}
}

std::vector<std::pair<NodePos, size_t>> getNodeposAln(const vg::Alignment& aln, const GfaGraph& graph, bool reverse)
{
	std::vector<std::pair<NodePos, size_t>> result;
	for (int i = 0; i < aln.path().mapping_size(); i++)
	{
		result.emplace_back();
		result.back().first.id = aln.path().mapping(i).position().node_id();
		result.back().first.end = aln.path().mapping(i).position().is_reverse();
		if (aln.path().mapping(i).position().node_id() == -1)
		{
			result.back().second = aln.path().mapping(i).edit(0).from_length();
		}
		else
		{
			result.back().second = graph.nodes.at(result.back().first.id).size() - graph.edgeOverlap;
		}
	}
	if (reverse)
	{
		std::reverse(result.begin(), result.end());
		for (size_t i = 0; i < result.size(); i++)
		{
			result[i].first = result[i].first.Reverse();
		}
	}
	return result;
}

size_t getOverlapSize(const std::vector<std::pair<NodePos, size_t>>& left, const std::vector<std::pair<NodePos, size_t>>& right)
{
	std::vector<std::vector<int>> matchscore;
	matchscore.resize(left.size()+1);
	for (size_t i = 0; i < matchscore.size(); i++)
	{
		matchscore[i].resize(right.size()+1, 0);
	}
	for (size_t j = 1; j < right.size(); j++)
	{
		int rightnodelen = right[j].second;
		matchscore[0][j] = matchscore[0][j-1] + rightnodelen * MismatchScore;
	}
	for (size_t i = 0; i < left.size(); i++)
	{
		for (size_t j = 0; j < right.size(); j++)
		{
			int leftnodelen = left[i].second;
			int rightnodelen = right[j].second;
			if (left[i].first == right[j].first && left[i].first.id != UnalignedNodeId)
			{
				assert(leftnodelen == rightnodelen);
				matchscore[i+1][j+1] = matchscore[i][j] + leftnodelen * MatchScore;
			}
			else
			{
				matchscore[i+1][j+1] = matchscore[i][j] + std::max(leftnodelen, rightnodelen) * MismatchScore;
			}
			matchscore[i+1][j+1] = std::max(matchscore[i+1][j+1], matchscore[i+1][j] + rightnodelen * MismatchScore);
			matchscore[i+1][j+1] = std::max(matchscore[i+1][j+1], matchscore[i][j+1] + leftnodelen * MismatchScore);
		}
	}
	for (size_t j = right.size() - 1; j > 0; j--)
	{
		if (matchscore.back()[j] > 0)
		{
			size_t x = left.size();
			size_t y = j;
			while (y != 0)
			{
				if (x == 0)
				{
					y = 0;
					break;
				}
				assert(x > 0);
				assert(y > 0);
				int leftnodelen = left[x-1].second;
				int rightnodelen = right[y-1].second;
				if (left[x-1].first == right[y-1].first && left[x-1].first.id != UnalignedNodeId)
				{
					assert(leftnodelen == rightnodelen);
					if (matchscore[x][y] == matchscore[x-1][y-1] + leftnodelen * MatchScore)
					{
						x--;
						y--;
						continue;
					}
				}
				else if (matchscore[x][y] == matchscore[x-1][y-1] + std::max(leftnodelen, rightnodelen) * MismatchScore)
				{
					x--;
					y--;
					continue;
				}
				else if (matchscore[x][y] == matchscore[x-1][y] + leftnodelen * MismatchScore)
				{
					x--;
					continue;
				}
				else if (matchscore[x][y] == matchscore[x][y-1] + rightnodelen * MismatchScore)
				{
					y--;
					continue;
				}
				assert(false);
			}
			size_t leftlen = 0;
			size_t rightlen = 0;
			for (size_t i = x; i < left.size(); i++)
			{
				leftlen += left[i].second;
			}
			for (size_t i = y; i <= j; i++)
			{
				rightlen += right[i].second;
			}
			return std::min(leftlen, rightlen);
		}
	}
	return 0;
}

std::vector<Overlap> calculateOverlaps(const GfaGraph& graph, const std::vector<vg::Alignment>& alns, const double overlapLenFraction)
{
	std::unordered_map<int, std::vector<size_t>> alnsCrossingNode;
	std::vector<std::vector<std::pair<NodePos, size_t>>> nodeposAlns;
	std::vector<std::vector<std::pair<NodePos, size_t>>> reverseNodeposAlns;
	for (size_t i = 0; i < alns.size(); i++)
	{
		for (int j = 0; j < alns[i].path().mapping_size(); j++)
		{
			if (alns[i].path().mapping(j).position().node_id() == UnalignedNodeId) continue;
			alnsCrossingNode[alns[i].path().mapping(j).position().node_id()].push_back(i);
		}
		nodeposAlns.push_back(getNodeposAln(alns[i], graph, false));
		reverseNodeposAlns.push_back(getNodeposAln(alns[i], graph, true));
	}
	std::vector<Overlap> result;
	for (size_t i = 0; i < alns.size(); i++)
	{
		if (i % 100 == 0) std::cerr << i << "/" << alns.size() << std::endl;
		size_t size = 0;
		std::map<size_t, size_t> possibleOverlaps;
		for (size_t j = 0; j < nodeposAlns[i].size(); j++)
		{
			size += nodeposAlns[i][j].second;
			if (nodeposAlns[i][j].first.id == UnalignedNodeId) continue;
			for (auto a : alnsCrossingNode[nodeposAlns[i][j].first.id])
			{
				if (a >= i) continue;
				possibleOverlaps[a] += nodeposAlns[i][j].second;
			}
		}
		size_t minOverlapSize = size * overlapLenFraction;
		for (auto pair : possibleOverlaps)
		{
			size_t a = pair.first;
			assert(a < i);
			if (pair.second < minOverlapSize) continue;
			auto overlap = getOverlapSize(nodeposAlns[i], nodeposAlns[a]);
			if (overlap >= minOverlapSize)
			{
				result.emplace_back();
				result.back().overlapSize = overlap;
				result.back().leftread = i;
				result.back().rightread = a;
				result.back().leftReverse = false;
				result.back().rightReverse = false;
			}
			overlap = getOverlapSize(nodeposAlns[i], reverseNodeposAlns[a]);
			if (overlap >= minOverlapSize)
			{
				result.emplace_back();
				result.back().overlapSize = overlap;
				result.back().leftread = i;
				result.back().rightread = a;
				result.back().leftReverse = false;
				result.back().rightReverse = true;
			}
			overlap = getOverlapSize(reverseNodeposAlns[i], nodeposAlns[a]);
			if (overlap >= minOverlapSize)
			{
				result.emplace_back();
				result.back().overlapSize = overlap;
				result.back().leftread = i;
				result.back().rightread = a;
				result.back().leftReverse = true;
				result.back().rightReverse = false;
			}
			overlap = getOverlapSize(reverseNodeposAlns[i], reverseNodeposAlns[a]);
			if (overlap >= minOverlapSize)
			{
				result.emplace_back();
				result.back().overlapSize = overlap;
				result.back().leftread = i;
				result.back().rightread = a;
				result.back().leftReverse = true;
				result.back().rightReverse = true;
			}
		}
	}
	return result;
}

std::vector<Overlap> filterOverlaps(const std::vector<Overlap>& raws, const std::vector<std::string>& reads, const double overlapLenFraction, const double secondaryLenFraction)
{
	std::unordered_map<NodePos, size_t> largestOverlap;
	for (auto overlap : raws)
	{
		NodePos left, right;
		left.id = overlap.leftread;
		left.end = overlap.leftReverse;
		right.id = overlap.rightread;
		right.end = overlap.rightReverse;
		largestOverlap[left] = std::max(largestOverlap[left], overlap.overlapSize);
		largestOverlap[right.Reverse()] = std::max(largestOverlap[right.Reverse()], overlap.overlapSize);
	}
	std::vector<Overlap> result;
	for (auto overlap : raws)
	{
		size_t leftlen = reads[overlap.leftread].size();
		size_t rightlen = reads[overlap.rightread].size();
		NodePos left, right;
		left.id = overlap.leftread;
		left.end = overlap.leftReverse;
		right.id = overlap.rightread;
		right.end = overlap.rightReverse;
		if (overlap.overlapSize < largestOverlap[left] * secondaryLenFraction) continue;
		if (overlap.overlapSize < largestOverlap[right.Reverse()] * secondaryLenFraction) continue;
		if (overlap.overlapSize < leftlen * overlapLenFraction) continue;
		if (overlap.overlapSize < rightlen * overlapLenFraction) continue;
		result.push_back(overlap);
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string inputAlns { argv[1] };
	std::string inputGraph { argv[2] };
	double overlapLenFraction = std::stod(argv[3]);
	double secondaryLenFraction = std::stod(argv[4]);
	std::string outputGraph { argv[5] };

	auto graph = GfaGraph::LoadFromFile(inputGraph);
	auto alns = CommonUtils::LoadVGAlignments(inputAlns);

	auto overlaps = calculateOverlaps(graph, alns, overlapLenFraction);
	auto corrected = getCorrectedReads(alns, graph);
	auto filteredOverlaps = filterOverlaps(overlaps, corrected, overlapLenFraction, secondaryLenFraction);
	
	writeGraph(outputGraph, corrected, filteredOverlaps);
}