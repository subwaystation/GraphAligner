#include <queue>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include "stream.hpp"
#include "CommonUtils.h"
#include "vg.pb.h"
#include "GfaGraph.h"
#include "Assemble.h"

template <typename T>
class Oriented2dVector
{
public:
	T& operator[](const std::pair<size_t, NodePos>& pos)
	{
		if (pos.second.end) return plus[pos.first][pos.second.id];
		return minus[pos.first][pos.second.id];
	}
	const T& operator[](const std::pair<size_t, NodePos>& pos) const
	{
		if (pos.second.end) return plus[pos.first][pos.second.id];
		return minus[pos.first][pos.second.id];
	}
	void resize(size_t size)
	{
		plus.resize(size);
		minus.resize(size);
	}
	void resizePart(size_t index, size_t size)
	{
		plus[index].resize(size);
		minus[index].resize(size);
	}
private:
	std::vector<std::vector<T>> plus;
	std::vector<std::vector<T>> minus;
};

struct TransitiveClosureMapping
{
	std::map<std::pair<size_t, NodePos>, size_t> mapping;
};

std::pair<NodePos, NodePos> canon(NodePos left, NodePos right)
{
	if (left.id == right.id)
	{
		if (!left.end && !right.end) return std::make_pair(right.Reverse(), left.Reverse());
		return std::make_pair(left, right);
	}
	if (left < right) return std::make_pair(left, right);
	assert(right.Reverse() < left.Reverse());
	return std::make_pair(right.Reverse(), left.Reverse());
}

struct ClosureEdges
{
	std::map<std::pair<NodePos, NodePos>, size_t> coverage;
	std::map<std::pair<NodePos, NodePos>, size_t> overlap;
};

struct DoublestrandedTransitiveClosureMapping
{
	std::map<std::pair<size_t, size_t>, NodePos> mapping;
};

template <typename T>
T find(std::map<T, T>& parent, T key)
{
	if (parent.count(key) == 0)
	{
		parent[key] = key;
		return key;
	}
	if (parent.at(key) == key)
	{
		return key;
	}
	auto result = find(parent, parent.at(key));
	parent[key] = result;
	return result;
}

template <typename T>
void set(std::map<T, T>& parent, T key, T target)
{
	auto found = find(parent, key);
	parent[found] = find(parent, target);
}

std::pair<size_t, NodePos> find(Oriented2dVector<std::pair<size_t, NodePos>>& parent, std::pair<size_t, NodePos> key)
{
	if (parent[key] == key)
	{
		return key;
	}
	auto result = find(parent, parent[key]);
	parent[key] = result;
	return result;
}

void set(Oriented2dVector<std::pair<size_t, NodePos>>& parent, std::pair<size_t, NodePos> key, std::pair<size_t, NodePos> target)
{
	auto found = find(parent, key);
	parent[found] = find(parent, target);
}

TransitiveClosureMapping getTransitiveClosures(const std::vector<Path>& paths, const std::set<std::pair<size_t, size_t>>& pickedAlns, std::string overlapFile)
{
	TransitiveClosureMapping result;
	Oriented2dVector<std::pair<size_t, NodePos>> parent;
	parent.resize(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		parent.resizePart(i, paths[i].position.size());
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			parent[std::pair<size_t, NodePos> { i, NodePos { j, true } }] = std::pair<size_t, NodePos> { i, NodePos { j, true } };
			parent[std::pair<size_t, NodePos> { i, NodePos { j, false } }] = std::pair<size_t, NodePos> { i, NodePos { j, false } };
		}
	}
	{
		StreamAlignments(overlapFile, [&parent, &pickedAlns](const Alignment& aln){
			bool picked = pickedAlns.count(std::pair<size_t, size_t>{aln.leftPath, aln.rightPath}) == 1;
			if (!picked) return;
			for (auto match : aln.alignedPairs)
			{
				std::pair<size_t, NodePos> leftKey { aln.leftPath, NodePos { match.leftIndex, match.leftReverse } };
				std::pair<size_t, NodePos> rightKey { aln.rightPath, NodePos { match.rightIndex, match.rightReverse } };
				set(parent, leftKey, rightKey);
				std::pair<size_t, NodePos> revLeftKey { aln.leftPath, NodePos { match.leftIndex, !match.leftReverse } };
				std::pair<size_t, NodePos> revRightKey { aln.rightPath, NodePos { match.rightIndex, !match.rightReverse } };
				set(parent, revLeftKey, revRightKey);
			}
		});
	}
	std::map<std::pair<size_t, NodePos>, size_t> closureNumber;
	size_t nextClosure = 1;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			auto key = std::pair<size_t, NodePos> { i, NodePos { j, true } };
			auto found = find(parent, key);
			if (closureNumber.count(found) == 0)
			{
				closureNumber[found] = nextClosure;
				nextClosure += 1;
			}
			result.mapping[key] = closureNumber.at(found);
			key = std::pair<size_t, NodePos> { i, NodePos { j, false } };
			found = find(parent, key);
			if (closureNumber.count(found) == 0)
			{
				closureNumber[found] = nextClosure;
				nextClosure += 1;
			}
			result.mapping[key] = closureNumber.at(found);
		}
	}
	std::cerr << (nextClosure-1) << " transitive closure sets" << std::endl;
	std::cerr << result.mapping.size() << " transitive closure items" << std::endl;
	return result;
}

GfaGraph getGraph(const DoublestrandedTransitiveClosureMapping& transitiveClosures, const ClosureEdges& edges, const std::vector<Path>& paths, const GfaGraph& graph)
{
	std::unordered_map<size_t, size_t> closureCoverage;
	for (auto pair : transitiveClosures.mapping)
	{
		closureCoverage[pair.second.id] += 1;
	}
	std::unordered_set<size_t> outputtedClosures;
	GfaGraph result;
	result.edgeOverlap = graph.edgeOverlap;
	for (auto pair : transitiveClosures.mapping)
	{
		if (outputtedClosures.count(pair.second.id) == 1) continue;
		NodePos pos = paths[pair.first.first].position[pair.first.second];
		auto seq = graph.nodes.at(pos.id);
		if (!pos.end) seq = CommonUtils::ReverseComplement(seq);
		if (!pair.second.end) seq = CommonUtils::ReverseComplement(seq);
		result.nodes[pair.second.id] = seq;
		result.tags[pair.second.id] = "LN:i:" + std::to_string(seq.size() - graph.edgeOverlap) + "\tRC:i:" + std::to_string((seq.size() - graph.edgeOverlap) * closureCoverage[pair.second.id]) + "\tkm:f:" + std::to_string(closureCoverage[pair.second.id]) + "\toi:Z:" + std::to_string(pos.id) + (pos.end ? "+" : "-");
		outputtedClosures.insert(pair.second.id);
	}
	std::cerr << outputtedClosures.size() << " outputted closures" << std::endl;
	for (auto pair : edges.coverage)
	{
		if (outputtedClosures.count(pair.first.first.id) == 0 || outputtedClosures.count(pair.first.second.id) == 0) continue;
		result.edges[pair.first.first].push_back(pair.first.second);
		result.edgeTags[std::make_pair(pair.first.first, pair.first.second)] = "RC:i:" + std::to_string(pair.second);
	}
	result.varyingOverlaps.insert(edges.overlap.begin(), edges.overlap.end());
	std::cerr << edges.coverage.size() << " outputted edges" << std::endl;
	return result;
}

DoublestrandedTransitiveClosureMapping mergeDoublestrandClosures(const std::vector<Path>& paths, const TransitiveClosureMapping& original)
{
	DoublestrandedTransitiveClosureMapping result;
	std::unordered_map<size_t, NodePos> mapping;
	int nextId = 1;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			std::pair<size_t, NodePos> fwKey { i, NodePos { j, true } };
			std::pair<size_t, NodePos> bwKey { i, NodePos { j, false } };
			assert(original.mapping.count(fwKey) == 1);
			assert(original.mapping.count(bwKey) == 1);
			size_t fwSet = original.mapping.at(fwKey);
			size_t bwSet = original.mapping.at(bwKey);
			assert(mapping.count(fwSet) == mapping.count(bwSet));
			if (mapping.count(fwSet) == 0)
			{
				if (fwSet == bwSet)
				{
					mapping[fwSet] = NodePos { nextId, true };
					assert(false);
				}
				else
				{
					mapping[fwSet] = NodePos { nextId, true };
					mapping[bwSet] = NodePos { nextId, false };
				}
				nextId += 1;
			}
			assert(mapping.count(fwSet) == 1);
			result.mapping[std::pair<size_t, size_t> { i, j }] = mapping.at(fwSet);
		}
	}
	std::cerr << (nextId-1) << " doublestranded transitive closure sets" << std::endl;
	return result;
}

std::vector<Alignment> doubleAlignments(const std::vector<Alignment>& alns)
{
	std::vector<Alignment> result = alns;
	result.reserve(alns.size() * 2);
	for (auto aln : alns)
	{
		result.emplace_back();
		result.back().alignmentLength = aln.alignmentLength;
		result.back().alignmentIdentity = aln.alignmentIdentity;
		result.back().leftPath = aln.leftPath;
		result.back().rightPath = aln.rightPath;
		result.back().alignedPairs = aln.alignedPairs;
		for (size_t i = 0; i < result.back().alignedPairs.size(); i++)
		{
			result.back().alignedPairs[i].leftReverse = !result.back().alignedPairs[i].leftReverse;
			result.back().alignedPairs[i].rightReverse = !result.back().alignedPairs[i].rightReverse;
		}
	}
	std::cerr << result.size() << " alignments after doubling" << std::endl;
	return result;
}

std::vector<Alignment> removeContained(const std::vector<Path>& paths, const std::vector<Alignment>& original)
{
	std::vector<std::vector<size_t>> continuousEnd;
	continuousEnd.resize(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		continuousEnd[i].resize(paths[i].position.size(), 0);
	}
	for (auto aln : original)
	{
		assert(aln.leftStart <= aln.leftEnd);
		assert(aln.rightStart <= aln.rightEnd);
		for (size_t i = aln.leftStart; i <= aln.leftEnd; i++)
		{
			continuousEnd[aln.leftPath][i] = std::max(continuousEnd[aln.leftPath][i], aln.leftEnd);
		}
		for (size_t i = aln.rightStart; i <= aln.rightEnd; i++)
		{
			continuousEnd[aln.rightPath][i] = std::max(continuousEnd[aln.rightPath][i], aln.rightEnd);
		}
	}
	std::vector<Alignment> result;
	for (auto aln : original)
	{
		if (continuousEnd[aln.leftPath][aln.leftStart] > aln.leftEnd) continue;
		if (aln.leftStart > 0 && continuousEnd[aln.leftPath][aln.leftStart-1] >= aln.leftEnd) continue;
		if (continuousEnd[aln.rightPath][aln.rightStart] > aln.rightEnd) continue;
		if (aln.rightStart > 0 && continuousEnd[aln.rightPath][aln.rightStart-1] >= aln.rightEnd) continue;
		result.push_back(aln);
	}
	std::cerr << result.size() << " alignments after removing contained" << std::endl;
	return result;
}

std::vector<Alignment> pickLowestErrorPerRead(const std::vector<Path>& paths, const std::vector<Alignment>& alns, size_t maxNum)
{
	std::vector<std::vector<Alignment>> alnsPerRead;
	alnsPerRead.resize(paths.size());
	for (auto aln : alns)
	{
		alnsPerRead[aln.leftPath].push_back(aln);
		alnsPerRead[aln.rightPath].push_back(aln);
	}
	std::vector<Alignment> result;
	for (size_t i = 0; i < alnsPerRead.size(); i++)
	{
		if (alnsPerRead[i].size() > maxNum)
		{
			std::sort(alnsPerRead[i].begin(), alnsPerRead[i].end(), [](const Alignment& left, const Alignment& right){ return left.alignmentIdentity < right.alignmentIdentity; });
			for (size_t j = alnsPerRead[i].size() - maxNum; j < alnsPerRead[i].size(); j++)
			{
				result.push_back(alnsPerRead[i][j]);
			}
		}
		else
		{
			result.insert(result.end(), alnsPerRead[i].begin(), alnsPerRead[i].end());
		}
	}
	std::cerr << result.size() << " alignments after picking lowest erro" << std::endl;
	return result;
}

std::set<std::pair<size_t, size_t>> pickLongestPerRead(const std::vector<Path>& paths, std::string alnFile, size_t maxNum)
{
	std::vector<Alignment> alns;
	StreamAlignments(alnFile, [&alns](const Alignment& aln){
		alns.emplace_back(aln);
		decltype(aln.alignedPairs) tmp;
		std::swap(tmp, alns.back().alignedPairs);
	});
	std::vector<std::vector<size_t>> leftAlnsPerRead;
	std::vector<std::vector<size_t>> rightAlnsPerRead;
	std::vector<int> picked;
	picked.resize(alns.size(), 0);
	leftAlnsPerRead.resize(paths.size());
	rightAlnsPerRead.resize(paths.size());
	for (size_t i = 0; i < alns.size(); i++)
	{
		assert(alns[i].leftPath < paths.size());
		assert(alns[i].rightPath < paths.size());
		assert(alns[i].leftEnd < paths[alns[i].leftPath].position.size());
		assert(alns[i].rightEnd < paths[alns[i].rightPath].position.size());
		if (alns[i].leftStart == 0) leftAlnsPerRead[alns[i].leftPath].push_back(i);
		if (alns[i].leftEnd == paths[alns[i].leftPath].position.size()-1) rightAlnsPerRead[alns[i].leftPath].push_back(i);
		if (alns[i].rightStart == 0) leftAlnsPerRead[alns[i].rightPath].push_back(i);
		if (alns[i].rightEnd == paths[alns[i].rightPath].position.size()-1) rightAlnsPerRead[alns[i].rightPath].push_back(i);
	}
	for (size_t i = 0; i < leftAlnsPerRead.size(); i++)
	{
		std::sort(leftAlnsPerRead[i].begin(), leftAlnsPerRead[i].end(), AlignmentComparerLT { alns });
		std::sort(rightAlnsPerRead[i].begin(), rightAlnsPerRead[i].end(), AlignmentComparerLT { alns });
		std::set<size_t> pickedHere;
		for (size_t j = leftAlnsPerRead[i].size() > maxNum ? (leftAlnsPerRead[i].size() - maxNum) : 0; j < leftAlnsPerRead[i].size(); j++)
		{
			pickedHere.insert(leftAlnsPerRead[i][j]);
		}
		for (size_t j = rightAlnsPerRead[i].size() > maxNum ? (rightAlnsPerRead[i].size() - maxNum) : 0; j < rightAlnsPerRead[i].size(); j++)
		{
			pickedHere.insert(rightAlnsPerRead[i][j]);
		}
		for (auto index : pickedHere)
		{
			picked[index] += 1;
		}
	}
	std::set<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < alns.size(); i++)
	{
		assert(picked[i] >= 0);
		assert(picked[i] <= 2);
		if (picked[i] == 2) result.emplace(alns[i].leftPath, alns[i].rightPath);
	}
	std::cerr << result.size() << " alignments after picking longest" << std::endl;
	return result;
}

DoublestrandedTransitiveClosureMapping removeOutsideCoverageClosures(const DoublestrandedTransitiveClosureMapping& closures, int minCoverage, int maxCoverage)
{
	std::unordered_map<size_t, size_t> coverage;
	for (auto pair : closures.mapping)
	{
		coverage[pair.second.id] += 1;
	}
	DoublestrandedTransitiveClosureMapping result;
	std::unordered_set<size_t> numbers;
	for (auto pair : closures.mapping)
	{
		if (coverage[pair.second.id] >= minCoverage && coverage[pair.second.id] <= maxCoverage)
		{
			result.mapping[pair.first] = pair.second;
			numbers.insert(pair.second.id);
		}
	}
	std::cerr << numbers.size() << " closures after removing low coverage" << std::endl;
	std::cerr << result.mapping.size() << " closure items after removing low coverage" << std::endl;
	return result;
}

DoublestrandedTransitiveClosureMapping insertMiddles(const DoublestrandedTransitiveClosureMapping& original, const std::vector<Path>& paths)
{
	DoublestrandedTransitiveClosureMapping result;
	result = original;
	int nextNum = 0;
	for (auto pair : original.mapping)
	{
		nextNum = std::max(nextNum, pair.second.id);
	}
	nextNum =+ 1;
	for (size_t i = 0; i < paths.size(); i++)
	{
		size_t firstExisting = paths[i].position.size();
		size_t lastExisting = paths[i].position.size();
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			std::pair<size_t, size_t> key { i, j };
			if (original.mapping.count(key) == 1)
			{
				if (firstExisting == paths[i].position.size()) firstExisting = j;
				lastExisting = j;
			}
		}
		for (size_t j = firstExisting; j < lastExisting; j++)
		{
			std::pair<size_t, size_t> key { i, j };
			if (original.mapping.count(key) == 1) continue;
			result.mapping[key] = NodePos { nextNum, true };
			nextNum += 1;
		}
	}
	std::cerr << nextNum << " transitive closure sets after inserting middles" << std::endl;
	std::cerr << result.mapping.size() << " transitive closure items after inserting middles" << std::endl;
	return result;
}

std::vector<Alignment> removeHighCoverageAlignments(const std::vector<Path>& paths, const std::vector<Alignment>& alns, size_t maxCoverage)
{
	std::vector<std::vector<size_t>> alnsPerRead;
	std::vector<bool> validAln;
	validAln.resize(alns.size(), true);
	alnsPerRead.resize(paths.size());
	for (size_t i = 0; i < alns.size(); i++)
	{
		alnsPerRead[alns[i].leftPath].push_back(i);
		alnsPerRead[alns[i].rightPath].push_back(i);
	}
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::vector<size_t> startCount;
		std::vector<size_t> endCount;
		startCount.resize(paths[i].position.size(), 0);
		endCount.resize(paths[i].position.size(), 0);
		for (auto alnIndex : alnsPerRead[i])
		{
			auto aln = alns[alnIndex];
			if (aln.leftPath == i)
			{
				startCount[aln.leftStart] += 1;
				endCount[aln.leftEnd] += 1;
			}
			else
			{
				startCount[aln.rightStart] += 1;
				endCount[aln.rightEnd] += 1;
			}
		}
		std::vector<size_t> coverage;
		coverage.resize(paths[i].position.size(), 0);
		coverage[0] = startCount[0];
		for (size_t j = 1; j < coverage.size(); j++)
		{
			coverage[j] = coverage[j-1] + startCount[j] - endCount[j-1];
		}
		for (auto alnIndex : alnsPerRead[i])
		{
			auto aln = alns[alnIndex];
			bool valid = false;
			size_t start, end;
			if (aln.leftPath == i)
			{
				start = aln.leftStart;
				end = aln.leftEnd;
			}
			else
			{
				start = aln.rightStart;
				end = aln.rightEnd;
			}
			for (size_t j = start; j <= end; j++)
			{
				if (coverage[j] <= maxCoverage)
				{
					valid = true;
					break;
				}
			}
			if (!valid)
			{
				validAln[alnIndex] = false;
			}
		}
	}
	std::vector<Alignment> result;
	for (size_t i = 0; i < validAln.size(); i++)
	{
		if (validAln[i]) result.push_back(alns[i]);
	}
	std::cerr << result.size() << " after removing high coverage alignments" << std::endl;
	return result;
}

std::vector<Alignment> removeNonDovetails(const std::vector<Path>& paths, const std::vector<Alignment>& alns)
{
	std::vector<Alignment> result;
	for (auto aln : alns)
	{
		if (aln.leftStart == 0) continue;
		if (aln.leftEnd != paths[aln.leftPath].position.size()-1) continue;
		if (aln.rightReverse)
		{
			if (aln.rightStart == 0) continue;
			if (aln.rightEnd != paths[aln.rightPath].position.size()-1) continue;
		}
		else
		{
			if (aln.rightStart != 0) continue;
			if (aln.rightEnd == paths[aln.rightPath].position.size()-1) continue;
		}
		result.push_back(aln);
	}
	std::cerr << result.size() << " alignments after removing non-dovetails" << std::endl;
	return result;
}

ClosureEdges getClosureEdges(const DoublestrandedTransitiveClosureMapping& closures, const std::vector<Path>& paths)
{
	ClosureEdges result;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 1; j < paths[i].position.size(); j++)
		{
			NodePos oldPos = closures.mapping.at(std::make_pair(i, j-1));
			NodePos newPos = closures.mapping.at(std::make_pair(i, j));
			result.coverage[canon(oldPos, newPos)] += 1;
		}
	}
	std::cerr << result.coverage.size() << " edges" << std::endl;
	return result;
}

ClosureEdges removeChimericEdges(const DoublestrandedTransitiveClosureMapping& closures, const ClosureEdges& edges, size_t maxRemovableCoverage, double fraction)
{
	std::unordered_map<NodePos, size_t> maxOutEdgeCoverage;
	for (auto edge : edges.coverage)
	{
		maxOutEdgeCoverage[edge.first.first] = std::max(maxOutEdgeCoverage[edge.first.first], edge.second);
		maxOutEdgeCoverage[edge.first.second.Reverse()] = std::max(maxOutEdgeCoverage[edge.first.second.Reverse()], edge.second);
	}
	ClosureEdges result;
	for (auto edge : edges.coverage)
	{
		if (edge.second <= maxRemovableCoverage)
		{
			if ((double)edge.second < (double)maxOutEdgeCoverage[edge.first.first] * fraction) continue;
			if ((double)edge.second < (double)maxOutEdgeCoverage[edge.first.second.Reverse()] * fraction) continue;
		}
		result.coverage[edge.first] = edge.second;
	}
	std::cerr << result.coverage.size() << " edges after chimeric removal" << std::endl;
	return result;
}

std::pair<DoublestrandedTransitiveClosureMapping, ClosureEdges> bridgeTips(const DoublestrandedTransitiveClosureMapping& closures, const ClosureEdges& edges, const std::vector<Path>& paths, size_t minCoverage)
{
	std::unordered_set<NodePos> isNotTip;
	for (auto pair : edges.coverage)
	{
		isNotTip.insert(pair.first.first);
		isNotTip.insert(pair.first.second.Reverse());
	}
	std::unordered_map<std::pair<NodePos, NodePos>, std::vector<std::tuple<size_t, size_t, size_t>>> pathsSupportingEdge;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::vector<size_t> gapStarts;
		for (size_t j = 1; j < paths[i].position.size(); j++)
		{
			auto currentKey = std::make_pair(i, j);
			auto previousKey = std::make_pair(i, j-1);
			if (closures.mapping.count(previousKey) == 1 && isNotTip.count(closures.mapping.at(previousKey)) == 0)
			{
				gapStarts.push_back(j-1);
			}
			if (closures.mapping.count(currentKey) == 1 && isNotTip.count(closures.mapping.at(currentKey).Reverse()) == 0)
			{
				auto endPos = closures.mapping.at(currentKey);
				for (auto start : gapStarts)
				{
					auto startPos = closures.mapping.at(std::make_pair(i, start));
					pathsSupportingEdge[canon(startPos, endPos)].emplace_back(i, start, j);
				}
			}
		}
	}
	DoublestrandedTransitiveClosureMapping resultClosures = closures;
	ClosureEdges resultEdges = edges;
	for (auto pair : pathsSupportingEdge)
	{
		std::set<size_t> readsSupportingPath;
		for (auto t : pair.second)
		{
			readsSupportingPath.insert(std::get<0>(t));
		}
		if (readsSupportingPath.size() >= minCoverage)
		{
			resultEdges.coverage[pair.first] = readsSupportingPath.size();
		}
	}
	std::cerr << resultEdges.coverage.size() << " edges after bridging tips" << std::endl;
	return std::make_pair(resultClosures, resultEdges);
}

size_t getLongestOverlap(const std::string& left, const std::string& right, size_t maxOverlap)
{
	assert(left.size() >= maxOverlap);
	assert(right.size() >= maxOverlap);
	for (size_t i = maxOverlap; i > 0; i--)
	{
		if (left.substr(left.size()-maxOverlap) == right.substr(0, maxOverlap)) return i;
	}
	return 0;
}

ClosureEdges determineClosureOverlaps(const std::vector<Path>& paths, const DoublestrandedTransitiveClosureMapping& closures, const ClosureEdges& edges, const GfaGraph& graph)
{
	ClosureEdges result;
	std::unordered_map<size_t, NodePos> closureRepresentsNode;
	for (auto pair : closures.mapping)
	{
		assert(pair.first.first < paths.size());
		assert(pair.first.second < paths[pair.first.first].position.size());
		NodePos pos = paths[pair.first.first].position[pair.first.second];
		assert(graph.nodes.count(pos.id) == 1);
		if (!pair.second.end) pos = pos.Reverse();
		assert(closureRepresentsNode.count(pair.second.id) == 0 || closureRepresentsNode.at(pair.second.id) == pos);
		closureRepresentsNode[pair.second.id] = pos;
	}
	for (auto pair : edges.coverage)
	{
		NodePos fromClosure = pair.first.first;
		NodePos toClosure = pair.first.second;
		if (closureRepresentsNode.count(fromClosure.id) == 0) continue;
		if (closureRepresentsNode.count(toClosure.id) == 0) continue;
		result.coverage[pair.first] = pair.second;
		auto key = std::make_pair(fromClosure, toClosure);
		if (graph.varyingOverlaps.count(key) == 1)
		{
			result.overlap[key] = graph.varyingOverlaps.at(key);
			continue;
		}
		assert(closureRepresentsNode.count(fromClosure.id) == 1);
		assert(closureRepresentsNode.count(toClosure.id) == 1);
		NodePos fromNode = closureRepresentsNode[fromClosure.id];
		if (!fromClosure.end) fromNode = fromNode.Reverse();
		NodePos toNode = closureRepresentsNode[toClosure.id];
		if (!toClosure.end) toNode = toNode.Reverse();
		bool hasEdge = false;
		if (graph.edges.count(fromNode) == 1)
		{
			for (auto target : graph.edges.at(fromNode))
			{
				if (target == toNode) hasEdge = true;
			}
		}
		if (hasEdge)
		{
			result.overlap[key] = graph.edgeOverlap;
			continue;
		}
		assert(graph.nodes.count(fromNode.id) == 1);
		assert(graph.nodes.count(toNode.id) == 1);
		std::string before = graph.nodes.at(fromNode.id);
		if (!fromNode.end) before = CommonUtils::ReverseComplement(before);
		std::string after = graph.nodes.at(toNode.id);
		if (!toNode.end) after = CommonUtils::ReverseComplement(after);
		result.overlap[key] = getLongestOverlap(before, after, graph.edgeOverlap);
	}
	return result;
}

void outputRemappedReads(std::string filename, const std::vector<Path>& paths, const DoublestrandedTransitiveClosureMapping& closures, const ClosureEdges& edges)
{
	std::vector<vg::Alignment> alns;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::vector<std::vector<NodePos>> validSubpaths;
		validSubpaths.emplace_back();
		for (size_t j = 0; j < paths[j].position.size(); j++)
		{
			auto key = std::make_pair(i, j);
			if (closures.mapping.count(key) == 0)
			{
				if (validSubpaths.back().size() != 0) validSubpaths.emplace_back();
				continue;
			}
			if (j > 0)
			{
				auto previousKey = std::make_pair(i, j-1);
				bool broken = false;
				if (closures.mapping.count(previousKey) == 0) broken = true;
				else if (edges.coverage.count(canon(closures.mapping.at(previousKey), closures.mapping.at(key))) == 0) broken = true;
				if (broken)
				{
					if (validSubpaths.back().size() != 0) validSubpaths.emplace_back();
				}
			}
			validSubpaths.back().emplace_back(closures.mapping.at(key));
		}
		size_t currentNum = 0;
		for (size_t j = 0; j < validSubpaths.size(); j++)
		{
			if (validSubpaths[j].size() == 0) continue;
			alns.emplace_back();
			alns.back().set_name(paths[i].name + "_" + std::to_string(currentNum));
			currentNum++;
			for (size_t k = 0; k < validSubpaths[j].size(); k++)
			{
				auto mapping = alns.back().mutable_path()->add_mapping();
				mapping->mutable_position()->set_node_id(validSubpaths[j][k].id);
				mapping->mutable_position()->set_is_reverse(!validSubpaths[j][k].end);
			}
		}
	}

	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alns, 0);
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	std::string inputOverlaps { argv[3] };
	std::string outputGraph { argv[4] };
	std::string outputPaths { argv[5] };
	int maxAlnCount = std::stoi(argv[6]);

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	std::vector<Path> paths;
	{
		auto nodeSizes = getNodeSizes(graph);
		std::cerr << "load paths" << std::endl;
		paths = loadAlignmentsAsPaths(inputAlns, 1000, nodeSizes);
		std::cerr << paths.size() << " paths after filtering by length" << std::endl;
	}
	std::cerr << "pick longest alignments" << std::endl;
	auto pickedAlns = pickLongestPerRead(paths, inputOverlaps, maxAlnCount);
	std::cerr << "get transitive closure" << std::endl;
	auto transitiveClosures = getTransitiveClosures(paths, pickedAlns, inputOverlaps);
	std::cerr << "deallocate picked" << std::endl;
	{
		decltype(pickedAlns) tmp;
		std::swap(pickedAlns, tmp);
	}
	std::cerr << "merge double strands" << std::endl;
	auto doubleStrandedClosures = mergeDoublestrandClosures(paths, transitiveClosures);
	std::cerr << "deallocate one-stranded closures" << std::endl;
	{
		decltype(transitiveClosures) tmp;
		std::swap(transitiveClosures, tmp);
	}
	std::cerr << "get closure edges" << std::endl;
	auto closureEdges = getClosureEdges(doubleStrandedClosures, paths);
	std::cerr << "remove wrong coverage closures" << std::endl;
	doubleStrandedClosures = removeOutsideCoverageClosures(doubleStrandedClosures, 3, 10000);
	std::cerr << "bridge tips" << std::endl;
	std::tie(doubleStrandedClosures, closureEdges) = bridgeTips(doubleStrandedClosures, closureEdges, paths, 2);
	// std::cerr << "insert middles" << std::endl;
	// doubleStrandedClosures = insertMiddles(doubleStrandedClosures, paths);
	std::cerr << "remove chimeric edges" << std::endl;
	closureEdges = removeChimericEdges(doubleStrandedClosures, closureEdges, 5, 0.2);
	std::cerr << "determine closure overlaps" << std::endl;
	closureEdges = determineClosureOverlaps(paths, doubleStrandedClosures, closureEdges, graph);
	std::cerr << "graphify" << std::endl;
	auto result = getGraph(doubleStrandedClosures, closureEdges, paths, graph);
	std::cerr << "output graph" << std::endl;
	result.SaveToFile(outputGraph);
	std::cerr << "output translated paths" << std::endl;
	outputRemappedReads(outputPaths, paths, doubleStrandedClosures, closureEdges);
}