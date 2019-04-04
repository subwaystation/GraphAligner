#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "GfaGraph.h"

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right)
{
	return left.first < right.first || (left.first == right.first && left.second < right.second);
}

struct Path
{
	std::vector<NodePos> position;
	std::vector<size_t> nodeSize;
	Path Reverse() const
	{
		Path result = *this;
		std::reverse(result.position.begin(), result.position.end());
		std::reverse(result.nodeSize.begin(), result.nodeSize.end());
		for (size_t i = 0; i < result.position.size(); i++)
		{
			result.position[i] = result.position[i].Reverse();
		}
		return result;
	}
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
};

struct DoublestrandedTransitiveClosureMapping
{
	std::map<std::pair<size_t, size_t>, NodePos> mapping;
};

std::vector<Path> loadAlignmentsAsPaths(std::string fileName)
{
	std::unordered_map<std::string, std::vector<vg::Alignment>> alnsPerRead;
	{
		auto vgAlns = CommonUtils::LoadVGAlignments(fileName);
		for (auto aln : vgAlns)
		{
			alnsPerRead[aln.name()].push_back(aln);
		}
	}
	std::vector<Path> result;
	for (auto pair : alnsPerRead)
	{
		assert(pair.second.size() > 0);
		std::sort(pair.second.begin(), pair.second.end(), [](const vg::Alignment& left, const vg::Alignment& right) { return left.query_position() < right.query_position(); });
		result.emplace_back();
		for (auto aln : pair.second)
		{
			for (int i = 0; i < aln.path().mapping_size(); i++)
			{
				result.back().position.emplace_back(aln.path().mapping(i).position().node_id(), !aln.path().mapping(i).position().is_reverse());
			}
		}
		if (result.back().position.size() == 0) result.pop_back();
	}
	std::cerr << result.size() << " paths" << std::endl;
	return result;
}

std::vector<Alignment> align(const std::vector<NodePos>& leftPath, const std::vector<NodePos>& rightPath, const std::vector<size_t>& leftNodeSize, const std::vector<size_t>& rightNodeSize, size_t left, size_t right, double mismatchPenalty)
{
	enum BacktraceType
	{
		Insertion,
		Deletion,
		Match,
		Mismatch,
		Start
	};
	static thread_local std::vector<std::vector<double>> DPscores;
	static thread_local std::vector<std::vector<BacktraceType>> DPtrace;
	if (DPscores.size() < leftPath.size()+1) DPscores.resize(leftPath.size()+1);
	if (DPtrace.size() < leftPath.size()+1) DPtrace.resize(leftPath.size()+1);
	if (DPscores.back().size() < rightPath.size()+1)
	{
		for (size_t i = 0; i < DPscores.size(); i++)
		{
			DPscores[i].resize(rightPath.size()+1, 0);
			DPtrace[i].resize(rightPath.size()+1, Start);
		}
	}
	std::set<std::pair<size_t, size_t>> maxima;
	for (size_t i = 0; i < leftPath.size(); i++)
	{
		for (size_t j = 0; j < rightPath.size(); j++)
		{
			DPscores[i+1][j+1] = 0;
			DPtrace[i+1][j+1] = Start;
			bool match = (leftPath[i] == rightPath[j]);
			size_t leftSize = leftNodeSize[i];
			size_t rightSize = rightNodeSize[j];
			double insertionCost = leftSize * mismatchPenalty;
			double deletionCost = rightSize * mismatchPenalty;
			double mismatchCost = std::max(insertionCost, deletionCost);
			double matchScore = leftSize;
			assert(!match || leftSize == rightSize);
			if (DPscores[i][j+1] - insertionCost > DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j+1] - insertionCost;
				DPtrace[i+1][j+1] = Insertion;
			}
			if (DPscores[i+1][j] - deletionCost > DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i+1][j] - deletionCost;
				DPtrace[i+1][j+1] = Deletion;
			}
			if (match && DPscores[i][j] + matchScore >= DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j] + matchScore;
				DPtrace[i+1][j+1] = Match;
				maxima.erase(std::make_pair(i, j));
				maxima.emplace(i+1, j+1);
			}
			if (!match && DPscores[i][j] - mismatchCost >= DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j] - mismatchCost;
				DPtrace[i+1][j+1] = Mismatch;
			}
		}
	}
	size_t matchLen = 0;
	size_t mismatchLen = 0;
	std::vector<Alignment> result;
	for (auto maximum : maxima)
	{
		size_t maxI = maximum.first;
		size_t maxJ = maximum.second;
		result.emplace_back();
		result.back().leftPath = left;
		result.back().rightPath = right;
		result.back().alignmentLength = 0;
		result.back().leftEnd = maxI-1;
		result.back().rightEnd = maxJ-1;
		while (DPtrace[maxI][maxJ] != Start)
		{
			assert(maxI > 0);
			assert(maxJ > 0);
			size_t leftSize = leftNodeSize[maxI-1];
			size_t rightSize = rightNodeSize[maxJ-1];
			result.back().leftStart = maxI-1;
			result.back().rightStart = maxJ-1;
			switch(DPtrace[maxI][maxJ])
			{
				case Insertion:
					mismatchLen += leftSize;
					maxI -= 1;
					continue;
				case Deletion:
					mismatchLen += rightSize;
					maxJ -= 1;
					continue;
				case Match:
					assert(leftSize == rightSize);
					result.back().alignedPairs.emplace_back();
					result.back().alignedPairs.back().leftIndex = maxI-1;
					result.back().alignedPairs.back().rightIndex = maxJ-1;
					matchLen += leftSize;
					maxI -= 1;
					maxJ -= 1;
					continue;
				case Mismatch:
					mismatchLen += std::max(leftSize, rightSize);
					maxI -= 1;
					maxJ -= 1;
					continue;
				case Start:
				default:
					assert(false);
			}
		}
		result.back().alignmentLength = matchLen + mismatchLen;
		if (result.back().alignmentLength == 0)
		{
			result.back().alignmentIdentity = 0;
		}
		else
		{
			result.back().alignmentIdentity = (double)matchLen / ((double)matchLen + (double)mismatchLen);
		}
	}
	return result;
}

std::vector<Alignment> induceOverlaps(const std::vector<Path>& paths, double mismatchPenalty, size_t minAlnLength, double minAlnIdentity, int numThreads)
{
	std::vector<Alignment> result;
	std::vector<Path> reversePaths;
	reversePaths.reserve(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		reversePaths.push_back(paths[i].Reverse());
	}
	std::unordered_map<size_t, std::vector<size_t>> crossesNode;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (auto node : paths[i].position)
		{
			crossesNode[node.id].push_back(i);
		}
	}
	std::vector<std::vector<Alignment>> resultsPerThread;
	resultsPerThread.resize(numThreads);
	std::vector<std::thread> threads;
	std::mutex nextReadMutex;
	size_t nextRead = 0;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&paths, &nextRead, &nextReadMutex, &resultsPerThread, thread, &crossesNode, minAlnIdentity, minAlnLength, mismatchPenalty, &reversePaths]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> guard { nextReadMutex };
					i = nextRead;
					nextRead += 1;
				}
				if (i >= paths.size()) break;
				std::cerr << i << "/" << paths.size() << std::endl;
				std::unordered_map<size_t, size_t> possibleMatches;
				for (size_t j = 0; j < paths[i].position.size(); j++)
				{
					auto node = paths[i].position[j];
					size_t nodeSize = paths[i].nodeSize[j];
					for (auto other : crossesNode[node.id])
					{
						if (other <= i) continue;
						possibleMatches[other] += nodeSize;
					}
				}
				for (auto pair : possibleMatches)
				{
					size_t j = pair.first;
					if (pair.second < minAlnLength) continue;
					if (i == j) continue;
					auto fwAlns = align(paths[i].position, paths[j].position, paths[i].nodeSize, paths[j].nodeSize, i, j, mismatchPenalty);
					for (auto aln : fwAlns)
					{
						if (aln.alignmentLength < minAlnLength) continue;
						if (aln.alignmentIdentity < minAlnIdentity) continue;
						for (size_t i = 0; i < aln.alignedPairs.size(); i++)
						{
							aln.alignedPairs[i].leftReverse = false;
							aln.alignedPairs[i].rightReverse = false;
						}
						aln.rightReverse = false;
						resultsPerThread[thread].push_back(aln);
					}
					auto bwAlns = align(paths[i].position, reversePaths[j].position, paths[i].nodeSize, reversePaths[j].nodeSize, i, j, mismatchPenalty);
					for (auto aln : bwAlns)
					{
						if (aln.alignmentLength < minAlnLength) continue;
						if (aln.alignmentIdentity < minAlnIdentity) continue;
						aln.rightStart = paths[j].position.size() - 1 - aln.rightStart;
						aln.rightEnd = paths[j].position.size() - 1 - aln.rightEnd;
						std::swap(aln.rightStart, aln.rightEnd);
						for (size_t i = 0; i < aln.alignedPairs.size(); i++)
						{
							aln.alignedPairs[i].leftReverse = false;
							aln.alignedPairs[i].rightReverse = true;
							aln.alignedPairs[i].rightIndex = paths[j].position.size() - 1 - aln.alignedPairs[i].rightIndex;
						}
						aln.rightReverse = true;
						resultsPerThread[thread].push_back(aln);
					}
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
		result.insert(result.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
	}
	std::cerr << result.size() << " induced alignments" << std::endl;
	return result;
}

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

TransitiveClosureMapping getTransitiveClosures(const std::vector<Path>& paths, const std::vector<Alignment>& alns)
{
	TransitiveClosureMapping result;
	std::map<std::pair<size_t, NodePos>, std::pair<size_t, NodePos>> parent;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			find(parent, std::pair<size_t, NodePos> { i, NodePos { j, true } });
			find(parent, std::pair<size_t, NodePos> { i, NodePos { j, false } });
		}
	}
	for (auto aln : alns)
	{
		for (auto pair : aln.alignedPairs)
		{
			std::pair<size_t, NodePos> leftKey { aln.leftPath, NodePos { pair.leftIndex, pair.leftReverse } };
			std::pair<size_t, NodePos> rightKey { aln.rightPath, NodePos { pair.rightIndex, pair.rightReverse } };
			set(parent, leftKey, rightKey);
		}
	}
	std::map<std::pair<size_t, NodePos>, size_t> closureNumber;
	size_t nextClosure = 0;
	for (auto key : parent)
	{
		auto found = find(parent, key.first);
		if (closureNumber.count(found) == 0)
		{
			closureNumber[found] = nextClosure;
			nextClosure += 1;
		}
		result.mapping[key.first] = closureNumber.at(found);
	}
	std::cerr << nextClosure << " transitive closure sets" << std::endl;
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
		result.tags[pair.second.id] = "LN:i:" + std::to_string(seq.size() - graph.edgeOverlap) + "\tKC:i:" + std::to_string((seq.size() - graph.edgeOverlap) * closureCoverage[pair.second.id]) + "\tkm:f:" + std::to_string(closureCoverage[pair.second.id]);
		outputtedClosures.insert(pair.second.id);
	}
	std::cerr << outputtedClosures.size() << " outputted closures" << std::endl;
	for (auto pair : edges.coverage)
	{
		result.edges[pair.first.first].push_back(pair.first.second);
	}
	std::cerr << edges.coverage.size() << " outputted edges" << std::endl;
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
	}
	return result;
}

DoublestrandedTransitiveClosureMapping mergeDoublestrandClosures(const std::vector<Path>& paths, const TransitiveClosureMapping& original)
{
	DoublestrandedTransitiveClosureMapping result;
	std::unordered_map<size_t, NodePos> mapping;
	int nextId = 0;
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
	std::cerr << nextId << " doublestranded transitive closure sets" << std::endl;
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

std::vector<Alignment> pickLongestPerRead(const std::vector<Path>& paths, const std::vector<Alignment>& alns, size_t maxNum)
{
	std::vector<std::vector<size_t>> alnsPerRead;
	std::vector<bool> picked;
	picked.resize(alns.size(), false);
	alnsPerRead.resize(paths.size());
	for (size_t i = 0; i < alns.size(); i++)
	{
		alnsPerRead[alns[i].leftPath].push_back(i);
		alnsPerRead[alns[i].rightPath].push_back(i);
	}
	for (size_t i = 0; i < alnsPerRead.size(); i++)
	{
		if (alnsPerRead[i].size() > maxNum)
		{
			std::sort(alnsPerRead[i].begin(), alnsPerRead[i].end(), [&alns](size_t left, size_t right){ return alns[left].alignmentLength * alns[left].alignmentIdentity < alns[right].alignmentLength * alns[right].alignmentIdentity; });
			// std::sort(alnsPerRead[i].begin(), alnsPerRead[i].end(), [](const Alignment& left, const Alignment& right){ return left.alignmentIdentity < right.alignmentIdentity; });
			// std::stable_sort(alnsPerRead[i].begin(), alnsPerRead[i].end(), [](const Alignment& left, const Alignment& right){ return left.alignmentLength < right.alignmentLength; });
			for (size_t j = alnsPerRead[i].size() - maxNum; j < alnsPerRead[i].size(); j++)
			{
				picked[alnsPerRead[i][j]] = true;
			}
		}
		else
		{
			for (auto aln : alnsPerRead[i])
			{
				picked[aln] = true;
			}
		}
	}
	std::vector<Alignment> result;
	for (size_t i = 0; i < alns.size(); i++)
	{
		if (picked[i]) result.push_back(alns[i]);
	}
	std::cerr << result.size() << " alignments after picking longest" << std::endl;
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
	std::cerr << result.size() << " alignments after filtering by length" << std::endl;
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
	for (auto pair : closures.mapping)
	{
		if (coverage[pair.second.id] >= minCoverage && coverage[pair.second.id] <= maxCoverage)
		{
			result.mapping[pair.first] = pair.second;
		}
	}
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

ClosureEdges removeChimericEdges(const DoublestrandedTransitiveClosureMapping& closures, const ClosureEdges& edges, size_t maxRemovableCoverage, size_t minSafeCoverage)
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
			if (maxOutEdgeCoverage[edge.first.first] >= minSafeCoverage) continue;
			if (maxOutEdgeCoverage[edge.first.second.Reverse()] >= minSafeCoverage) continue;
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
			auto previousKey = std::make_pair(i, j);
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

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	size_t alignmentEndCut = std::stol(argv[3]);
	double mismatchPenalty = std::stod(argv[4]);
	size_t minAlnLength = std::stol(argv[5]);
	double minAlnIdentity = std::stod(argv[6]);
	std::string outputGraph { argv[7] };
	int numThreads = std::stoi(argv[8]);
	int maxAlnCount = std::stoi(argv[9]);

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	std::cerr << "load paths" << std::endl;
	auto paths = loadAlignmentsAsPaths(inputAlns);
	std::cerr << "add node lengths" << std::endl;
	paths = addNodeLengths(paths, graph);
	std::cerr << "filter paths on length" << std::endl;
	paths = filterByLength(paths, 1000);
	std::cerr << "induce overlaps" << std::endl;
	auto alns = induceOverlaps(paths, mismatchPenalty, minAlnLength, minAlnIdentity, numThreads);
	std::cerr << "remove non-dovetail alignments" << std::endl;
	alns = removeNonDovetails(paths, alns);
	std::cerr << "pick longest alignments" << std::endl;
	alns = pickLongestPerRead(paths, alns, maxAlnCount);
	// std::cerr << "double alignments" << std::endl;
	// alns = doubleAlignments(alns);
	// std::cerr << "remove contained alignments" << std::endl;
	// alns = removeContained(paths, alns);
	// std::cerr << "double alignments" << std::endl;
	// alns = doubleAlignments(alns);
	// std::cerr << "remove high coverage alignments" << std::endl;
	// alns = removeHighCoverageAlignments(paths, alns, 40);
	// std::cerr << "pick lowest error alignments" << std::endl;
	// alns = pickLowestErrorPerRead(paths, alns, 3);
	std::cerr << "double alignments" << std::endl;
	alns = doubleAlignments(alns);
	std::cerr << "get transitive closure" << std::endl;
	auto transitiveClosures = getTransitiveClosures(paths, alns);
	std::cerr << "merge double strands" << std::endl;
	auto doubleStrandedClosures = mergeDoublestrandClosures(paths, transitiveClosures);
	std::cerr << "get closure edges" << std::endl;
	auto closureEdges = getClosureEdges(doubleStrandedClosures, paths);
	std::cerr << "remove wrong coverage closures" << std::endl;
	doubleStrandedClosures = removeOutsideCoverageClosures(doubleStrandedClosures, 3, 10000);
	// std::cerr << "insert middles" << std::endl;
	// doubleStrandedClosures = insertMiddles(doubleStrandedClosures, paths);
	std::cerr << "remove chimeric edges" << std::endl;
	closureEdges = removeChimericEdges(doubleStrandedClosures, closureEdges, 3, 10);
	// std::cerr << "bridge tips" << std::endl;
	// std::tie(doubleStrandedClosures, closureEdges) = bridgeTips(doubleStrandedClosures, closureEdges, paths, 2);
	std::cerr << "graphify" << std::endl;
	auto result = getGraph(doubleStrandedClosures, closureEdges, paths, graph);
	std::cerr << "output" << std::endl;
	result.SaveToFile(outputGraph);
}