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

bool operator<(const std::pair<size_t, NodePos>& left, const std::pair<size_t, NodePos>& right)
{
	return left.first < right.first || (left.first == right.first && left.second < right.second);
}

struct Path
{
	std::string name;
	std::vector<NodePos> position;
	std::vector<size_t> nodeSize;
	std::vector<size_t> cumulativePrefixLength;
	std::unordered_map<NodePos, std::vector<size_t>> occurrences;
	Path Reverse() const
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
	void calculateCumulativePrefixLength()
	{
		cumulativePrefixLength.resize(nodeSize.size()+1, 0);
		for (size_t i = 1; i < cumulativePrefixLength.size(); i++)
		{
			cumulativePrefixLength[i] = cumulativePrefixLength[i-1] + nodeSize[i-1];
		}
	}
	void calculateOccurrences()
	{
		occurrences.clear();
		for (size_t i = 0; i < position.size(); i++)
		{
			occurrences[position[i]].push_back(i);
		}
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

bool AlignmentQualityCompareLT(const Alignment& left, const Alignment& right)
{
	return left.alignmentLength * left.alignmentIdentity < right.alignmentLength * right.alignmentIdentity;
}

struct AlignmentComparerLT
{
public:
	AlignmentComparerLT(const std::vector<Alignment>& paths) :
	paths(&paths)
	{}
	AlignmentComparerLT() :
	paths(nullptr)
	{}
	bool operator()(const Alignment& left, const Alignment& right) const
	{
		return AlignmentQualityCompareLT(left, right);
	}
	bool operator()(size_t left, size_t right) const
	{
		assert(paths != nullptr);
		return AlignmentQualityCompareLT(paths->at(left), paths->at(right));
	}
private:
	const std::vector<Alignment>* const paths;
};

bool AlignmentQualityCompareGT(const Alignment& left, const Alignment& right)
{
	return left.alignmentLength * left.alignmentIdentity > right.alignmentLength * right.alignmentIdentity;
}

struct AlignmentComparerGT
{
public:
	AlignmentComparerGT(const std::vector<Alignment>& paths) :
	paths(&paths)
	{}
	AlignmentComparerGT() :
	paths(nullptr)
	{}
	bool operator()(const Alignment& left, const Alignment& right) const
	{
		return AlignmentQualityCompareLT(left, right);
	}
	bool operator()(size_t left, size_t right) const
	{
		assert(paths != nullptr);
		return AlignmentQualityCompareGT(paths->at(left), paths->at(right));
	}
private:
	const std::vector<Alignment>* const paths;
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
	std::cerr << result.size() << " paths" << std::endl;
	return result;
}

Alignment alignSparse(const Path& leftPath, const Path& rightPath, size_t left, size_t right, double mismatchPenalty)
{
	std::vector<std::tuple<size_t, size_t, double, size_t>> trace;
	size_t bestTraceStart = -1;
	double bestStartScore = 0;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightPath.occurrences.count(leftPath.position[i]) == 1)
		{
			for (auto j : rightPath.occurrences.at(leftPath.position[i]))
			{
				trace.emplace_back(i, j, 0, (size_t)-1);
				std::get<2>(trace.back()) = -(double)std::min(leftPath.cumulativePrefixLength[i], rightPath.cumulativePrefixLength[j]) * mismatchPenalty;
				assert(std::get<2>(trace.back()) < leftPath.cumulativePrefixLength.back() && std::get<2>(trace.back()) < rightPath.cumulativePrefixLength.back());
				assert(std::get<2>(trace.back()) > -(double)leftPath.cumulativePrefixLength.back() * mismatchPenalty || std::get<2>(trace.back()) > -(double)rightPath.cumulativePrefixLength.back() * mismatchPenalty);
				for (size_t k = 0; k < trace.size()-1; k++)
				{
					if (std::get<0>(trace[k]) < i && std::get<1>(trace[k]) < j)
					{
						double insertions = leftPath.cumulativePrefixLength[i] - leftPath.cumulativePrefixLength[std::get<0>(trace[k])+1];
						double deletions = rightPath.cumulativePrefixLength[j] - rightPath.cumulativePrefixLength[std::get<1>(trace[k])+1];
						assert(insertions >= 0);
						assert(deletions >= 0);
						assert(insertions <= leftPath.cumulativePrefixLength.back());
						assert(deletions <= rightPath.cumulativePrefixLength.back());
						double btScore = std::get<2>(trace[k]) - std::max(insertions, deletions) * mismatchPenalty;
						if (btScore > std::get<2>(trace.back()))
						{
							std::get<2>(trace.back()) = btScore;
							std::get<3>(trace.back()) = k;
						}
					}
				}
				assert(leftPath.position[i] == rightPath.position[j]);
				std::get<2>(trace.back()) += leftPath.nodeSize[i];
				double startScoreHere = std::get<2>(trace.back()) - (double)std::min(leftPath.cumulativePrefixLength.back() - leftPath.cumulativePrefixLength[i+1], rightPath.cumulativePrefixLength.back() - rightPath.cumulativePrefixLength[j+1]) * mismatchPenalty;
				if (startScoreHere > bestStartScore || bestTraceStart == -1)
				{
					bestTraceStart = trace.size()-1;
					bestStartScore = startScoreHere;
				}
			}
		}
	}
	Alignment result;
	if (trace.size() == 0)
	{
		result.alignmentIdentity = 0;
		return result;
	}
	size_t pos = bestTraceStart;
	size_t i = std::get<0>(trace[pos]);
	size_t j = std::get<1>(trace[pos]);
	size_t matchLen = 0;
	size_t mismatchLen = 0;
	size_t indelSize = std::min(leftPath.cumulativePrefixLength.back() - leftPath.cumulativePrefixLength[i+1], rightPath.cumulativePrefixLength.back() - rightPath.cumulativePrefixLength[j+1]);
	assert(indelSize >= 0);
	assert(indelSize <= leftPath.cumulativePrefixLength.back() || indelSize <= rightPath.cumulativePrefixLength.back());
	result.leftPath = left;
	result.rightPath = right;
	result.leftEnd = i;
	result.rightEnd = j;
	assert(result.leftEnd < leftPath.position.size());
	assert(result.rightEnd < rightPath.position.size());
	// one past the real end
	while (result.leftEnd < leftPath.position.size() && leftPath.cumulativePrefixLength[result.leftEnd+1] - leftPath.cumulativePrefixLength[i+1] <= indelSize) result.leftEnd++;
	while (result.rightEnd < rightPath.position.size() && rightPath.cumulativePrefixLength[result.rightEnd+1] - rightPath.cumulativePrefixLength[j+1] <= indelSize) result.rightEnd++;
	// fix to correct position
	result.leftEnd--;
	result.rightEnd--;
	assert(result.leftEnd == leftPath.position.size()-1 || result.rightEnd == rightPath.position.size()-1);
	mismatchLen = indelSize;
	while (std::get<3>(trace[pos]) != (size_t)-1)
	{
		assert(std::get<3>(trace[pos]) < pos);
		i = std::get<0>(trace[pos]);
		j = std::get<1>(trace[pos]);
		matchLen += leftPath.nodeSize[i];
		result.alignedPairs.emplace_back();
		result.alignedPairs.back().leftIndex = i;
		result.alignedPairs.back().rightIndex = j;
		pos = std::get<3>(trace[pos]);
		size_t nextI = std::get<0>(trace[pos]);
		size_t nextJ = std::get<1>(trace[pos]);
		indelSize = std::max(leftPath.cumulativePrefixLength[i] - leftPath.cumulativePrefixLength[nextI+1], rightPath.cumulativePrefixLength[j] - rightPath.cumulativePrefixLength[nextJ+1]);
		assert(indelSize >= 0);
		assert(indelSize <= leftPath.cumulativePrefixLength.back() || indelSize <= rightPath.cumulativePrefixLength.back());
		mismatchLen += indelSize;
	}
	i = std::get<0>(trace[pos]);
	j = std::get<1>(trace[pos]);
	result.alignedPairs.emplace_back();
	result.alignedPairs.back().leftIndex = i;
	result.alignedPairs.back().rightIndex = j;
	matchLen += leftPath.nodeSize[i];
	indelSize = std::min(leftPath.cumulativePrefixLength[i], rightPath.cumulativePrefixLength[j]);
	assert(indelSize >= 0);
	assert(indelSize <= leftPath.cumulativePrefixLength.back() || indelSize <= rightPath.cumulativePrefixLength.back());
	mismatchLen += indelSize;
	result.leftStart = i;
	result.rightStart = j;
	while (result.leftStart > 0 && leftPath.cumulativePrefixLength[i] - leftPath.cumulativePrefixLength[result.leftStart-1] <= indelSize) result.leftStart--;
	while (result.rightStart > 0 && rightPath.cumulativePrefixLength[j] - rightPath.cumulativePrefixLength[result.rightStart-1] <= indelSize) result.rightStart--;
	assert(result.leftStart == 0 || result.rightStart == 0);
	result.alignmentLength = matchLen + mismatchLen;
	result.alignmentIdentity = (double)matchLen / ((double)mismatchLen + (double)matchLen);
	return result;
}

Alignment align(const std::vector<NodePos>& leftPath, const std::vector<NodePos>& rightPath, const std::vector<size_t>& leftNodeSize, const std::vector<size_t>& rightNodeSize, size_t left, size_t right, double mismatchPenalty)
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
	static thread_local std::vector<std::vector<size_t>> matches;
	if (DPscores.size() < leftPath.size()+1) DPscores.resize(leftPath.size()+1);
	if (DPtrace.size() < leftPath.size()+1) DPtrace.resize(leftPath.size()+1);
	if (matches.size() < leftPath.size()+1) matches.resize(leftPath.size()+1);
	if (DPscores.back().size() < rightPath.size()+1)
	{
		for (size_t i = 0; i < DPscores.size(); i++)
		{
			DPscores[i].resize(rightPath.size()+1, 0);
			DPtrace[i].resize(rightPath.size()+1, Start);
			matches[i].resize(rightPath.size()+1, 0);
		}
	}
	size_t maxI = 0;
	size_t maxJ = 0;
	for (size_t i = 0; i < leftPath.size(); i++)
	{
		for (size_t j = 0; j < rightPath.size(); j++)
		{
			matches[i+1][j+1] = 0;
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
			// if (DPscores[i][j+1] - insertionCost > DPscores[i+1][j+1])
			// {
				DPscores[i+1][j+1] = DPscores[i][j+1] - insertionCost;
				DPtrace[i+1][j+1] = Insertion;
				matches[i+1][j+1] = matches[i][j+1];
			// }
			if (DPscores[i+1][j] - deletionCost > DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i+1][j] - deletionCost;
				DPtrace[i+1][j+1] = Deletion;
				matches[i+1][j+1] = matches[i+1][j];
			}
			if (match && DPscores[i][j] + matchScore >= DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j] + matchScore;
				DPtrace[i+1][j+1] = Match;
				// maxima.erase(std::make_pair(i, j));
				// maxima.emplace(i+1, j+1);
				matches[i+1][j+1] = matches[i][j] + matchScore;
			}
			if (!match && DPscores[i][j] - mismatchCost >= DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j] - mismatchCost;
				DPtrace[i+1][j+1] = Mismatch;
				matches[i+1][j+1] = matches[i][j];
			}
			if ((i == leftPath.size()-1 || j == rightPath.size() - 1) && DPscores[i+1][j+1] >= DPscores[maxI][maxJ])
			{
				maxI = i+1;
				maxJ = j+1;
			}
		}
	}
	Alignment result;
	if (maxI == 0 && maxJ == 0)
	{
		result.alignmentIdentity = 0;
		return result;
	}
	size_t matchLen = 0;
	size_t mismatchLen = 0;
	result.leftPath = left;
	result.rightPath = right;
	result.alignmentLength = 0;
	result.leftEnd = maxI-1;
	result.rightEnd = maxJ-1;
	while (DPtrace[maxI][maxJ] != Start)
	{
		assert(maxI > 0);
		assert(maxJ > 0);
		size_t leftSize = leftNodeSize[maxI-1];
		size_t rightSize = rightNodeSize[maxJ-1];
		result.leftStart = maxI-1;
		result.rightStart = maxJ-1;
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
				result.alignedPairs.emplace_back();
				result.alignedPairs.back().leftIndex = maxI-1;
				result.alignedPairs.back().rightIndex = maxJ-1;
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
	result.alignmentLength = matchLen + mismatchLen;
	if (result.alignmentLength == 0)
	{
		result.alignmentIdentity = 0;
	}
	else
	{
		result.alignmentIdentity = (double)matchLen / ((double)matchLen + (double)mismatchLen);
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
					Alignment fwAln;
					if (pair.second > paths[i].cumulativePrefixLength.back() || pair.second > paths[j].cumulativePrefixLength.back())
					{
						fwAln = align(paths[i].position, paths[j].position, paths[i].nodeSize, paths[j].nodeSize, i, j, mismatchPenalty);
					}
					else
					{
						fwAln = alignSparse(paths[i], paths[j], i, j, mismatchPenalty);
					}
					if (fwAln.alignmentLength >= minAlnLength && fwAln.alignmentIdentity >= minAlnIdentity)
					{
						for (size_t i = 0; i < fwAln.alignedPairs.size(); i++)
						{
							fwAln.alignedPairs[i].leftReverse = false;
							fwAln.alignedPairs[i].rightReverse = false;
						}
						fwAln.rightReverse = false;
						resultsPerThread[thread].push_back(fwAln);
					}
					Alignment bwAln;
					if (pair.second > paths[i].cumulativePrefixLength.back() || pair.second > paths[j].cumulativePrefixLength.back())
					{
						bwAln = align(paths[i].position, reversePaths[j].position, paths[i].nodeSize, reversePaths[j].nodeSize, i, j, mismatchPenalty);
					}
					else
					{
						bwAln = alignSparse(paths[i], reversePaths[j], i, j, mismatchPenalty);
					}
					if (bwAln.alignmentLength >= minAlnLength && bwAln.alignmentIdentity >= minAlnIdentity)
					{
						bwAln.rightStart = paths[j].position.size() - 1 - bwAln.rightStart;
						bwAln.rightEnd = paths[j].position.size() - 1 - bwAln.rightEnd;
						std::swap(bwAln.rightStart, bwAln.rightEnd);
						for (size_t i = 0; i < bwAln.alignedPairs.size(); i++)
						{
							bwAln.alignedPairs[i].leftReverse = false;
							bwAln.alignedPairs[i].rightReverse = true;
							bwAln.alignedPairs[i].rightIndex = paths[j].position.size() - 1 - bwAln.alignedPairs[i].rightIndex;
						}
						bwAln.rightReverse = true;
						resultsPerThread[thread].push_back(bwAln);
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

TransitiveClosureMapping getTransitiveClosures(const std::vector<Path>& paths, const std::vector<Alignment>& alns)
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
	for (auto aln : alns)
	{
		for (auto pair : aln.alignedPairs)
		{
			std::pair<size_t, NodePos> leftKey { aln.leftPath, NodePos { pair.leftIndex, pair.leftReverse } };
			std::pair<size_t, NodePos> rightKey { aln.rightPath, NodePos { pair.rightIndex, pair.rightReverse } };
			set(parent, leftKey, rightKey);
			std::pair<size_t, NodePos> revLeftKey { aln.leftPath, NodePos { pair.leftIndex, !pair.leftReverse } };
			std::pair<size_t, NodePos> revRightKey { aln.rightPath, NodePos { pair.rightIndex, !pair.rightReverse } };
			set(parent, revLeftKey, revRightKey);
		}
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

std::vector<Alignment> pickLongestPerRead(const std::vector<Path>& paths, const std::vector<Alignment>& alns, size_t maxNum)
{
	std::vector<std::vector<size_t>> leftAlnsPerRead;
	std::vector<std::vector<size_t>> rightAlnsPerRead;
	std::vector<int> picked;
	picked.resize(alns.size(), 0);
	leftAlnsPerRead.resize(paths.size());
	rightAlnsPerRead.resize(paths.size());
	for (size_t i = 0; i < alns.size(); i++)
	{
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
	std::vector<Alignment> result;
	for (size_t i = 0; i < alns.size(); i++)
	{
		assert(picked[i] >= 0);
		assert(picked[i] <= 2);
		if (picked[i] == 2) result.push_back(alns[i]);
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
	std::string outputGraph { argv[3] };
	std::string outputPaths { argv[4] };
	size_t minAlnLength = std::stol(argv[5]);
	double minAlnIdentity = std::stod(argv[6]);
	int maxAlnCount = std::stoi(argv[7]);
	int numThreads = std::stoi(argv[8]);

	double mismatchPenalty = 10000;
	if (minAlnIdentity < 1.0)
	{
		mismatchPenalty = 1.0 / (1.0 - minAlnIdentity);
	}

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
	// std::cerr << "remove non-dovetail alignments" << std::endl;
	// alns = removeNonDovetails(paths, alns);
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
	// std::cerr << "double alignments" << std::endl;
	// alns = doubleAlignments(alns);
	std::cerr << "get transitive closure" << std::endl;
	auto transitiveClosures = getTransitiveClosures(paths, alns);
	std::cerr << "deallocate alignments" << std::endl;
	{
		decltype(alns) tmp;
		std::swap(alns, tmp);
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