#include <mutex>
#include <concurrentqueue.h>
#include <thread>
#include <iostream>
#include "Assemble.h"

size_t getDiagonalLength(int offset, int width, int height)
{
	int corner = width - height;
	if (offset >= 0 && offset >= corner)
	{
		assert(offset - corner <= height);
		return height - (offset - corner);
	}
	if (offset >= 0 && offset <= corner)
	{
		return height;
	}
	if (offset <= 0 && offset >= corner)
	{
		return width;
	}
	if (offset <= 0 && offset <= corner)
	{
		assert(corner - offset <= width);
		return width - (corner - offset);
	}
	assert(false);
}

Alignment approxAlign(const Path& leftPath, const Path& rightPath, const std::unordered_map<int, size_t>& nodeSizes, const std::vector<size_t>& leftCumulativePrefixLength, const std::vector<size_t>& rightCumulativePrefixLength, const std::unordered_map<NodePos, std::vector<size_t>>& rightOccurrences, size_t left, size_t right, double mismatchPenalty, int bandwidth)
{
	std::vector<std::pair<int, size_t>> diagonalMatches;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightOccurrences.count(leftPath.position[i]) == 0) continue;
		size_t size = nodeSizes.at(leftPath.position[i].id);
		for (auto j : rightOccurrences.at(leftPath.position[i]))
		{
			diagonalMatches.emplace_back((int)leftCumulativePrefixLength[i] - (int)rightCumulativePrefixLength[j], size);
		}
	}
	std::sort(diagonalMatches.begin(), diagonalMatches.end(), [](const std::pair<int, size_t>& left, const std::pair<int, size_t>& right) { return left.first < right.first; });
	if (diagonalMatches.size() == 0)
	{
		Alignment result;
		result.alignmentIdentity = 0;
		result.alignmentLength = 0;
		return result;
	}
	size_t start = 0;
	size_t end = 0;
	size_t matchesHere = diagonalMatches[0].second;
	size_t bestStart = -1;
	size_t bestEnd = -1;
	double bestScore = 0;
	while (start < diagonalMatches.size())
	{
		while (end < diagonalMatches.size()-1 && diagonalMatches[end+1].first <= diagonalMatches[start].first + bandwidth)
		{
			end++;
			matchesHere += diagonalMatches[end].second;
		}
		assert(end >= start);
		int middle = (diagonalMatches[start].first + diagonalMatches[end].first) / 2.0;
		size_t alignmentLength = getDiagonalLength(middle, leftCumulativePrefixLength.back(), rightCumulativePrefixLength.back());
		double estimatedScoreHere = matchesHere;
		if (matchesHere < alignmentLength) estimatedScoreHere -= (double)(alignmentLength - matchesHere) * mismatchPenalty;
		if (estimatedScoreHere > bestScore)
		{
			bestStart = start;
			bestEnd = end;
			bestScore = estimatedScoreHere;
		}
		while (start < diagonalMatches.size()-1 && diagonalMatches[start+1].first == diagonalMatches[start].first)
		{
			matchesHere -= diagonalMatches[start].second;
			start++;
		}
		assert(start < diagonalMatches.size());
		matchesHere -= diagonalMatches[start].second;
		start++;
	}
	if (bestStart == -1 || bestEnd == -1)
	{
		Alignment result;
		result.alignmentIdentity = 0;
		result.alignmentLength = 0;
		return result;
	}
	assert(bestEnd >= bestStart);
	Alignment result;
	result.alignmentLength = 0;
	result.leftPath = left;
	result.rightPath = right;
	size_t matches = 0;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightOccurrences.count(leftPath.position[i]) == 0) continue;
		size_t size = nodeSizes.at(leftPath.position[i].id);
		for (auto j : rightOccurrences.at(leftPath.position[i]))
		{
			int diagonalHere = (int)leftCumulativePrefixLength[i] - (int)rightCumulativePrefixLength[j];
			if (diagonalHere < diagonalMatches[bestStart].first) continue;
			if (diagonalHere > diagonalMatches[bestEnd].first) continue;
			if (result.alignedPairs.size() == 0 || (i > result.alignedPairs.back().leftIndex && j > result.alignedPairs.back().rightIndex))
			{
				result.alignedPairs.emplace_back();
				result.alignedPairs.back().leftIndex = i;
				result.alignedPairs.back().rightIndex = j;
				matches += nodeSizes.at(leftPath.position[i].id);
			}
		}
	}
	assert(result.alignedPairs.size() <= bestEnd - bestStart + 1);
	assert(result.alignedPairs.size() > 0);
	result.leftStart = result.alignedPairs[0].leftIndex;
	result.rightStart = result.alignedPairs[0].rightIndex;
	result.leftEnd = result.alignedPairs.back().leftIndex;
	result.rightEnd = result.alignedPairs.back().rightIndex;
	size_t startIndelSize = std::min(leftCumulativePrefixLength[result.leftStart], rightCumulativePrefixLength[result.rightStart]);
	size_t i = result.leftStart;
	size_t j = result.rightStart;
	while (result.leftStart > 0 && leftCumulativePrefixLength[i] - leftCumulativePrefixLength[result.leftStart-1] <= startIndelSize) result.leftStart--;
	while (result.rightStart > 0 && rightCumulativePrefixLength[j] - rightCumulativePrefixLength[result.rightStart-1] <= startIndelSize) result.rightStart--;
	assert(result.leftStart == 0 || result.rightStart == 0);
	i = result.leftEnd;
	j = result.rightEnd;
	size_t endIndelSize = std::min(leftCumulativePrefixLength.back() - leftCumulativePrefixLength[i+1], rightCumulativePrefixLength.back() - rightCumulativePrefixLength[j+1]);
	// one past the real end
	while (result.leftEnd < leftPath.position.size() && leftCumulativePrefixLength[result.leftEnd+1] - leftCumulativePrefixLength[i+1] <= endIndelSize) result.leftEnd++;
	while (result.rightEnd < rightPath.position.size() && rightCumulativePrefixLength[result.rightEnd+1] - rightCumulativePrefixLength[j+1] <= endIndelSize) result.rightEnd++;
	// fix to correct position
	result.leftEnd--;
	result.rightEnd--;
	assert(result.leftEnd == leftPath.position.size()-1 || result.rightEnd == rightPath.position.size()-1);
	result.alignmentLength = std::max(leftCumulativePrefixLength[result.leftEnd+1] - leftCumulativePrefixLength[result.leftStart], rightCumulativePrefixLength[result.rightEnd+1] - rightCumulativePrefixLength[result.rightStart]);
	assert(result.alignmentLength >= matches);
	result.alignmentIdentity = (double)matches / (double)result.alignmentLength;
	return result;
}

std::vector<size_t> getCumulativePrefixLength(const Path& path, const std::unordered_map<int, size_t>& nodeSizes)
{
	std::vector<size_t> result;
	result.resize(path.position.size()+1, 0);
	for (int i = 0; i < path.position.size(); i++)
	{
		result[i+1] = result[i] + nodeSizes.at(path.position[i].id);
	}
	return result;
}

std::unordered_map<NodePos, std::vector<size_t>> getOccurrences(const Path& path)
{
	std::unordered_map<NodePos, std::vector<size_t>> result;
	for (size_t i = 0; i < path.position.size(); i++)
	{
		result[path.position[i]].push_back(i);
	}
	return result;
}

void induceOverlaps(const std::vector<Path>& paths, const std::unordered_map<int, size_t>& nodeSizes, double mismatchPenalty, size_t minAlnLength, double minAlnIdentity, int numThreads, std::string tempAlnFileName, int bandwidth)
{
	std::unordered_map<NodePos, std::vector<size_t>> crossesNode;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (auto node : paths[i].position)
		{
			crossesNode[node].push_back(i);
		}
	}
	std::vector<std::vector<size_t>> cumulativePrefixLengths;
	cumulativePrefixLengths.resize(paths.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		cumulativePrefixLengths[i] = getCumulativePrefixLength(paths[i], nodeSizes);
	}
	moodycamel::ConcurrentQueue<Alignment> writequeue;
	std::atomic<bool> overlapsFinished;
	overlapsFinished = false;
	size_t alnCount = 0;
	std::thread overlapWriter { [&alnCount, tempAlnFileName, &overlapsFinished, &writequeue](){
		std::ofstream outfile { tempAlnFileName, std::ios::out | std::ios::binary };
		while (true)
		{
			Alignment alns[100] {};
			size_t gotOverlaps = writequeue.try_dequeue_bulk(alns, 100);
			if (gotOverlaps == 0)
			{
				if (!writequeue.try_dequeue(alns[0]))
				{
					if (overlapsFinished) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				gotOverlaps = 1;
			}
			for (size_t i = 0; i < gotOverlaps; i++)
			{
				WriteAlignment(outfile, alns[i]);
			}
			alnCount += gotOverlaps;
		}
	}};
	std::vector<std::thread> threads;
	std::mutex nextReadMutex;
	size_t nextRead = 0;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([bandwidth, &writequeue, &nodeSizes, &cumulativePrefixLengths, &paths, &nextRead, &nextReadMutex, thread, &crossesNode, minAlnIdentity, minAlnLength, mismatchPenalty]()
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
				std::unordered_map<size_t, size_t> possibleFwMatches;
				std::unordered_map<size_t, size_t> possibleBwMatches;
				auto reversePath = paths[i].Reverse();
				auto occurrences = getOccurrences(paths[i]);
				auto reverseOccurrences = getOccurrences(reversePath);
				auto reverseCumulativePrefixLengths = getCumulativePrefixLength(reversePath, nodeSizes);
				for (size_t j = 0; j < paths[i].position.size(); j++)
				{
					auto node = paths[i].position[j];
					size_t nodeSize = nodeSizes.at(paths[i].position[j].id);
					for (auto other : crossesNode[node])
					{
						if (other <= i) continue;
						possibleFwMatches[other] += nodeSize;
					}
				}
				for (size_t j = 0; j < reversePath.position.size(); j++)
				{
					auto node = reversePath.position[j];
					size_t nodeSize = nodeSizes.at(reversePath.position[j].id);
					for (auto other : crossesNode[node])
					{
						if (other <= i) continue;
						possibleBwMatches[other] += nodeSize;
					}
				}
				for (auto pair : possibleFwMatches)
				{
					size_t j = pair.first;
					if (pair.second < minAlnLength) continue;
					if (i == j) continue;
					Alignment fwAln;
					fwAln = approxAlign(paths[j], paths[i], nodeSizes, cumulativePrefixLengths[j], cumulativePrefixLengths[i], occurrences, j, i, mismatchPenalty, bandwidth);
					if (fwAln.alignmentLength >= minAlnLength && fwAln.alignmentIdentity >= minAlnIdentity)
					{
						for (size_t k = 0; k < fwAln.alignedPairs.size(); k++)
						{
							fwAln.alignedPairs[k].leftReverse = false;
							fwAln.alignedPairs[k].rightReverse = false;
						}
						fwAln.rightReverse = false;
						writequeue.enqueue(fwAln);
					}
				}
				for (auto pair : possibleBwMatches)
				{
					size_t j = pair.first;
					if (pair.second < minAlnLength) continue;
					if (i == j) continue;
					Alignment bwAln;
					bwAln = approxAlign(paths[j], reversePath, nodeSizes, cumulativePrefixLengths[j], reverseCumulativePrefixLengths, reverseOccurrences, j, i, mismatchPenalty, bandwidth);
					if (bwAln.alignmentLength >= minAlnLength && bwAln.alignmentIdentity >= minAlnIdentity)
					{
						bwAln.rightStart = paths[i].position.size() - 1 - bwAln.rightStart;
						bwAln.rightEnd = paths[i].position.size() - 1 - bwAln.rightEnd;
						std::swap(bwAln.rightStart, bwAln.rightEnd);
						for (size_t k = 0; k < bwAln.alignedPairs.size(); k++)
						{
							bwAln.alignedPairs[k].leftReverse = false;
							bwAln.alignedPairs[k].rightReverse = true;
							bwAln.alignedPairs[k].rightIndex = paths[i].position.size() - 1 - bwAln.alignedPairs[k].rightIndex;
						}
						bwAln.rightReverse = true;
						writequeue.enqueue(bwAln);
					}
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	overlapsFinished = true;
	overlapWriter.join();
	std::cerr << alnCount << " induced alignments" << std::endl;
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	size_t minAlnLength = std::stol(argv[3]);
	double minAlnIdentity = std::stod(argv[4]);
	int numThreads = std::stoi(argv[5]);
	int bandwidth = std::stoi(argv[6]);
	std::string outputOverlaps { argv[7] };


	double mismatchPenalty = 10000;
	if (minAlnIdentity < 1.0)
	{
		mismatchPenalty = 1.0 / (1.0 - minAlnIdentity);
	}

	std::unordered_map<int, size_t> nodeSizes;
	{
		std::cerr << "load graph" << std::endl;
		auto graph = GfaGraph::LoadFromFile(inputGraph);
		graph.confirmDoublesidedEdges();
		nodeSizes = getNodeSizes(graph);
	}
	std::cerr << "load paths" << std::endl;
	auto paths = loadAlignmentsAsPaths(inputAlns, 1000, nodeSizes);
	std::cerr << paths.size() << " paths after filtering by length" << std::endl;
	std::cerr << "induce overlaps" << std::endl;
	induceOverlaps(paths, nodeSizes, mismatchPenalty, minAlnLength, minAlnIdentity, numThreads, outputOverlaps, bandwidth);
}