#include <mutex>
#include <concurrentqueue.h>
#include <thread>
#include <iostream>
#include "Assemble.h"

size_t getDiagonalLength(int offset, int width, int height)
{
	assert(offset <= width);
	assert(offset >= -height);
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
	std::vector<size_t> matchesInParallelogram;
	size_t numParallelograms = (leftCumulativePrefixLength.back() + rightCumulativePrefixLength.back() + bandwidth-1) / bandwidth;
	matchesInParallelogram.resize(numParallelograms, 0);
	std::vector<std::tuple<int, size_t, size_t, size_t>> diagonalMatches;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightOccurrences.count(leftPath.position[i]) == 0) continue;
		size_t size = nodeSizes.at(leftPath.position[i].id);
		for (auto j : rightOccurrences.at(leftPath.position[i]))
		{
			assert(j < rightPath.position.size());
			int diagonal = (int)leftCumulativePrefixLength[i] - (int)rightCumulativePrefixLength[j];
			assert(diagonal >= -(int)rightCumulativePrefixLength.back());
			size_t parallelogram = (diagonal + (int)rightCumulativePrefixLength.back()) / bandwidth;
			assert(parallelogram < matchesInParallelogram.size());
			matchesInParallelogram[parallelogram] += size;
		}
	}
	size_t maxParallelogram = 0;
	double bestScore = 0;
	for (size_t i = 0; i < matchesInParallelogram.size(); i++)
	{
		if (matchesInParallelogram[i] == 0) continue;
		int diagonal = (int)i * bandwidth - (int)rightCumulativePrefixLength.back() + bandwidth * 0.5;
		if (i == matchesInParallelogram.size()-1)
		{
			diagonal = (int)i * bandwidth - (int)rightCumulativePrefixLength.back();
		}
		size_t alignmentLength = getDiagonalLength(diagonal, leftCumulativePrefixLength.back(), rightCumulativePrefixLength.back());
		double estimatedScoreHere = matchesInParallelogram[i];
		if (matchesInParallelogram[i] < alignmentLength) estimatedScoreHere -= (double)(alignmentLength - matchesInParallelogram[i]) * mismatchPenalty;
		if (estimatedScoreHere > bestScore)
		{
			bestScore = estimatedScoreHere;
			maxParallelogram = i;
		}
	}
	if (matchesInParallelogram[maxParallelogram] == 0)
	{
		Alignment result;
		result.alignmentLength = 0;
		result.alignmentIdentity = 0;
		return result;
	}
	int diagonalMin = ((int)maxParallelogram - 1) * bandwidth - (int)rightCumulativePrefixLength.back();
	int diagonalMax = ((int)maxParallelogram + 2) * bandwidth - (int)rightCumulativePrefixLength.back();
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
			int diagonal = (int)leftCumulativePrefixLength[i] - (int)rightCumulativePrefixLength[j];
			if (diagonal < diagonalMin || diagonal > diagonalMax) continue;
			if (result.alignedPairs.size() == 0 || (i > result.alignedPairs.back().leftIndex && j > result.alignedPairs.back().rightIndex))
			{
				result.alignedPairs.emplace_back();
				result.alignedPairs.back().leftIndex = i;
				result.alignedPairs.back().rightIndex = j;
				matches += size;
			}
		}
	}
	if (result.alignedPairs.size() == 0)
	{
		result.alignmentLength = 0;
		result.alignmentIdentity = 0;
		return result;
	}
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
		std::unordered_set<NodePos> nodes;
		for (auto node : paths[i].position)
		{
			nodes.insert(node);
		}
		for (auto node : nodes)
		{
			crossesNode[node].push_back(i);
		}
	}
	std::cerr << crossesNode.size() << " crossnodes before repeat filtering" << std::endl;
	{
		std::vector<std::pair<NodePos, size_t>> abundances;
		for (auto pair : crossesNode)
		{
			abundances.emplace_back(pair.first, pair.second.size());
		}
		std::sort(abundances.begin(), abundances.end(), [](const std::pair<NodePos, size_t>& left, const std::pair<NodePos, size_t>& right) { return left.second < right.second; });
		for (auto pair : abundances)
		{
			if (pair.second <= abundances[abundances.size()/2].second * 5) continue;
			crossesNode.erase(pair.first);
		}
	}
	std::cerr << crossesNode.size() << " crossnodes after repeat filtering" << std::endl;
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
		size_t nextID = 0;
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
				alns[i].alignmentID = nextID;
				nextID += 1;
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
				std::unordered_set<size_t> possibleFwMatches;
				std::unordered_set<size_t> possibleBwMatches;
				auto reversePath = paths[i].Reverse();
				for (size_t j = 0; j < paths[i].position.size(); j++)
				{
					auto node = paths[i].position[j];
					size_t nodeSize = nodeSizes.at(paths[i].position[j].id);
					for (auto other : crossesNode[node])
					{
						if (other <= i) continue;
						possibleFwMatches.insert(other);
					}
				}
				for (size_t j = 0; j < reversePath.position.size(); j++)
				{
					auto node = reversePath.position[j];
					size_t nodeSize = nodeSizes.at(reversePath.position[j].id);
					for (auto other : crossesNode[node])
					{
						if (other <= i) continue;
						possibleBwMatches.insert(other);
					}
				}
				if (possibleFwMatches.size() > 0)
				{
					auto occurrences = getOccurrences(paths[i]);
					for (auto j : possibleFwMatches)
					{
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
				}
				if (possibleBwMatches.size() > 0)
				{
					auto reverseOccurrences = getOccurrences(reversePath);
					auto reverseCumulativePrefixLengths = getCumulativePrefixLength(reversePath, nodeSizes);
					for (auto j : possibleBwMatches)
					{
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
		mismatchPenalty = minAlnIdentity / (1.0 - minAlnIdentity);
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