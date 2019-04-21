#include <mutex>
#include <concurrentqueue.h>
#include <thread>
#include <iostream>
#include "Assemble.h"

template <typename F>
void iterateMostRecentAncestors(std::unordered_set<size_t>& tips, const std::vector<std::vector<size_t>>& inNeighbors, F f)
{
	std::unordered_set<size_t> visited;
	std::vector<size_t> visitStack;
	for (auto tip : tips)
	{
		if (!f(tip))
		{
			for (auto neighbor : inNeighbors[tip])
			{
				visitStack.push_back(neighbor);
			}
		}
	}
	while (visitStack.size() > 0)
	{
		size_t k = visitStack.back();
		visitStack.pop_back();
		if (visited.count(k) == 1) return;
		visited.insert(k);
		if (!f(k))
		{
			for (auto neighbor : inNeighbors[k])
			{
				visitStack.push_back(neighbor);
			}
		}
	}
}

Alignment align(const Path& leftPath, const Path& rightPath, size_t left, size_t right, double mismatchPenalty, int diagonalMin, int diagonalMax)
{
	enum BacktraceType
	{
		Insertion,
		Deletion,
		Match,
		Mismatch,
		Start
	};
	std::vector<std::vector<double>> DPscores;
	std::vector<std::vector<BacktraceType>> DPtrace;
	std::vector<std::vector<size_t>> matches;
	std::vector<size_t> rowOffset;
	DPscores.resize(rightPath.position.size());
	DPtrace.resize(rightPath.position.size());
	matches.resize(rightPath.position.size());
	size_t rowStart = 0;
	size_t rowEnd = 0;
	size_t maxI = -1;
	size_t maxJ = -1;
	for (size_t i = 0; i < rightPath.position.size(); i++)
	{
		while (rowStart < leftPath.position.size() && (int)leftPath.cumulativePrefixLength[rowStart+1] - (int)rightPath.cumulativePrefixLength[i+1] < diagonalMin)
		{
			rowStart++;
		}
		while (rowEnd < leftPath.position.size() && (int)leftPath.cumulativePrefixLength[rowEnd+1] - (int)rightPath.cumulativePrefixLength[i+1] < diagonalMax)
		{
			rowEnd++;
		}
		assert(rowStart == 0);
		assert(rowEnd == leftPath.position.size());
		assert(rowEnd <= leftPath.position.size());
		DPscores[i].resize(rowEnd - rowStart);
		DPtrace[i].resize(rowEnd - rowStart);
		matches[i].resize(rowEnd - rowStart);
		for (size_t j = 0; j < rowEnd - rowStart; j++)
		{
			matches[i][j] = 0;
			DPscores[i][j] = -((double)leftPath.cumulativePrefixLength.back() + (double)rightPath.cumulativePrefixLength.back()) * mismatchPenalty;
			DPtrace[i][j] = Start;
			assert(rowStart+j < leftPath.nodeSize.size());
			assert(i < rightPath.nodeSize.size());
			bool match = (leftPath.position[rowStart+j] == rightPath.position[i]);
			size_t leftSize = leftPath.nodeSize[rowStart+j];
			size_t rightSize = rightPath.nodeSize[i];
			double insertionCost = rightSize * mismatchPenalty;
			double deletionCost = leftSize * mismatchPenalty;
			double mismatchCost = std::max(insertionCost, deletionCost);
			double matchScore = leftSize;
			assert(!match || leftSize == rightSize);
			if (i == 0 || j+rowStart == 0)
			{
				if (match)
				{
					matches[i][j] = matchScore;
					DPscores[i][j] = matchScore;
					DPtrace[i][j] = Match;
				}
				else
				{
					matches[i][j] = 0;
					DPscores[i][j] = -mismatchCost;
					DPtrace[i][j] = Mismatch;
				}
			}
			assert(i == 0 || rowStart >= rowOffset.back());
			size_t previousJ = 0;
			if (i > 0) previousJ = j + (rowStart - rowOffset.back());
			if (i > 0 && previousJ < DPscores[i-1].size() && DPscores[i-1][previousJ] - insertionCost > DPscores[i][j])
			{
				DPscores[i][j] = DPscores[i-1][previousJ] - insertionCost;
				DPtrace[i][j] = Insertion;
				matches[i][j] = matches[i-1][previousJ];
			}
			if (j > 0 && DPscores[i][j-1] - deletionCost > DPscores[i][j])
			{
				DPscores[i][j] = DPscores[i][j-1] - deletionCost;
				DPtrace[i][j] = Deletion;
				matches[i][j] = matches[i][j-1];
			}
			if (i > 0 && previousJ-1 < DPscores[i-1].size() && match && DPscores[i-1][previousJ-1] + matchScore >= DPscores[i][j])
			{
				DPscores[i][j] = DPscores[i-1][previousJ-1] + matchScore;
				DPtrace[i][j] = Match;
				matches[i][j] = matches[i-1][previousJ-1] + matchScore;
			}
			if (i > 0 && previousJ-1 < DPscores[i-1].size() && !match && DPscores[i-1][previousJ-1] - mismatchCost >= DPscores[i][j])
			{
				DPscores[i][j] = DPscores[i-1][previousJ-1] - mismatchCost;
				DPtrace[i][j] = Mismatch;
				matches[i][j] = matches[i-1][previousJ-1];
			}
			if ((i == rightPath.position.size()-1 || rowStart+j == leftPath.position.size() - 1) && (maxI == -1 || DPscores[i][j] >= DPscores[maxI][maxJ]))
			{
				maxI = i;
				maxJ = j;
			}
			assert(DPtrace[i][j] != Start);
		}
		rowOffset.push_back(rowStart);
	}
	Alignment result;
	if (maxI == -1 && maxJ == -1)
	{
		result.alignmentIdentity = 0;
		result.alignmentLength = 0;
		return result;
	}
	size_t matchLen = 0;
	size_t mismatchLen = 0;
	result.leftPath = left;
	result.rightPath = right;
	result.alignmentLength = 0;
	result.leftEnd = rowOffset[maxI] + maxJ;
	result.rightEnd = maxI;
	assert(result.leftEnd == leftPath.position.size()-1 || result.rightEnd == rightPath.position.size()-1);
	bool end = false;
	while (!end)
	{
		assert(maxI < DPscores.size());
		assert(maxI >= 0);
		assert(maxJ >= 0);
		assert(maxJ < DPscores[maxI].size());
		assert(rowOffset[maxI] + maxJ < leftPath.nodeSize.size());
		assert(maxI < rightPath.nodeSize.size());
		size_t leftSize = leftPath.nodeSize[rowOffset[maxI] + maxJ];
		size_t rightSize = rightPath.nodeSize[maxI];
		result.leftStart = rowOffset[maxI] + maxJ;
		result.rightStart = maxI;
		switch(DPtrace[maxI][maxJ])
		{
			case Insertion:
				mismatchLen += rightSize;
				assert(maxI > 0);
				maxI -= 1;
				maxJ += rowOffset[maxI+1] - rowOffset[maxI];
				break;
			case Deletion:
				mismatchLen += leftSize;
				assert(maxJ > 0);
				maxJ -= 1;
				break;
			case Match:
				assert(leftSize == rightSize);
				result.alignedPairs.emplace_back();
				result.alignedPairs.back().leftIndex = rowOffset[maxI] + maxJ;
				result.alignedPairs.back().rightIndex = maxI;
				assert(leftPath.position[result.alignedPairs.back().leftIndex] == rightPath.position[result.alignedPairs.back().rightIndex]);
				matchLen += leftSize;
				if (maxI == 0 || maxJ + rowOffset[maxI] == 0)
				{
					end = true;
					break;
				}
				maxI -= 1;
				maxJ += rowOffset[maxI+1] - rowOffset[maxI];
				maxJ -= 1;
				break;
			case Mismatch:
				mismatchLen += std::max(leftSize, rightSize);
				if (maxI == 0 || maxJ + rowOffset[maxI] == 0)
				{
					end = true;
					break;
				}
				maxI -= 1;
				maxJ += rowOffset[maxI+1] - rowOffset[maxI];
				maxJ -= 1;
				break;
			case Start:
			default:
				assert(false);
		}
	}
	assert(result.leftStart == 0 || result.rightStart == 0);
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

Alignment alignSparse(const Path& leftPath, const Path& rightPath, size_t left, size_t right, double mismatchPenalty, int diagonalMin, int diagonalMax)
{
	std::vector<std::tuple<size_t, size_t, double, size_t>> trace;
	std::unordered_set<size_t> tips;
	std::vector<std::vector<size_t>> inNeighbors;
	size_t bestTraceStart = -1;
	double bestStartScore = 0;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightPath.occurrences.count(leftPath.position[i]) == 1)
		{
			for (auto j : rightPath.occurrences.at(leftPath.position[i]))
			{
				int diagonal = leftPath.cumulativePrefixLength[i] - rightPath.cumulativePrefixLength[j];
				if (diagonal < diagonalMin || diagonal > diagonalMax) continue;
				size_t index = trace.size();
				tips.insert(index);
				inNeighbors.emplace_back();
				trace.emplace_back(i, j, 0, (size_t)-1);
				std::get<2>(trace.back()) = -(double)std::min(leftPath.cumulativePrefixLength[i], rightPath.cumulativePrefixLength[j]) * mismatchPenalty;
				assert(std::get<2>(trace.back()) < leftPath.cumulativePrefixLength.back() && std::get<2>(trace.back()) < rightPath.cumulativePrefixLength.back());
				assert(std::get<2>(trace.back()) > -(double)leftPath.cumulativePrefixLength.back() * mismatchPenalty || std::get<2>(trace.back()) > -(double)rightPath.cumulativePrefixLength.back() * mismatchPenalty);
				std::vector<size_t> removedTips;
				iterateMostRecentAncestors(tips, inNeighbors, [mismatchPenalty, &leftPath, &rightPath, &trace, &inNeighbors, &removedTips, i, j, index](size_t k)
				{
					if (k == index) return false;
					if (std::get<0>(trace[k]) >= i || std::get<1>(trace[k]) >= j) return false;
					removedTips.push_back(k);
					inNeighbors[index].push_back(k);
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
					return true;
				});
				for (auto k : removedTips)
				{
					tips.erase(k);
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

Alignment alignBandSparse(const Path& leftPath, const Path& rightPath, size_t left, size_t right, double mismatchPenalty, size_t minAlnLength, double minAlnIdentity, size_t bandWidth)
{
	std::vector<std::pair<size_t, size_t>> diagonalMatches;
	for (size_t i = 0; i < leftPath.position.size(); i++)
	{
		if (rightPath.occurrences.count(leftPath.position[i]) == 0) continue;
		for (auto j : rightPath.occurrences.at(leftPath.position[i]))
		{
			int diagonal = leftPath.cumulativePrefixLength[i] - rightPath.cumulativePrefixLength[j];
			diagonalMatches.emplace_back(diagonal, leftPath.nodeSize[i]);
		}
	}
	std::sort(diagonalMatches.begin(), diagonalMatches.end(), [](const std::pair<size_t, size_t>& left, const std::pair<size_t, size_t>& right) { return left.first < right.first; });
	if (diagonalMatches.size() == 0)
	{
		Alignment result;
		result.alignmentIdentity = 0;
		result.alignmentLength = 0;
		return result;
	}
	size_t start = 0;
	size_t end = 0;
	size_t currentMatchSize = diagonalMatches[0].second;
	size_t bestStart = 0;
	size_t bestEnd = 0;
	size_t bestMatchSize = 0;
	while (start < diagonalMatches.size())
	{
		while (end < diagonalMatches.size()-1 && diagonalMatches[end+1].first <= diagonalMatches[start].first + bandWidth)
		{
			end++;
			currentMatchSize += diagonalMatches[end].second;
		}
		size_t diagonalSize = std::min(getDiagonalLength(diagonalMatches[start].first, leftPath.cumulativePrefixLength.back(), rightPath.cumulativePrefixLength.back()), getDiagonalLength(diagonalMatches[end].first, leftPath.cumulativePrefixLength.back(), rightPath.cumulativePrefixLength.back()));
		double potentialMatchScoreHere = currentMatchSize;
		if (diagonalSize > currentMatchSize) potentialMatchScoreHere -= (diagonalSize - currentMatchSize) * mismatchPenalty;
		double potentialIdentityHere = (double)currentMatchSize / (double)diagonalSize;
		if (potentialIdentityHere >= minAlnIdentity && currentMatchSize >= bestMatchSize)
		{
			bestMatchSize = currentMatchSize;
			bestStart = start;
			bestEnd = end;
		}
		while (start < diagonalMatches.size()-1 && diagonalMatches[start+1].first == diagonalMatches[start].first)
		{
			currentMatchSize -= diagonalMatches[start].second;
			start++;
		}
		currentMatchSize -= diagonalMatches[start].second;
		start++;
	}
	if (bestMatchSize < minAlnLength)
	{
		Alignment result;
		result.alignmentIdentity = 0;
		result.alignmentLength = 0;
		return result;
	}
	int diagonalMiddle = (diagonalMatches[bestEnd].first - diagonalMatches[bestStart].first) / 2;
	int diagonalMin = diagonalMiddle - (int)bandWidth / 2;
	int diagonalMax = diagonalMiddle + (int)bandWidth / 2;
	if ((bestEnd - bestStart) * (bestEnd - bestStart) < rightPath.position.size())
	{
		return alignSparse(leftPath, rightPath, left, right, mismatchPenalty, diagonalMin, diagonalMax);
	}
	else
	{
		return align(leftPath, rightPath, left, right, mismatchPenalty, diagonalMin, diagonalMax);
	}
}

void induceOverlaps(const std::vector<Path>& paths, double mismatchPenalty, size_t minAlnLength, double minAlnIdentity, int numThreads, std::string tempAlnFileName, size_t bandWidth)
{
	std::unordered_map<size_t, std::vector<size_t>> crossesNode;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (auto node : paths[i].position)
		{
			crossesNode[node.id].push_back(i);
		}
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
		threads.emplace_back([bandWidth, &writequeue, &paths, &nextRead, &nextReadMutex, thread, &crossesNode, minAlnIdentity, minAlnLength, mismatchPenalty]()
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
				auto reversePath = paths[i].Reverse();
				for (auto pair : possibleMatches)
				{
					size_t j = pair.first;
					if (pair.second < minAlnLength) continue;
					if (i == j) continue;
					Alignment fwAln = alignBandSparse(paths[j], paths[i], j, i, mismatchPenalty, minAlnLength, minAlnIdentity, bandWidth);
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
					Alignment bwAln = alignBandSparse(paths[j], reversePath, j, i, mismatchPenalty, minAlnLength, minAlnIdentity, bandWidth);
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
	size_t bandWidth = std::stol(argv[6]);
	std::string outputOverlaps { argv[7] };

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
	std::cerr << paths.size() << " paths" << std::endl;
	std::cerr << "add node lengths" << std::endl;
	paths = addNodeLengths(paths, graph);
	std::cerr << "filter paths on length" << std::endl;
	paths = filterByLength(paths, 1000);
	std::cerr << paths.size() << " alignments after filtering by length" << std::endl;
	std::cerr << "induce overlaps" << std::endl;
	induceOverlaps(paths, mismatchPenalty, minAlnLength, minAlnIdentity, numThreads, outputOverlaps, bandWidth);
}