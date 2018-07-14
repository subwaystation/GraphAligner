#include <unordered_map>
#include <iostream>
#include <chrono>
#include <cmath>
#include "MinimizerGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

constexpr int MAX_GRAPH_TOPO_DISTANCE = 32;
constexpr int MAX_SEQ_MINMER_DISTANCE = 32;

size_t charNum(char c)
{
	switch(c)
	{
		case 'A':
		case 'a':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'G':
		case 'g':
			return 2;
		case 'T':
		case 't':
			return 3;
	}
	assert(false);
	return 0;
}

std::string getHPC(const std::string& str)
{
	std::string result;
	result += str[0];
	for (size_t i = 1; i < str.size(); i++)
	{
		if (str[i] != str[i-1]) result += str[i];
	}
	return result;
}

MinimizerGraph::MinimizerGraph(size_t k, size_t w, const AlignmentGraph& graph) :
k(k),
w(w)
{
	assert(k < sizeof(size_t) * 8 / 2);
	initRandomOrdering();
	initTopology(graph);
}

void MinimizerGraph::initRandomOrdering()
{
	minmerOrdering.resize(pow(2, k * 2));
	for (size_t i = 0; i < minmerOrdering.size(); i++)
	{
		minmerOrdering[i] = i;
	}
	for (size_t i = 0; i < minmerOrdering.size(); i++)
	{
		std::swap(minmerOrdering[i], minmerOrdering[i+(rand() % (minmerOrdering.size()-i))]);
	}
}

bool MinimizerGraph::minmerCompare(size_t left, size_t right) const
{
	return minmerOrdering[left] < minmerOrdering[right];
}

bool addUnique(std::vector<size_t>& vec, size_t minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return false;
	}
	vec.emplace_back(minmer);
	return true;
}

bool addUnique(std::vector<GraphPos>& vec, GraphPos minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return false;
	}
	vec.emplace_back(minmer);
	return true;
}

bool addUnique(std::vector<std::pair<size_t, size_t>>& vec, size_t minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return false;
	}
	vec.emplace_back(minmer, length);
	return true;
}

bool addUnique(std::vector<std::pair<GraphPos, size_t>>& vec, GraphPos minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return false;
	}
	vec.emplace_back(minmer, length);
	return true;
}

void MinimizerGraph::addMinmerIfUnique(const std::vector<size_t>& posStack, size_t minmer, size_t pos, std::unordered_map<size_t, std::vector<size_t>>& endPositions)
{
	assert(pos >= k-1);
	for (auto index : endPositions[posStack[pos]])
	{
		if (minmers[index] == minmer) return;
	}
	endPositions[posStack[pos]].push_back(minmers.size());
	minmerIndex[minmer].push_back(minmers.size());
	minmers.push_back(minmer);
}

void MinimizerGraph::addMinmersRec(size_t pos, std::string& currentMer, std::vector<size_t>& posStack, std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph)
{
	assert(posStack.size() == currentMer.size());
	posStack.push_back(pos);
	currentMer += graph.NodeSequences(pos);
	if (posStack.size() == w+k)
	{
		assert(posStack.size() == currentMer.size());
		assert(posStack.size() == w+k);
		auto found = findOneMinimizer(currentMer);
		addMinmerIfUnique(posStack, found.first, found.second, endPositions);
	}
	else
	{
		if (pos < graph.NodeEnd(graph.IndexToNode(pos))-1)
		{
			addMinmersRec(pos+1, currentMer, posStack, endPositions, graph);
		}
		else
		{
			for (auto neighbor : graph.outNeighbors[graph.IndexToNode(pos)])
			{
				addMinmersRec(graph.NodeStart(neighbor), currentMer, posStack, endPositions, graph);
			}
		}
	}
	currentMer.pop_back();
	posStack.pop_back();
}

void MinimizerGraph::initTopoRec(size_t minmerIndex, size_t pos, size_t len, const std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph)
{
	assert(len < (k + w) * 2);
	if (endPositions.count(pos) == 1)
	{
		assert(endPositions.at(pos).size() > 0);
		for (auto previousIndex : endPositions.at(pos))
		{
			addUnique(topology[minmerIndex], previousIndex, len);
		}
	}
	else
	{
		if (pos == graph.NodeStart(graph.IndexToNode(pos)))
		{
			for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
			{
				initTopoRec(minmerIndex, graph.NodeEnd(neighbor)-1, len+1, endPositions, graph);
			}
		}
		else
		{
			initTopoRec(minmerIndex, pos-1, len+1, endPositions, graph);
		}
	}
}

void MinimizerGraph::expandTopoRec(size_t start, size_t current, size_t len)
{
	for (auto inNeighbor : topology[current])
	{
		assert(inNeighbor.second > 0);
		if (len + inNeighbor.second <= MAX_GRAPH_TOPO_DISTANCE)
		{
			if (addUnique(topology[start], inNeighbor.first, len + inNeighbor.second))
			{
				expandTopoRec(start, inNeighbor.first, len + inNeighbor.second);
			}
		}
	}
}

void MinimizerGraph::initTopology(const AlignmentGraph& graph)
{
	minmerIndex.resize(pow(2, k * 2));
	std::string kmerStack;
	std::vector<size_t> posStack;
	std::cerr << "get minimizers" << std::endl;
	std::unordered_map<size_t, std::vector<size_t>> endPositions;
	size_t processed = 0;
	for (size_t nodeIndex = 0; nodeIndex < graph.nodeStart.size(); nodeIndex++)
	{
		for (size_t i = graph.NodeStart(nodeIndex); i < graph.NodeEnd(nodeIndex); i++)
		{
			addMinmersRec(i, kmerStack, posStack, endPositions, graph);
			assert(posStack.size() == 0);
			assert(kmerStack.size() == 0);
		}
		if (processed % 1000 == 0) std::cerr << "node " << processed << "/" << graph.nodeStart.size() << " minmers " << minmers.size() << " est " << (graph.nodeStart.size() * minmers.size()) / (processed+1) << std::endl;
		processed++;
	}
	topology.resize(minmers.size());
	std::cerr << "get topology" << std::endl;
	processed = 0;
	for (auto pos : endPositions)
	{
		for (auto minmerIndex : pos.second)
		{
			if (pos.first == graph.NodeStart(graph.IndexToNode(pos.first)))
			{
				for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos.first)])
				{
					initTopoRec(minmerIndex, graph.NodeEnd(neighbor)-1, 1, endPositions, graph);
				}
			}
			else
			{
				initTopoRec(minmerIndex, pos.first-1, 1, endPositions, graph);
			}
		}
		if (processed % 1000 == 0) std::cerr << "init topo " << processed << "/" << endPositions.size() << std::endl;
		processed++;
	}
	processed = 0;
	for (size_t i = 0; i < minmers.size(); i++)
	{
		expandTopoRec(i, i, 0);
		if (processed % 1000 == 0) std::cerr << "expand topo " << processed << "/" << minmers.size() << std::endl;
		processed++;
	}
}

size_t MinimizerGraph::minmerize(std::string kmer) const
{
	assert(kmer.size() == k);
	size_t result = 0;
	for (size_t i = 0; i < kmer.size(); i++)
	{
		result |= charNum(kmer[i]) << (i * 2);
	}
	assert(result < pow(4, k));
	return result;
}

size_t MinimizerGraph::nextminmer(size_t minmer, char newChar) const
{
	minmer >>= 2;
	minmer |= charNum(newChar) << ((k-1) * 2);
	assert(minmer < pow(4, k));
	return minmer;
}

std::vector<std::pair<GraphPos, size_t>> MinimizerGraph::inNeighbors(GraphPos minmer) const
{
	std::vector<std::pair<GraphPos, size_t>> result;
	// for (auto mm : minmersPerEdgeGroup[minmer.node])
	// {
	// 	if (mm.pos < minmer.pos) result.emplace_back(mm, minmer.pos - mm.pos);
	// }
	// for (auto predecessor : topology[minmer.node])
	// {
	// 	for (auto mm : minmersPerEdgeGroup[predecessor.first])
	// 	{
	// 		if (mm.pos < minmer.pos + predecessor.second) result.emplace_back(mm, minmer.pos + predecessor.second - mm.pos);
	// 	}
	// }
	return result;
}

void MinimizerGraph::align(const std::string& originalSeq) const
{}

// void MinimizerGraph::align(const std::string& originalSeq) const
// {
// 	// std::string seq = getHPC(originalSeq);
// 	std::string seq = originalSeq;
// 	assert(seq.size() > k + w);
// // #ifdef PRINTLENS
// 	auto timeStart = std::chrono::system_clock::now();
// 	std::vector<size_t> lens;
// 	lens.resize(seq.size(), 0);
// // #endif
// 	std::unordered_map<GraphPos, size_t> lastMinmerSeqPos;
// 	std::unordered_map<GraphPos, size_t> lastMinmerMatchLen;
// 	doForMinimizers(seq, [this, &lens, &lastMinmerSeqPos, &lastMinmerMatchLen](size_t minimizer, size_t seqPos)
// 	{
// 		for (auto graphPos : minmerIndex[minimizer])
// 		{
// 			auto neighbors = inNeighbors(graphPos);
// 			size_t currentLen = 0;
// // #ifdef PRINTLENS
// 			size_t lastPos = 0;
// // #endif
// 			for (size_t i = 0; i < neighbors.size(); i++)
// 			{
// 				size_t lastLen = lastMinmerMatchLen[neighbors[i].first];
// 				int seqDist = seqPos - lastMinmerSeqPos[neighbors[i].first];
// 				int topoDistance = neighbors[i].second;
// 				// assert(topoDistance <= k+w);
// 				bool valid = true;
// 				if (seqDist > 32 || topoDistance > 32)
// 				{
// 					if (seqDist < topoDistance * 0.6 || seqDist > topoDistance / 0.6) valid = false;
// 				}
// 				// if (seqDist > MAX_SEQ_MINMER_DISTANCE) valid = false;
// 				if (valid && topoDistance + lastLen > currentLen)
// 				{
// 					currentLen = topoDistance + lastLen;
// // #ifdef PRINTLENS
// 					lastPos = lastMinmerSeqPos[neighbors[i].first];
// // #endif
// 				}
// 			}
// 			lastMinmerSeqPos[graphPos] = seqPos;
// 			lastMinmerMatchLen[graphPos] = currentLen;
// // #ifdef PRINTLENS
// 			for (size_t i = lastPos; i < seqPos; i++)
// 			{
// 				lens[i] = std::max(lens[i], currentLen);
// 			}
// // #endif
// 		}
// 	});
// // #ifdef PRINTLENS
// 	auto timeEnd = std::chrono::system_clock::now();
// 	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
// 	std::cerr << "time: " << time << std::endl;
// 	size_t maxlen = 0;
// 	for (size_t i = 0; i < lens.size(); i++)
// 	{
// 		maxlen = std::max(maxlen, lens[i]);
// 	}
// 	std::cerr << "maxlen: " << maxlen << std::endl;
// 	std::cerr << "lens: " << std::endl;
// 	for (size_t i = 0; i < lens.size(); i++)
// 	{
// 		std::cerr << lens[i] << std::endl;
// 	}
// // #endif
// }

std::pair<size_t, size_t> MinimizerGraph::findOneMinimizer(const std::string& seq) const
{
	size_t current = minmerize(seq.substr(0, k));
	size_t smallest = current;
	size_t pos = k-1;
	for (size_t i = k; i < seq.size(); i++)
	{
		current = nextminmer(current, seq[i]);
		if (minmerCompare(current, smallest))
		{
			smallest = current;
			pos = i;
		}
	}
	return std::make_pair(smallest, pos);
}

GraphPos::GraphPos(int node, size_t pos) :
node(node),
pos(pos)
{
}

bool GraphPos::operator==(const GraphPos& other) const
{
	return node == other.node && pos == other.pos;
}

bool GraphPos::operator!=(const GraphPos& other) const
{
	return !(*this == other);
}