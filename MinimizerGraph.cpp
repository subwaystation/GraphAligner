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

MinimizerGraph::MinimizerGraph(size_t k, size_t w, const GfaGraph& graph) :
k(k),
w(w)
{
	assert(k < sizeof(size_t) * 8 / 2);
	assert(graph.edgeOverlap > k+w);
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

void add(std::vector<size_t>& vec, size_t minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return;
	}
	vec.emplace_back(minmer);
}

void add(std::vector<GraphPos>& vec, GraphPos minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return;
	}
	vec.emplace_back(minmer);
}

void add(std::vector<std::pair<size_t, size_t>>& vec, size_t minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return;
	}
	vec.emplace_back(minmer, length);
}

void add(std::vector<std::pair<GraphPos, size_t>>& vec, GraphPos minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return;
	}
	vec.emplace_back(minmer, length);
}

void MinimizerGraph::initTopology(const GfaGraph& graph)
{
	minmerIndex.resize(pow(2, k * 2));
	minmerIndex.reserve(graph.nodes.size());
	topology.reserve(graph.nodes.size());
	minmersPerEdgeGroup.reserve(graph.nodes.size());

	std::unordered_map<std::string, size_t> edgeIds;
	for (auto node : graph.nodes)
	{
		std::string startEdge = node.second.substr(0, graph.edgeOverlap);
		std::string endEdge = node.second.substr(node.second.size() - graph.edgeOverlap);
		size_t startId = 0;
		size_t endId = 0;
		if (edgeIds.count(startEdge) == 1)
		{
			startId = edgeIds[startEdge];
		}
		else
		{
			minmerIndex.emplace_back();
			topology.emplace_back();
			minmersPerEdgeGroup.emplace_back();

			startId = edgeIds.size();
			doForMinimizers(startEdge, [this, startId](size_t minmer, size_t pos)
			{
				add(minmersPerEdgeGroup[startId], GraphPos { startId, pos, false });
				add(minmerIndex[minmer], GraphPos { startId, pos, false });
			});
			edgeIds[startEdge] = startId;
		}
		if (edgeIds.count(endEdge) == 1)
		{
			endId = edgeIds[endEdge];
		}
		else
		{
			minmerIndex.emplace_back();
			topology.emplace_back();
			minmersPerEdgeGroup.emplace_back();

			endId = edgeIds.size();
			doForMinimizers(endEdge, [this, endId](size_t minmer, size_t pos)
			{
				add(minmersPerEdgeGroup[endId], GraphPos { endId, pos, false });
				add(minmerIndex[minmer], GraphPos { endId, pos, false });
			});
			edgeIds[endEdge] = endId;
		}
		add(topology[endId], startId, node.second.size() - graph.edgeOverlap);
		if (topology.size() % 1000 == 0) std::cerr << minmerIndex.size() << " " << topology.size() << std::endl;
	}
	// size_t numEdges = 0;
	// for (auto pair : minmerTopology)
	// {
	// 	numEdges += pair.second.size();
	// }
	// std::cerr << "number of minimizers: " << minmerTopology.size() << std::endl;
	// std::cerr << "number of edges: " << numEdges << std::endl;
	// std::cerr << "graph density: " << ((double)numEdges / (double)minmerTopology.size() / (double)minmerTopology.size()) << std::endl;
}

size_t MinimizerGraph::minmerize(std::string kmer) const
{
	assert(kmer.size() == k);
	size_t result = 0;
	for (size_t i = 0; i < kmer.size(); i++)
	{
		result |= charNum(kmer[i]) << (i * 2);
	}
	return result;
}

size_t MinimizerGraph::nextminmer(size_t minmer, char newChar) const
{
	minmer >>= 2;
	minmer |= charNum(newChar) << (k * 2);
	return minmer;
}

std::vector<std::pair<GraphPos, size_t>> MinimizerGraph::inNeighbors(GraphPos minmer) const
{
	std::vector<std::pair<GraphPos, size_t>> result;
	for (auto mm : minmersPerEdgeGroup[minmer.node])
	{
		if (mm.pos < minmer.pos) result.emplace_back(mm, minmer.pos - mm.pos);
	}
	for (auto predecessor : topology[minmer.node])
	{
		for (auto mm : minmersPerEdgeGroup[predecessor.first])
		{
			if (mm.pos < minmer.pos + predecessor.second) result.emplace_back(mm, minmer.pos + predecessor.second - mm.pos);
		}
	}
	return result;
}

void MinimizerGraph::align(const std::string& originalSeq) const
{
	// std::string seq = getHPC(originalSeq);
	std::string seq = originalSeq;
	assert(seq.size() > k + w);
// #ifdef PRINTLENS
	auto timeStart = std::chrono::system_clock::now();
	std::vector<size_t> lens;
	lens.resize(seq.size(), 0);
// #endif
	std::unordered_map<GraphPos, size_t> lastMinmerSeqPos;
	std::unordered_map<GraphPos, size_t> lastMinmerMatchLen;
	doForMinimizers(seq, [this, &lens, &lastMinmerSeqPos, &lastMinmerMatchLen](size_t minimizer, size_t seqPos)
	{
		for (auto graphPos : minmerIndex[minimizer])
		{
			auto neighbors = inNeighbors(graphPos);
			size_t currentLen = 0;
// #ifdef PRINTLENS
			size_t lastPos = 0;
// #endif
			for (size_t i = 0; i < neighbors.size(); i++)
			{
				size_t lastLen = lastMinmerMatchLen[neighbors[i].first];
				int seqDist = seqPos - lastMinmerSeqPos[neighbors[i].first];
				int topoDistance = neighbors[i].second;
				// assert(topoDistance <= k+w);
				bool valid = true;
				if (seqDist > 32 || topoDistance > 32)
				{
					if (seqDist < topoDistance * 0.6 || seqDist > topoDistance / 0.6) valid = false;
				}
				// if (seqDist > MAX_SEQ_MINMER_DISTANCE) valid = false;
				if (valid && topoDistance + lastLen > currentLen)
				{
					currentLen = topoDistance + lastLen;
// #ifdef PRINTLENS
					lastPos = lastMinmerSeqPos[neighbors[i].first];
// #endif
				}
			}
			lastMinmerSeqPos[graphPos] = seqPos;
			lastMinmerMatchLen[graphPos] = currentLen;
// #ifdef PRINTLENS
			for (size_t i = lastPos; i < seqPos; i++)
			{
				lens[i] = std::max(lens[i], currentLen);
			}
// #endif
		}
	});
// #ifdef PRINTLENS
	auto timeEnd = std::chrono::system_clock::now();
	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::cerr << "time: " << time << std::endl;
	size_t maxlen = 0;
	for (size_t i = 0; i < lens.size(); i++)
	{
		maxlen = std::max(maxlen, lens[i]);
	}
	std::cerr << "maxlen: " << maxlen << std::endl;
	std::cerr << "lens: " << std::endl;
	for (size_t i = 0; i < lens.size(); i++)
	{
		std::cerr << lens[i] << std::endl;
	}
// #endif
}

GraphPos::GraphPos(int node, size_t pos, bool reverse) :
node(node),
pos(pos),
reverse(reverse)
{
}

bool GraphPos::operator==(const GraphPos& other) const
{
	return node == other.node && pos == other.pos && reverse == other.reverse;
}