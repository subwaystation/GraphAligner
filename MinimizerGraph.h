#ifndef MinimizerGraph_h
#define MinimizerGraph_h

#include <unordered_map>
#include "GfaGraph.h"

struct GraphPos
{
	GraphPos(int node, size_t pos, bool reverse);
	bool operator==(const GraphPos& other) const;
	int node;
	size_t pos;
	bool reverse;
};

namespace std 
{
	template <> 
	struct hash<GraphPos>
	{
		size_t operator()(const GraphPos& x) const
		{
			return hash<int>()(x.node) ^ hash<size_t>()(x.pos) ^ hash<bool>()(x.reverse);
		}
	};
}

class MinimizerGraph
{
public:
	MinimizerGraph(size_t k, size_t w, const GfaGraph& graph);
	size_t minmerize(std::string kmer) const;
	size_t nextminmer(size_t minmer, char newChar) const;
	std::vector<std::pair<GraphPos, size_t>> inNeighbors(GraphPos minmer) const;
	void align(const std::string& seq) const;
private:
	template <typename F>
	void doForMinimizers(const std::string& seq, F f) const
	{
		std::vector<size_t> window;
		window.resize(w);
		window[0] = minmerize(seq.substr(0, k));
		size_t smallestPos = 0;
		for (size_t i = 1; i < w; i++)
		{
			window[i] = nextminmer(window[i-1], seq[i]);
			if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
		}
		f(window[smallestPos], smallestPos);
		for (size_t seqPos = k+w; seqPos < seq.size(); seqPos++)
		{
			size_t windowpos = (seqPos - k) % w;
			size_t newMinmer = nextminmer(window[(windowpos + w - 1) % w], seq[seqPos]);
			window[windowpos] = newMinmer;
			if (smallestPos == windowpos)
			{
				f(newMinmer, seqPos);
				for (size_t i = 0; i < w; i++)
				{
					if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
				}
			}
			else if (minmerCompare(newMinmer, window[smallestPos]))
			{
				f(newMinmer, seqPos);
			}
		}
	}
	bool minmerCompare(size_t left, size_t right) const;
	void initTopology(const GfaGraph& graph);
	void initRandomOrdering();
	size_t k;
	size_t w;
	std::vector<std::vector<GraphPos>> minmerIndex;
	std::vector<std::vector<std::pair<size_t, size_t>>> topology;
	std::vector<std::vector<GraphPos>> minmersPerEdgeGroup;
	std::vector<size_t> minmerOrdering;
};

#endif
