#ifndef MinimizerGraph_h
#define MinimizerGraph_h

#include <unordered_map>
#include "AlignmentGraph.h"

struct GraphPos
{
	GraphPos(int node, size_t pos);
	bool operator==(const GraphPos& other) const;
	bool operator!=(const GraphPos& other) const;
	int node;
	size_t pos;
};

namespace std 
{
	template <> 
	struct hash<GraphPos>
	{
		size_t operator()(const GraphPos& x) const
		{
			return hash<int>()(x.node) ^ hash<size_t>()(x.pos);
		}
	};
}

class MinimizerGraph
{
public:
	MinimizerGraph(size_t k, size_t w, const AlignmentGraph& graph);
	size_t minmerize(std::string kmer) const;
	size_t nextminmer(size_t minmer, char newChar) const;
	std::vector<std::pair<size_t, size_t>> inNeighbors(size_t minmer) const;
	void align(const std::string& seq) const;
private:
	std::pair<size_t, size_t> findOneMinimizer(const std::string& seq) const;
	template <typename F>
	void doForMinimizers(const std::string& seq, F f) const
	{
		std::vector<size_t> window;
		window.resize(w, 0);
		window[0] = minmerize(seq.substr(0, k));
		size_t smallestPos = 0;
		for (size_t i = 1; i < w; i++)
		{
			window[i] = nextminmer(window[i-1], seq[i]);
			if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
		}
		assert(minmerOrdering[window[smallestPos]] != std::numeric_limits<size_t>::max());
		f(window[smallestPos], k + smallestPos - 1);
		for (size_t seqPos = k+w; seqPos < seq.size(); seqPos++)
		{
			size_t windowpos = (seqPos - k) % w;
			size_t newMinmer = nextminmer(window[(windowpos + w - 1) % w], seq[seqPos]);
			window[windowpos] = newMinmer;
			if (smallestPos == windowpos)
			{
				assert(minmerOrdering[newMinmer] != std::numeric_limits<size_t>::max());
				f(newMinmer, seqPos);
				for (size_t i = 0; i < w; i++)
				{
					if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
				}
			}
			else if (minmerCompare(newMinmer, window[smallestPos]))
			{
				assert(minmerOrdering[newMinmer] != std::numeric_limits<size_t>::max());
				f(newMinmer, seqPos);
			}
		}
	}
	void expandTopoRec(size_t start, size_t current, size_t len);
	void initTopoRec(size_t minmerIndex, size_t pos, size_t len, const std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph);
	void addMinmerIfUnique(const std::vector<size_t>& posStack, size_t minmer, size_t pos, std::unordered_map<size_t, std::vector<size_t>>& endPositions);
	void addMinmersRec(size_t pos, std::string& currentMer, std::vector<size_t>& posStack, std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph);
	void addTopoRec(size_t startNode, std::vector<GraphPos>& posStack, const std::unordered_map<GraphPos, std::vector<size_t>>& endPositions, const AlignmentGraph& graph);
	bool minmerCompare(size_t left, size_t right) const;
	void initTopology(const AlignmentGraph& graph);
	void initRandomOrdering();
	size_t k;
	size_t w;
	std::vector<size_t> minmers;
	std::vector<std::vector<size_t>> minmerIndex;
	std::vector<std::vector<std::pair<size_t, size_t>>> topology;
	std::vector<size_t> minmerOrdering;
};

#endif
