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
	struct PathTree
	{
		std::vector<size_t> nodeLabels;
		std::vector<std::unordered_map<size_t, size_t>> children;
		std::vector<size_t> nodeDepths;
		std::vector<std::vector<std::pair<size_t, size_t>>> leafMinmers;
	};
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
		window.resize(w-k+1, 0);
		assert(w-k+1 < sizeof(size_t) * 8);
		size_t printed = 0;
		window[0] = minmerize(seq.substr(0, k));
		size_t smallestPos = k-1;
		for (size_t i = 1; i < w-k+1; i++)
		{
			window[i] = nextminmer(window[i-1], seq[k-1+i]);
			if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
		}
		assert(minmerOrdering[window[smallestPos]] != std::numeric_limits<size_t>::max());
		auto test = findOneMinimizer(seq.substr(0, w));
		assert(test.first == window[smallestPos]);
		printed |= 1 << smallestPos;
		f(window[smallestPos], k + smallestPos - 1);
		for (size_t seqPos = w; seqPos < seq.size(); seqPos++)
		{
			size_t windowpos = (seqPos - k + 1) % window.size();
			size_t newMinmer = nextminmer(window[(windowpos + window.size() - 1) % window.size()], seq[seqPos]);
			window[windowpos] = newMinmer;
			printed &= ~(1 << windowpos);
			if (smallestPos == windowpos)
			{
				assert(minmerOrdering[newMinmer] != std::numeric_limits<size_t>::max());
				for (size_t i = 0; i < window.size(); i++)
				{
					if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
				}
				std::pair<size_t, size_t> test = findOneMinimizer(seq.substr(seqPos - w + 1, w));
				assert(test.first == window[smallestPos]);
				for (size_t i = windowpos+1; i <= windowpos+window.size(); i++)
				{
					if (window[i % window.size()] == window[smallestPos] && ((printed & (1 << smallestPos)) == 0))
					{
						assert(seqPos + i >= windowpos + window.size());
						printed |= 1 << smallestPos;
						f(window[smallestPos], seqPos + i - (windowpos + window.size()));
					}
				}
			}
			else if (minmerCompare(newMinmer, window[smallestPos]) || newMinmer == window[smallestPos])
			{
				smallestPos = windowpos;
				assert(minmerOrdering[newMinmer] != std::numeric_limits<size_t>::max());
				std::pair<size_t, size_t> test = findOneMinimizer(seq.substr(seqPos - w + 1, w));
				assert(test.first == newMinmer);
				printed |= 1 << smallestPos;
				f(newMinmer, seqPos);
			}
		}
	}
	void expandTopoRec(size_t start, size_t current, size_t len);
	void initTopoRec(size_t minmerIndex, size_t pos, size_t len, const std::unordered_map<size_t, std::vector<size_t>>& endPositions);
	void addMinmerIfUnique(const std::vector<size_t>& posStack, size_t minmer, size_t pos, std::unordered_map<size_t, std::vector<size_t>>& endPositions);
	void addMinmersRec(size_t pos, std::string& currentMer, std::vector<size_t>& posStack, std::unordered_map<size_t, std::vector<size_t>>& endPositions);
	void addTopoRec(size_t startNode, std::vector<GraphPos>& posStack, const std::unordered_map<GraphPos, std::vector<size_t>>& endPositions);
	bool minmerCompare(size_t left, size_t right) const;
	void initTopology();
	void initRandomOrdering();
	void initTree();
	void buildTreeRec(size_t parent, const std::vector<std::pair<size_t, size_t>>& activeIndices, const std::vector<size_t>& uniqueMinmers);
	void getTotalPathsRec(size_t node, std::vector<size_t>& paths) const;
	void verifyMinmerPositions(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const;
	void verifyExpandedTopology(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const;
	void verifyDirectTopology(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const;
	bool checkMinmerExistsRec(const std::string& minmer, size_t pos, size_t minmerPos) const;
	void checkExpandedTopologyRec(const std::unordered_map<size_t, std::vector<size_t>>& endPositions, size_t pos, size_t lengthLeft, std::set<size_t>& foundNeighbors) const;
	void checkDirectTopologyRec(const std::unordered_map<size_t, std::vector<size_t>>& endPositions, size_t pos, size_t lengthLeft, std::set<size_t>& foundNeighbors) const;
	bool findPath(const std::string& seq, std::vector<size_t>& posStack, const AlignmentGraph& graph) const;
	bool tryLoadOrdering(std::string filename);
	bool tryLoadTopology(std::string filename);
	bool tryLoadTree(std::string filename);
	void storeTree(std::string filename);
	void storeTopology(std::string filename);
	void storeOrdering(std::string filename);
	size_t k;
	size_t w;
	std::vector<size_t> minmers;
	std::vector<std::vector<size_t>> minmerIndex;
	std::vector<size_t> minmerPosition;
	std::vector<std::vector<std::pair<size_t, size_t>>> topology;
	std::vector<std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>> reverseTopology;
	std::vector<size_t> minmerOrdering;
	PathTree tree;
	const AlignmentGraph& graph;
};

#endif
