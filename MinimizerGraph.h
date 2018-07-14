#ifndef MinimizerGraph_h
#define MinimizerGraph_h

#include <unordered_map>
#include "GfaGraph.h"

class MinimizerGraph
{
public:
	MinimizerGraph(size_t k, size_t w, const GfaGraph& graph);
	size_t minmerize(std::string kmer) const;
	size_t nextminmer(size_t minmer, char newChar) const;
	std::vector<std::pair<size_t, size_t>> inNeighbors(size_t minmer) const;
	void align(const std::string& seq) const;
private:
	bool minmerCompare(size_t left, size_t right) const;
	void initTopology(const GfaGraph& graph);
	size_t k;
	size_t w;
	std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> minmerTopology;
};

#endif
