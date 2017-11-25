#include <set>
#include <iostream>
#include <map>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <array>
#include <cassert>

template <typename T, size_t N>
class SizeArray : public std::array<T, N>
{
public:
	SizeArray() : elems(0) {};
	void push_back(T val)
	{
		assert(elems < N);
		(*this)[elems] = val;
		elems++;
	}
	void push_unique(T val)
	{
		for (size_t i = 0; i < size(); i++)
		{
			if ((*this)[i] == val) return;
		}
		push_back(val);
	}
	size_t size() const
	{
		return elems;
	}
private:
	char elems;
};

class OriginalPosition
{
public:
	int id;
	int offset;
	bool reverse;
	bool operator==(const OriginalPosition& other) const
	{
		return id == other.id && offset == other.offset && reverse == other.reverse;
	}
	bool operator<(const OriginalPosition& other) const
	{
		return id < other.id || (id == other.id && offset < other.offset) || (id == other.id && offset == other.offset && !reverse && other.reverse);
	}
};

class EdgeGroupingItem
{
public:
	SizeArray<size_t, 4> leftCliques;
	SizeArray<size_t, 4> leftEdges;
	SizeArray<size_t, 4> rightCliques;
	SizeArray<size_t, 4> rightEdges;
	OriginalPosition original;
};

class GroupingBluntify
{
public:
	GroupingBluntify() : overlap(-1), items() {};
	int overlap;
	std::vector<EdgeGroupingItem> items;
};

class BluntedCliques
{
public:
	std::vector<OriginalPosition> originalNodeAndOffset;
	std::vector<std::pair<size_t, size_t>> edges;
};

class NodeClique
{
public:
	std::set<size_t> left;
	std::set<size_t> right;
};

GroupingBluntify getSmallerGraph(const BluntedCliques& cliques, int oldOverlap)
{
	assert(oldOverlap > 0);
	GroupingBluntify result;
	result.overlap = oldOverlap-1;
	result.items.resize(cliques.originalNodeAndOffset.size());
	for (size_t i = 0; i < cliques.originalNodeAndOffset.size(); i++)
	{
		result.items[i].original = cliques.originalNodeAndOffset[i];
	}
	for (size_t i = 0; i < cliques.edges.size(); i++)
	{
		auto from = cliques.edges[i].first;
		auto to = cliques.edges[i].second;
		assert(from < result.items.size());
		assert(to < result.items.size());
		result.items[from].leftEdges.push_back(to);
		result.items[to].rightEdges.push_back(from);
	}
	return result;
}

void fillComponent(std::vector<std::pair<size_t, bool>>& component, const GroupingBluntify& graph, size_t pos, bool side)
{
	if (std::any_of(component.begin(), component.end(), [pos, side](auto pair) { return pair.first == pos && pair.second == side; }))
	{
		return;
	}
	component.emplace_back(pos, side);
	if (side)
	{
		for (size_t i = 0; i < graph.items[pos].rightEdges.size(); i++)
		{
			fillComponent(component, graph, graph.items[pos].rightEdges[i], false);
		}
	}
	else
	{
		for (size_t i = 0; i < graph.items[pos].leftEdges.size(); i++)
		{
			fillComponent(component, graph, graph.items[pos].leftEdges[i], true);
		}
	}
}

std::vector<std::pair<size_t, bool>> getComponent(const GroupingBluntify& graph, size_t pos, bool side)
{
	std::vector<std::pair<size_t, bool>> result;
	fillComponent(result, graph, pos, side);
	return result;
}

bool isClique(int leftSize, int rightSize, uint32_t assignments, int connections)
{
	int assigned[16];
	for (int i = 0; i < 16; i++)
	{
		assigned[i] = (assignments >> (i * 2)) & 0x3;
	}
	for (int clique = 0; clique < 4; clique++)
	{
		SizeArray<int, 4> L;
		SizeArray<int, 4> R;
		for (int edge = 0; edge < 16; edge++)
		{
			if (!(connections & (1 << edge))) continue;
			if (assigned[edge] != clique) continue;
			int i = edge / 4;
			int j = edge % 4;
			if (i < leftSize) L.push_unique(i);
			if (j < rightSize) R.push_unique(j);
		}
		for (int i = 0; i < L.size(); i++)
		{
			for (int j = 0; j < R.size(); j++)
			{
				if (!(connections & (1 << (4 * L[i] + R[j])))) return false;
			}
		}
	}
	return true;
}

//https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSet64
//at most 16 set bits
int popcount(int v)
{
	assert(v <= 65536);
	uint64_t c =  ((v & 0xfff) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f;
	c += (((v & 0xfff000) >> 12) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f;
	return c;
}

uint32_t getLexicographicallyFirstWithNumber(int numCliques)
{
	assert(numCliques == 2 || numCliques == 3);
	if (numCliques == 2) return 1;
	return 4 + 2;
}

uint32_t nextSolution(int leftSize, int rightSize, uint32_t previous, int numCliques)
{
	assert(numCliques == 2 || numCliques == 3);
	previous++;
	for (int i = 0; i < 16; i++)
	{
		if ((previous >> (i * 2)) % 4 == numCliques)
		{
			previous -= numCliques << (i * 2);
			previous += 1 << (i * 2 + 2);
		}
	}
	if (previous >> 30 == 1) return 0;
	return previous;
}

uint32_t getTrivialSolution(int leftSize, int rightSize, int connections)
{
	uint32_t result = 0;
	if (rightSize < leftSize)
	{
		for (int i = 0; i < 16; i++)
		{
			auto val = i % 4;
			if (val >= rightSize) val = 0;
			result += val << (i * 2);
		}
	}
	else
	{
		for (int i = 0; i < 16; i++)
		{
			auto val = i / 4;
			if (val >= leftSize) val = 0;
			result += val << (i * 2);
		}
	}
	return result;
}

uint32_t solveOneClique(int leftSize, int rightSize, int connections)
{
	int maxSize = std::min(leftSize, rightSize);
	if (isClique(leftSize, rightSize, 0, connections)) return 0;
	assert(popcount(connections) >= maxSize);
	for (int numCliques = 2; numCliques < maxSize; numCliques++)
	{
		uint32_t firstSolution = 0;
		if (isClique(leftSize, rightSize, firstSolution, connections)) return firstSolution;
		uint32_t solution = nextSolution(leftSize, rightSize, firstSolution, numCliques);
		while (solution != 0)
		{
			if (isClique(leftSize, rightSize, solution, connections)) return solution;
			solution = nextSolution(leftSize, rightSize, solution, numCliques);
		}
	}
	uint32_t solution = getTrivialSolution(leftSize, rightSize, connections);
	assert(isClique(leftSize, rightSize, solution, connections));
	return solution;
}

std::vector<uint32_t> presolvedCliques()
{
	std::vector<uint32_t> result;
	result.reserve(4 * 4 * 65536);
	for (int leftSize = 1; leftSize <= 4; leftSize++)
	{
		for (int rightSize = 1; rightSize <= 4; rightSize++)
		{
			for (int connections = 0; connections < 65536; connections++)
			{
				result.push_back(solveOneClique(leftSize, rightSize, connections));
			}
		}
	}
	return result;
}

uint32_t getCliqueSolution(size_t problemNumber)
{
	static std::unordered_map<size_t, uint32_t> cache;
	auto found = cache.find(problemNumber);
	if (found == cache.end())
	{
		cache[problemNumber] = solveOneClique(problemNumber / (4 * 65536) + 1, (problemNumber / 65536) % 4 + 1, problemNumber % 65536);
		found = cache.find(problemNumber);
	}
	return found->second;
}

std::vector<NodeClique> solveCliques(const std::vector<std::pair<size_t, bool>>& component, const GroupingBluntify& graph)
{
	SizeArray<size_t, 4> left;
	SizeArray<size_t, 4> right;
	std::array<bool, 16> connected {0};
	for (auto pair : component)
	{
		if (pair.second)
		{
			right.push_back(pair.first);
		}
		else
		{
			left.push_back(pair.first);
		}
	}
	if (left.size() == 0)
	{
		assert(right.size() == 1);
		std::vector<NodeClique> result;
		result.emplace_back();
		result.back().right.insert(right[0]);
		return result;
	}
	if (right.size() == 0)
	{
		assert(left.size() == 1);
		std::vector<NodeClique> result;
		result.emplace_back();
		result.back().left.insert(left[0]);
		return result;
	}
	assert(left.size() >= 1);
	assert(right.size() >= 1);
	for (size_t i = 0; i < left.size(); i++)
	{
		for (size_t j = 0; j < right.size(); j++)
		{
			bool found = false;
			for (size_t k = 0; k < graph.items[left[i]].leftEdges.size(); k++)
			{
				if (graph.items[left[i]].leftEdges[k] == right[j])
				{
					found = true;
				}
			}
			if (found)
			{
				connected[i*4+j] = true;
			}
		}
	}
	size_t problemNumber = 0;
	problemNumber = (left.size()-1) * 65536 * 4 + (right.size()-1) * 65536;
	for (size_t i = 0; i < 16; i++)
	{
		if (connected[i]) problemNumber += (1 << i);
	}
	uint32_t solution = getCliqueSolution(problemNumber);
	std::vector<NodeClique> result;
	result.resize(4);
	for (int i = 0; i < left.size(); i++)
	{
		for (int j = 0; j < right.size(); j++)
		{
			if (!connected[i * 4 + j]) continue;
			auto cliqueNum = (solution >> ((i * 4 + j) * 2)) % 4;
			result[cliqueNum].left.insert(left[i]);
			result[cliqueNum].right.insert(right[j]);
		}
	}
	while (result.size() > 0 && result.back().left.size() == 0 && result.back().right.size() == 0) result.pop_back();
	assert(result.size() > 0);

#ifndef NDEBUG

	for (size_t i = 0; i < result.size(); i++)
	{
		for (auto leftnode : result[i].left)
		{
			for (auto rightnode : result[i].right)
			{
				bool found = false;
				for (size_t k = 0; k < graph.items[leftnode].leftEdges.size(); k++)
				{
					if (graph.items[leftnode].leftEdges[k] == rightnode) found = true;
				}
				assert(found);
			}
		}
	}

#endif

	return result;
}

void fillEdgeCliques(GroupingBluntify& graph)
{
	std::vector<bool> leftHandled;
	std::vector<bool> rightHandled;
	leftHandled.resize(graph.items.size(), false);
	rightHandled.resize(graph.items.size(), false);
	size_t totalCliques = 0;
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		if (!leftHandled[i])
		{
			auto component = getComponent(graph, i, false);
			auto cliques = solveCliques(component, graph);
			for (auto clique : cliques)
			{
				for (auto node : clique.left)
				{
					assert(!leftHandled[node]);
					assert(node < graph.items.size());
					graph.items[node].leftCliques.push_back(totalCliques);
				}
				for (auto node : clique.right)
				{
					assert(!rightHandled[node]);
					assert(node < graph.items.size());
					graph.items[node].rightCliques.push_back(totalCliques);
				}
				totalCliques++;
			}
			for (auto clique : cliques)
			{
				for (auto node : clique.left)
				{
					leftHandled[node] = true;
				}
				for (auto node : clique.right)
				{
					rightHandled[node] = true;
				}
			}
		}
		if (!rightHandled[i])
		{
			auto component = getComponent(graph, i, true);
			auto cliques = solveCliques(component, graph);
			for (auto clique : cliques)
			{
				for (auto node : clique.left)
				{
					assert(!leftHandled[node]);
					assert(node < graph.items.size());
					graph.items[node].leftCliques.push_back(totalCliques);
				}
				for (auto node : clique.right)
				{
					assert(!rightHandled[node]);
					assert(node < graph.items.size());
					graph.items[node].rightCliques.push_back(totalCliques);
				}
				totalCliques++;
			}
			for (auto clique : cliques)
			{
				for (auto node : clique.left)
				{
					leftHandled[node] = true;
				}
				for (auto node : clique.right)
				{
					rightHandled[node] = true;
				}
			}
		}
	}
}

BluntedCliques getBluntedCliques(const GroupingBluntify& grouping)
{
	BluntedCliques result;
	size_t size = 0;
	for (size_t i = 0; i < grouping.items.size(); i++)
	{
		for (size_t j = 0; j < grouping.items[i].leftCliques.size(); j++)
		{
			size = std::max(size, grouping.items[i].leftCliques[j] + 1);
		}
		for (size_t j = 0; j < grouping.items[i].rightCliques.size(); j++)
		{
			size = std::max(size, grouping.items[i].rightCliques[j] + 1);
		}
	}
	result.originalNodeAndOffset.resize(size, {0, -1, false});
	for (size_t i = 0; i < grouping.items.size(); i++)
	{
		for (size_t j = 0; j < grouping.items[i].rightCliques.size(); j++)
		{
			for (size_t k = 0; k < grouping.items[i].leftCliques.size(); k++)
			{
				result.edges.emplace_back(grouping.items[i].rightCliques[j], grouping.items[i].leftCliques[k]);
			}
		}
		for (size_t j = 0; j < grouping.items[i].rightCliques.size(); j++)
		{
			auto clique = grouping.items[i].rightCliques[j];
			assert(clique < size);
			result.originalNodeAndOffset[clique] = grouping.items[i].original;
		}
		for (size_t j = 0; j < grouping.items[i].leftCliques.size(); j++)
		{
			auto clique = grouping.items[i].leftCliques[j];
			assert(clique < size);
			result.originalNodeAndOffset[clique] = grouping.items[i].original;
			result.originalNodeAndOffset[clique].offset += grouping.items[i].original.reverse ? -1 : 1;
		}
	}
#ifndef NDEBUG
	for (size_t i = 0; i < result.originalNodeAndOffset.size(); i++)
	{
		assert(result.originalNodeAndOffset[i].offset != -1);
	}
#endif
	return result;
}

GroupingBluntify loadUnbluntGFA(std::string filename)
{
	size_t resultSize = 0;
	size_t numNodes = 0;
	std::unordered_map<int, size_t> originalNodeStarts;
	std::unordered_map<int, size_t> originalNodeReverse;
	std::unordered_map<int, size_t> originalNodeEnds;
	GroupingBluntify result;
	{
		std::ifstream file { filename };
		while (file.good())
		{
			std::string linestr;
			std::getline(file, linestr);
			std::stringstream line { linestr };
			if (linestr[0] == 'L')
			{
				std::string dummy, fromstart, toend, overlap;
				int fromid, toid;
				line >> dummy >> fromid >> fromstart >> toid >> toend >> overlap;
				assert(fromstart == "+" || fromstart == "-");
				assert(toend == "+" || fromstart == "-");
				int overlapint = std::stoi(overlap.substr(0, overlap.size()-1));
				if (result.overlap == -1) result.overlap = overlapint;
				assert(result.overlap == overlapint);
				break;
			}
		}
	}
	{
		std::ifstream file { filename };
		while (file.good())
		{
			std::string linestr;
			std::getline(file, linestr);
			std::stringstream line { linestr };
			if (linestr[0] == 'S')
			{
				std::string dummy, seq;
				int id;
				line >> dummy >> id >> seq;
				assert(seq.size() > result.overlap);
				auto nodeSize = seq.size() - result.overlap;
				originalNodeStarts[id] = resultSize;
				originalNodeReverse[id] = resultSize + nodeSize;
				resultSize += nodeSize * 2;
				numNodes += 1;
				originalNodeEnds[id] = resultSize - 1;
			}
		}
	}
	result.items.resize(resultSize);
	{
		size_t pos = 0;
		std::ifstream file { filename };
		while (file.good())
		{
			std::string linestr;
			std::getline(file, linestr);
			std::stringstream line { linestr };
			if (linestr[0] == 'S')
			{
				std::string dummy, seq;
				int id;
				line >> dummy >> id >> seq;
				result.items[pos].original.id = id;
				result.items[pos].original.offset = 0;
				result.items[pos].original.reverse = false;
				assert(seq.size() > result.overlap);
				auto nodeSize = seq.size() - result.overlap;
				for (size_t i = 1; i < nodeSize; i++)
				{
					result.items[pos+i-1].leftEdges.push_back(pos+i);
					result.items[pos+i].rightEdges.push_back(pos+i-1);
					result.items[pos+i].original.id = id;
					result.items[pos+i].original.offset = i;
					result.items[pos+i].original.reverse = false;
				}
				pos += nodeSize;
				result.items[pos].original.id = id;
				result.items[pos].original.offset = seq.size() - 1;
				result.items[pos].original.reverse = true;
				for (size_t i = 1; i < nodeSize; i++)
				{
					result.items[pos+i-1].leftEdges.push_back(pos+i);
					result.items[pos+i].rightEdges.push_back(pos+i-1);
					result.items[pos+i].original.id = id;
					result.items[pos+i].original.offset = seq.size() - 1 - i;
					result.items[pos+i].original.reverse = true;
				}
				pos += nodeSize;
			}
			if (linestr[0] == 'L')
			{
				std::string dummy, fromstart, toend, overlap;
				int fromid, toid;
				line >> dummy >> fromid >> fromstart >> toid >> toend >> overlap;
				assert(fromstart == "+" || fromstart == "-");
				assert(toend == "+" || toend == "-");
				auto start1 = originalNodeReverse[fromid] - 1;
				auto end1 = originalNodeStarts[toid];
				auto start2 = originalNodeEnds[toid];
				auto end2 = originalNodeReverse[fromid];
				if (fromstart == "-")
				{
					start1 = originalNodeEnds[fromid];
					end2 = originalNodeStarts[fromid];
				}
				if (toend == "-")
				{
					end1 = originalNodeReverse[toid];
					start2 = originalNodeReverse[toid] - 1;
				}
				result.items[start1].leftEdges.push_unique(end1);
				result.items[start2].leftEdges.push_unique(end2);
				result.items[end1].rightEdges.push_unique(start1);
				result.items[end2].rightEdges.push_unique(start2);
			}
		}
		assert(pos == resultSize);
	}
	assert(result.items.size() == resultSize);
	assert(result.overlap != -1);
	return result;
}

class DoublestrandedResult
{
public:
	std::vector<std::pair<int, int>> nodes;
	std::vector<std::tuple<size_t, bool, size_t, bool>> edges;
};

DoublestrandedResult dontMergeReverses(const GroupingBluntify& graph)
{
	DoublestrandedResult result;
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		result.nodes.emplace_back(graph.items[i].original.id, graph.items[i].original.offset);
		for (size_t j = 0; j < graph.items[i].leftEdges.size(); j++)
		{
			auto target = graph.items[i].leftEdges[j];
			result.edges.emplace_back(i, true, target, true);
		}
	}
	return result;
}

DoublestrandedResult mergeReverses(const GroupingBluntify& graph)
{
	std::map<OriginalPosition, std::set<OriginalPosition>> edges;
	std::map<std::pair<int, int>, size_t> resultPositions;
	DoublestrandedResult result;
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		for (size_t j = 0; j < graph.items[i].leftEdges.size(); j++)
		{
			auto from = graph.items[i].original;
			auto to = graph.items[graph.items[i].leftEdges[j]].original;
			if (from < to)
			{
				std::swap(from, to);
				from.reverse = !from.reverse;
				to.reverse = !to.reverse;
			}
			edges[from].insert(to);
		}
		auto resultPosTuple = std::make_pair(graph.items[i].original.id, graph.items[i].original.offset);
		if (resultPositions.count(resultPosTuple) == 0)
		{
			resultPositions[resultPosTuple] = result.nodes.size();
			result.nodes.emplace_back(resultPosTuple);
		}
	}
	assert(result.nodes.size() == resultPositions.size());
	result.nodes.shrink_to_fit();
	for (auto pair : edges)
	{
		auto fwpostuple = pair.first;
		auto resultpos = resultPositions[std::make_pair(fwpostuple.id, fwpostuple.offset)];
		for (auto target : pair.second)
		{
			auto targetpostuple = std::make_pair(target.id, target.offset);
			result.edges.emplace_back(resultpos, !fwpostuple.reverse, resultPositions[targetpostuple], !target.reverse);
		}
	}
	result.edges.shrink_to_fit();
	return result;
}

void writeResultGfa(const DoublestrandedResult& graph, std::string originalGraph, std::string filename)
{
	std::unordered_map<int, std::string> originalSequences;
	{
		std::ifstream file { originalGraph };
		while (file.good())
		{
			std::string linestr;
			std::getline(file, linestr);
			std::stringstream line { linestr };
			if (linestr[0] == 'S')
			{
				std::string dummy, seq;
				int id;
				line >> dummy >> id >> seq;
				originalSequences[id] = seq;
			}
		}
	}
	{
		std::ofstream file { filename };
		for (size_t i = 0; i < graph.nodes.size(); i++)
		{
			auto node = graph.nodes[i].first;
			auto pos = graph.nodes[i].second;
			assert(originalSequences.count(node) == 1);
			assert(originalSequences[node].size() > pos);
			char seq = originalSequences[node][pos];
			file << "S\t" << i << "\t" << seq << std::endl;
		}
		for (auto edge : graph.edges)
		{
			file << "L\t" << std::get<0>(edge) << "\t" << (std::get<1>(edge) ? "+" : "-") << "\t" << std::get<2>(edge) << "\t" << (std::get<3>(edge) ? "+" : "-") << "\t0M" << std::endl;
		}
	}
}

void printGraph(const GroupingBluntify& graph)
{
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		std::cerr << i << ":" << graph.items[i].original.id << "," << graph.items[i].original.offset << (graph.items[i].original.reverse ? "-" : "+");
		std::cerr << ":";
		for (size_t j = 0; j < graph.items[i].leftEdges.size(); j++)
		{
			std::cerr << graph.items[i].leftEdges[j] << ",";
		}
		std::cerr << std::endl;
	}
}

void verifyGraph(const GroupingBluntify& graph)
{
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		for (size_t j = 0; j < graph.items[i].leftEdges.size(); j++)
		{
			auto target = graph.items[i].leftEdges[j];
			bool found = false;
			for (size_t k = 0; k < graph.items[target].rightEdges.size(); k++)
			{
				if (graph.items[target].rightEdges[k] == i) found = true;
			}
			assert(found);
		}
		for (size_t j = 0; j < graph.items[i].rightEdges.size(); j++)
		{
			auto target = graph.items[i].rightEdges[j];
			bool found = false;
			for (size_t k = 0; k < graph.items[target].leftEdges.size(); k++)
			{
				if (graph.items[target].leftEdges[k] == i) found = true;
			}
			assert(found);
		}
	}
}

void printCliques(const GroupingBluntify& graph)
{
	std::cerr << "cliques" << std::endl;
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		std::cerr << i << ":";
		for (size_t j = 0; j < graph.items[i].leftCliques.size(); j++)
		{
			std::cerr << graph.items[i].leftCliques[j] << ",";
		}
		std::cerr << ":";
		for (size_t j = 0; j < graph.items[i].rightCliques.size(); j++)
		{
			std::cerr << graph.items[i].rightCliques[j] << ",";
		}
		std::cerr << std::endl;
	}
}

int main(int argc, char** argv)
{
	auto graph = loadUnbluntGFA(argv[1]);
	while (graph.overlap > 0)
	{
		std::cerr << "overlap " << graph.overlap << " size " << graph.items.size() << std::endl;
		// printGraph(graph);
		verifyGraph(graph);
		fillEdgeCliques(graph);
		// printCliques(graph);
		auto cliques = getBluntedCliques(graph);
		graph = getSmallerGraph(cliques, graph.overlap);
	}
	std::cerr << "overlap " << graph.overlap << " size " << graph.items.size() << std::endl;
	// printGraph(graph);
	verifyGraph(graph);
	auto doublestranded = dontMergeReverses(graph);
	writeResultGfa(doublestranded, argv[1], argv[2]);
}