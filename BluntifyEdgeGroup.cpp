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

typedef uint32_t NodeType;
enum NodeEnd {End, Start};
typedef std::pair<NodeType, NodeEnd> NodePos;

template <typename T>
class PosVector
{
public:
	size_t size() const
	{
		return first.size();
	}
	void resize(size_t size)
	{
		first.resize(size);
		second.resize(size);
	}
	void resize(size_t size, T val)
	{
		first.resize(size, val);
		second.resize(size, val);
	}
	T& operator[](NodePos pos)
	{
		if (pos.second == NodeEnd::End) return second[pos.first];
		return first[pos.first];
	}
	const T& operator[](NodePos pos) const
	{
		if (pos.second == NodeEnd::End) return second[pos.first];
		return first[pos.first];
	}
private:
	std::vector<T> first;
	std::vector<T> second;
};

template <>
class PosVector<bool>
{
public:
	size_t size() const
	{
		return first.size();
	}
	void resize(size_t size)
	{
		first.resize(size);
		second.resize(size);
	}
	void resize(size_t size, bool val)
	{
		first.resize(size, val);
		second.resize(size, val);
	}
	auto operator[](NodePos pos)
	{
		if (pos.second == NodeEnd::End) return second[pos.first];
		return first[pos.first];
	}
	auto operator[](NodePos pos) const
	{
		if (pos.second == NodeEnd::End) return second[pos.first];
		return first[pos.first];
	}
private:
	std::vector<bool> first;
	std::vector<bool> second;
};

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
	template <typename... Args>
	void push_unique(Args... args)
	{
		T val {args...};
		push_unique(val);
	}
	size_t size() const
	{
		return elems;
	}
	typename std::array<T, N>::iterator end()
	{
		return std::array<T, N>::begin() + elems;
	}
	typename std::array<T, N>::const_iterator end() const
	{
		return std::array<T, N>::begin() + elems;
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
	SizeArray<NodePos, 4> edges;
	SizeArray<NodePos, 4> cliques;
};

class GroupingBluntify
{
public:
	GroupingBluntify() : overlap(-1), items() {};
	int overlap;
	PosVector<EdgeGroupingItem> items;
	std::vector<OriginalPosition> originals;
	SizeArray<NodePos, 4>& neighbors(NodeType pos, NodeEnd side)
	{
		return items[{pos, side}].edges;
	}
	SizeArray<NodePos, 4>& neighbors(NodePos pos)
	{
		return items[pos].edges;
	}
	const SizeArray<NodePos, 4>& neighbors(NodeType pos, NodeEnd side) const
	{
		return items[{pos, side}].edges;
	}
	const SizeArray<NodePos, 4>& neighbors(NodePos pos) const
	{
		return items[pos].edges;
	}
};

class BluntedCliques
{
public:
	std::vector<OriginalPosition> originalNodeAndOffset;
	std::vector<std::pair<NodePos, NodePos>> edges;
};

class NodeClique
{
public:
	std::set<NodePos> left;
	std::set<NodePos> right;
};

GroupingBluntify getSmallerGraph(const BluntedCliques& cliques, int oldOverlap)
{
	assert(oldOverlap > 0);
	GroupingBluntify result;
	result.overlap = oldOverlap-1;
	result.items.resize(cliques.originalNodeAndOffset.size());
	result.originals.resize(cliques.originalNodeAndOffset.size());
	for (size_t i = 0; i < cliques.originalNodeAndOffset.size(); i++)
	{
		result.originals[i] = cliques.originalNodeAndOffset[i];
	}
	for (size_t i = 0; i < cliques.edges.size(); i++)
	{
		auto from = cliques.edges[i].first;
		auto to = cliques.edges[i].second;
		assert(from.first < result.items.size());
		assert(to.first < result.items.size());
		result.items[from].edges.push_unique(to);
		result.items[to].edges.push_unique(from);
	}
	return result;
}

void partitionRec(NodePos pos, bool side, PosVector<bool>& partitioned, PosVector<bool>& handled, const GroupingBluntify& graph)
{
	assert(!handled[pos]);
	handled[pos] = true;
	partitioned[pos] = side;
	for (auto neighbor : graph.neighbors(pos))
	{
		if (handled[neighbor] && neighbor != pos)
		{
			assert(partitioned[neighbor] != side);
		}
		else
		{
			partitionRec(neighbor, !side, partitioned, handled, graph);
		}
	}
}

void partition(PosVector<bool>& partitioned, const GroupingBluntify& graph)
{
	PosVector<bool> handled;
	handled.resize(graph.items.size(), false);
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		if (!handled[{i, NodeEnd::Start}])
		{
			partitionRec({i, NodeEnd::Start}, false, partitioned, handled, graph);
		}
		if (!handled[{i, NodeEnd::End}])
		{
			partitionRec({i, NodeEnd::End}, true, partitioned, handled, graph);
		}
	}
}

void getComponentRec(NodePos pos, PosVector<bool>& partitioned, PosVector<bool>& componentHandled, const GroupingBluntify& graph, std::pair<std::vector<NodePos>, std::vector<NodePos>>& result)
{
	assert(!componentHandled[pos]);
	componentHandled[pos] = true;
	if (partitioned[pos]) result.second.push_back(pos); else result.first.push_back(pos);
	for (auto neighbor : graph.neighbors(pos))
	{
		if (!componentHandled[neighbor]) getComponentRec(neighbor, partitioned, componentHandled, graph, result);
	}
}

std::pair<std::vector<NodePos>, std::vector<NodePos>> getComponent(NodePos pos, PosVector<bool>& partitioned, PosVector<bool>& componentHandled, const GroupingBluntify& graph)
{
	std::pair<std::vector<NodePos>, std::vector<NodePos>> result;
	getComponentRec(pos, partitioned, componentHandled, graph, result);
	return result;
}

bool connected(NodePos first, NodePos second, const GroupingBluntify& graph)
{
	for (auto neighbor : graph.neighbors(first))
	{
		if (neighbor == second) return true;
	}
	return false;
}

bool edgeEqual(NodePos first, NodePos second, const GroupingBluntify& graph)
{
	auto firstNeighbors = graph.neighbors(first);
	auto secondNeighbors = graph.neighbors(second);
	if (firstNeighbors.size() != secondNeighbors.size()) return false;
	for (auto firstNeighbor : firstNeighbors)
	{
		bool found = false;
		for (auto secondNeighbor : secondNeighbors)
		{
			if (secondNeighbor == firstNeighbor) found = true;
		}
		if (!found) return false;
	}
	for (auto secondNeighbor : secondNeighbors)
	{
		bool found = false;
		for (auto firstNeighbor : firstNeighbors)
		{
			if (firstNeighbor == secondNeighbor) found = true;
		}
		if (!found) return false;
	}
	return true;
}

std::vector<NodeClique> solveCliques(const std::pair<std::vector<NodePos>, std::vector<NodePos>>& component, const GroupingBluntify& graph)
{
	if (component.first.size() == 0)
	{
		assert(component.second.size() == 1);
		std::vector<NodeClique> result;
		result.emplace_back();
		result.back().right.insert(component.second[0]);
		return result;
	}
	if (component.second.size() == 0)
	{
		assert(component.first.size() == 1);
		std::vector<NodeClique> result;
		result.emplace_back();
		result.back().left.insert(component.first[0]);
		return result;
	}
	std::vector<size_t> equalToLeft;
	std::vector<size_t> equalToRight;
	equalToLeft.resize(component.first.size());
	equalToRight.resize(component.second.size());
	size_t leftDistinct = 0;
	size_t rightDistinct = 0;
	for (size_t i = 0; i < equalToLeft.size(); i++)
	{
		bool equal = false;
		for (size_t j = 0; j < i; j++)
		{
			if (edgeEqual(component.first[i], component.first[j], graph))
			{
				equal = true;
				equalToLeft[i] = j;
				break;
			}
		}
		if (!equal)
		{
			equalToLeft[i] = i;
			leftDistinct++;
		}
	}
	for (size_t i = 0; i < equalToRight.size(); i++)
	{
		bool equal = false;
		for (size_t j = 0; j < i; j++)
		{
			if (edgeEqual(component.second[i], component.second[j], graph))
			{
				equal = true;
				equalToRight[i] = j;
				break;
			}
		}
		if (!equal)
		{
			equalToRight[i] = i;
			rightDistinct++;
		}
	}
	std::vector<NodeClique> result;
	if (leftDistinct < rightDistinct)
	{
		std::vector<size_t> leftInClique;
		leftInClique.resize(component.first.size());
		for (size_t i = 0; i < component.first.size(); i++)
		{
			if (equalToLeft[i] == i)
			{
				leftInClique[i] = result.size();
				result.emplace_back();
				result.back().left.insert(component.first[i]);
				for (auto other : component.second)
				{
					if (connected(component.first[i], other, graph)) result.back().right.insert(other);
				}
			}
			else
			{
				assert(equalToLeft[i] < i);
				result[equalToLeft[i]].left.insert(component.first[i]);
				leftInClique[i] = leftInClique[equalToLeft[i]];
			}
		}
		assert(result.size() == leftDistinct);
	}
	else
	{
		std::vector<size_t> rightInClique;
		rightInClique.resize(component.second.size());
		for (size_t i = 0; i < component.second.size(); i++)
		{
			if (equalToRight[i] == i)
			{
				rightInClique[i] = result.size();
				result.emplace_back();
				result.back().right.insert(component.second[i]);
				for (auto other : component.first)
				{
					if (connected(component.second[i], other, graph)) result.back().left.insert(other);
				}
			}
			else
			{
				assert(equalToRight[i] < i);
				result[equalToRight[i]].right.insert(component.second[i]);
				rightInClique[i] = rightInClique[equalToRight[i]];
			}
		}
		assert(result.size() == rightDistinct);
	}
	return result;
}

std::vector<NodeClique> getCliques(const GroupingBluntify& graph)
{
	PosVector<bool> partitioned;
	partitioned.resize(graph.items.size(), false);
	partition(partitioned, graph);
	PosVector<bool> componentHandled;
	componentHandled.resize(graph.items.size(), false);
	std::vector<NodeClique> result;
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		if (!componentHandled[{i, NodeEnd::Start}])
		{
			auto component = getComponent({i, NodeEnd::Start}, partitioned, componentHandled, graph);
			auto cliques = solveCliques(component, graph);
			result.insert(result.end(), cliques.begin(), cliques.end());
		}
		if (!componentHandled[{i, NodeEnd::End}])
		{
			auto component = getComponent({i, NodeEnd::End}, partitioned, componentHandled, graph);
			auto cliques = solveCliques(component, graph);
			result.insert(result.end(), cliques.begin(), cliques.end());
		}
	}
	return result;
}

void fillEdgeCliques(GroupingBluntify& graph)
{
	auto cliques = getCliques(graph);
	for (size_t i = 0; i < cliques.size(); i++)
	{
		for (auto left : cliques[i].left)
		{
			graph.items[left].cliques.push_unique(i, NodeEnd::End);
		}
		for (auto right : cliques[i].right)
		{
			graph.items[right].cliques.push_unique(i, NodeEnd::Start);
		}
	}
}

BluntedCliques getBluntedCliques(const GroupingBluntify& grouping)
{
	BluntedCliques result;
	NodeType size = 0;
	for (size_t i = 0; i < grouping.items.size(); i++)
	{
		for (size_t j = 0; j < grouping.items[{i, NodeEnd::Start}].cliques.size(); j++)
		{
			size = std::max(size, grouping.items[{i, NodeEnd::Start}].cliques[j].first + 1);
		}
		for (size_t j = 0; j < grouping.items[{i, NodeEnd::End}].cliques.size(); j++)
		{
			size = std::max(size, grouping.items[{i, NodeEnd::End}].cliques[j].first + 1);
		}
	}
	result.originalNodeAndOffset.resize(size, {0, -1, false});
	for (size_t i = 0; i < grouping.items.size(); i++)
	{
		for (size_t j = 0; j < grouping.items[{i, NodeEnd::End}].cliques.size(); j++)
		{
			for (size_t k = 0; k < grouping.items[{i, NodeEnd::Start}].cliques.size(); k++)
			{ //order doesn't matter (symmetric)
				result.edges.emplace_back(grouping.items[{i, NodeEnd::Start}].cliques[k], grouping.items[{i, NodeEnd::End}].cliques[j]);
			}
		}
		for (size_t j = 0; j < grouping.items[{i, NodeEnd::Start}].cliques.size(); j++)
		{ //todo check
			auto clique = grouping.items[{i, NodeEnd::Start}].cliques[j];
			assert(clique.first < size);
			result.originalNodeAndOffset[clique.first] = grouping.originals[i];
			if (result.originalNodeAndOffset[clique.first].reverse) result.originalNodeAndOffset[clique.first].offset += 1;
			if (clique.second == NodeEnd::Start) result.originalNodeAndOffset[clique.first].reverse = !grouping.originals[i].reverse;
		}
		for (size_t j = 0; j < grouping.items[{i, NodeEnd::End}].cliques.size(); j++)
		{ //todo check
			auto clique = grouping.items[{i, NodeEnd::End}].cliques[j];
			assert(clique.first < size);
			result.originalNodeAndOffset[clique.first] = grouping.originals[i];
			if (!result.originalNodeAndOffset[clique.first].reverse) result.originalNodeAndOffset[clique.first].offset += 1;
			if (clique.second == NodeEnd::End) result.originalNodeAndOffset[clique.first].reverse = !grouping.originals[i].reverse;
			// result.originalNodeAndOffset[clique.first].offset += grouping.originals[i].reverse ? -1 : 1;
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
	std::unordered_map<int, NodeType> originalNodeStart;
	std::unordered_map<int, NodeType> originalNodeEnd;
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
				assert(toend == "+" || toend == "-");
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
				originalNodeStart[id] = resultSize;
				originalNodeEnd[id] = resultSize + nodeSize - 1;
				resultSize += nodeSize;
			}
		}
	}
	assert(resultSize < std::numeric_limits<NodeType>::max());
	result.items.resize(resultSize);
	result.originals.resize(resultSize);
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
				assert(pos < result.items.size());
				result.originals[pos].id = id;
				result.originals[pos].offset = 0;
				result.originals[pos].reverse = false;
				assert(seq.size() > result.overlap);
				auto nodeSize = seq.size() - result.overlap;
				assert(pos + nodeSize <= result.items.size());
				for (size_t i = 1; i < nodeSize; i++)
				{
					result.items[{pos+i-1, NodeEnd::End}].edges.push_unique(pos+i, NodeEnd::Start);
					result.items[{pos+i, NodeEnd::Start}].edges.push_unique(pos+i-1, NodeEnd::End);
					result.originals[pos+i].id = id;
					result.originals[pos+i].offset = i;
					result.originals[pos+i].reverse = false;
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
				NodePos start { originalNodeEnd[fromid], NodeEnd::End };
				NodePos end { originalNodeStart[toid], NodeEnd::Start };
				if (fromstart == "-")
				{
					start = { originalNodeStart[fromid], NodeEnd::Start };
				}
				if (toend == "-")
				{
					end = { originalNodeEnd[toid], NodeEnd::End };
				}
				assert(start.first < result.items.size());
				assert(end.first < result.items.size());
				result.items[start].edges.push_unique(end);
				result.items[end].edges.push_unique(start);
			}
		}
		assert(pos == resultSize);
	}
	assert(result.items.size() == resultSize);
	assert(result.overlap != -1);
	return result;
}

void writeResultGfa(const GroupingBluntify& graph, std::string originalGraph, std::string filename)
{
	assert(graph.overlap == 0);
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
		for (size_t i = 0; i < graph.items.size(); i++)
		{
			auto node = graph.originals[i].id;
			auto pos = graph.originals[i].offset;
			assert(originalSequences.count(node) == 1);
			assert(originalSequences[node].size() > pos);
			char seq = originalSequences[node][pos];
			if (graph.originals[i].reverse)
			{
				switch(seq)
				{
					case 'A':
						seq = 'T';
						break;
					case 'T':
						seq = 'A';
						break;
					case 'C':
						seq = 'G';
						break;
					case 'G':
						seq = 'C';
						break;
					default:
						assert(false);
						std::abort();
				}
			}
			file << "S\t" << (i + 1) << "\t" << seq << std::endl;
			auto neighbors = graph.neighbors({i, NodeEnd::End});
			for (auto neighbor : neighbors)
			{
				file << "L\t" << (i+1) << "\t" << "+" << "\t" << (neighbor.first + 1) << "\t" << (neighbor.second == NodeEnd::End ? "-" : "+") << "\t0M" << std::endl;
			}
			neighbors = graph.neighbors({i, NodeEnd::Start});
			for (auto neighbor : neighbors)
			{
				file << "L\t" << (i+1) << "\t" << "-" << "\t" << (neighbor.first + 1) << "\t" << (neighbor.second == NodeEnd::End ? "-" : "+") << "\t0M" << std::endl;
			}
		}
	}
}

void printGraph(const GroupingBluntify& graph)
{
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		std::cerr << i << ":" << graph.originals[i].id << "," << graph.originals[i].offset << (graph.originals[i].reverse ? "-" : "+");
		std::cerr << ":";
		for (size_t j = 0; j < graph.items[{i, NodeEnd::End}].edges.size(); j++)
		{
			std::cerr << graph.items[{i, NodeEnd::End}].edges[j].first << (graph.items[{i, NodeEnd::End}].edges[j].second == NodeEnd::Start ? "+" : "-") << ",";
		}
		std::cerr << ":";
		for (size_t j = 0; j < graph.items[{i, NodeEnd::Start}].edges.size(); j++)
		{
			std::cerr << graph.items[{i, NodeEnd::Start}].edges[j].first << (graph.items[{i, NodeEnd::Start}].edges[j].second == NodeEnd::Start ? "+" : "-") << ",";
		}
		std::cerr << std::endl;
	}
}

void verifyGraph(const GroupingBluntify& graph)
{
	for (size_t i = 0; i < graph.items.size(); i++)
	{
		for (size_t j = 0; j < graph.items[{i, NodeEnd::Start}].edges.size(); j++)
		{
			auto target = graph.items[{i, NodeEnd::Start}].edges[j];
			bool found = false;
			for (size_t k = 0; k < graph.items[target].edges.size(); k++)
			{
				if (graph.items[target].edges[k].first == i && graph.items[target].edges[k].second == NodeEnd::Start) found = true;
			}
			assert(found);
		}
		for (size_t j = 0; j < graph.items[{i, NodeEnd::End}].edges.size(); j++)
		{
			auto target = graph.items[{i, NodeEnd::End}].edges[j];
			bool found = false;
			for (size_t k = 0; k < graph.items[target].edges.size(); k++)
			{
				if (graph.items[target].edges[k].first == i && graph.items[target].edges[k].second == NodeEnd::End) found = true;
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
		for (size_t j = 0; j < graph.items[{i, NodeEnd::End}].cliques.size(); j++)
		{
			std::cerr << graph.items[{i, NodeEnd::End}].cliques[j].first << (graph.items[{i, NodeEnd::End}].cliques[j].second == NodeEnd::Start ? "+" : "-") << ",";
		}
		std::cerr << ":";
		for (size_t j = 0; j < graph.items[{i, NodeEnd::Start}].cliques.size(); j++)
		{
			std::cerr << graph.items[{i, NodeEnd::Start}].cliques[j].first << (graph.items[{i, NodeEnd::Start}].cliques[j].second == NodeEnd::Start ? "+" : "-") << ",";
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
	writeResultGfa(graph, argv[1], argv[2]);
}