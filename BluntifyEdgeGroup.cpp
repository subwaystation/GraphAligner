#include <iostream>
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
	size_t size() const
	{
		return elems;
	}
private:
	char elems;
};

class EdgeGroupingItem
{
public:
	SizeArray<size_t, 4> leftCliques;
	SizeArray<size_t, 4> leftEdges;
	SizeArray<size_t, 4> rightCliques;
	SizeArray<size_t, 4> rightEdges;
	int originalNodeId;
	int originalNodeOffset;
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
	std::vector<std::pair<size_t, int>> originalNodeAndOffset;
	std::vector<std::pair<size_t, size_t>> edges;
};

class NodeClique
{
public:
	std::vector<size_t> left;
	std::vector<size_t> right;
};

GroupingBluntify getSmallerGraph(const BluntedCliques& cliques, int oldOverlap)
{
	assert(oldOverlap > 0);
	GroupingBluntify result;
	result.overlap = oldOverlap-1;
	result.items.resize(cliques.originalNodeAndOffset.size());
	for (size_t i = 0; i < cliques.originalNodeAndOffset.size(); i++)
	{
		result.items[i].originalNodeId = cliques.originalNodeAndOffset[i].first;
		result.items[i].originalNodeOffset = cliques.originalNodeAndOffset[i].second;
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

bool isClique(int leftSize, int rightSize, uint16_t assignments, int connections)
{
	int L[4], R[4];
	for (int i = 0; i < 4; i++)
	{
		L[i] = (assignments >> (i * 2)) & 0x3;
		R[i] = (assignments >> (8 + i * 2)) & 0x3;
	}
	for (int i = 0; i < leftSize; i++)
	{
		for (int j = 0; j < rightSize; j++)
		{
			if (L[i] != R[j]) continue;
			if (!(connections & (1 << (i*4 + j)))) return false;
		}
	}
	return true;
}

uint16_t getLexicographicallyFirstWithNumber(int numCliques)
{
	return numCliques-1;
}

uint16_t nextSolution(int leftSize, int rightSize, uint16_t previous, int numCliques)
{
	int R[] {0, 0, 0, 0};
	int L[] {0, 0, 0, 0};
	bool has[] {false, false, false, false};
	for (int i = 0; i < 4; i++)
	{
		L[i] = (previous >> (i * 2)) & 0x3;
		R[i] = (previous >> (8 + i * 2)) & 0x3;
	}
	L[0]++;
	for (int i = 0; i < leftSize && i < 3; i++)
	{
		if (L[i] == numCliques)
		{
			L[i] = 0;
			L[i+1]++;
			has[L[i]] = true;
		}
	}
	if (L[3] == numCliques)
	{
		L[3] = 0;
		R[0]++;
	}
	has[L[3]] = true;
	for (int i = 0; i < rightSize && i < 3; i++)
	{
		if (R[i] == numCliques)
		{
			R[i] = 0;
			R[i+1]++;
			has[L[i]] = true;
		}
	}
	if (R[3] == numCliques)
	{
		return 0;
	}
	has[R[3]] = true;
	uint16_t result = 0;
	for (int i = 0; i < 4; i++)
	{
		result += L[i] << (i * 2);
		result += R[i] << (8 + i * 2);
	}
	for (int i = 0; i < numCliques; i++)
	{
		if (!has[i]) return nextSolution(leftSize, rightSize, result, numCliques);
	}
	return result;
}

uint16_t solveOneClique(int leftSize, int rightSize, int connections)
{
	for (int numCliques = 1; numCliques <= 4; numCliques++)
	{
		uint16_t firstSolution = getLexicographicallyFirstWithNumber(numCliques);
		if (isClique(leftSize, rightSize, firstSolution, connections)) return firstSolution;
		uint16_t solution = nextSolution(leftSize, rightSize, firstSolution, numCliques);
		while (solution != 0)
		{
			if (isClique(leftSize, rightSize, solution, connections)) return solution;
			solution = nextSolution(leftSize, rightSize, solution, numCliques);
		}
	}
	assert(false);
}

std::vector<uint16_t> presolvedCliques()
{
	std::vector<uint16_t> result;
	result.resize(4 * 4 * 65536);
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

std::vector<NodeClique> solveCliques(const std::vector<std::pair<size_t, bool>>& component, const GroupingBluntify& graph)
{
	static auto presolved = presolvedCliques();
	SizeArray<size_t, 4> left;
	SizeArray<size_t, 4> right;
	std::array<bool, 16> connected;
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
		result.back().right.push_back(right[0]);
		return result;
	}
	if (right.size() == 0)
	{
		assert(left.size() == 1);
		std::vector<NodeClique> result;
		result.emplace_back();
		result.back().left.push_back(left[0]);
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
	uint16_t solution = presolved[problemNumber];
	std::vector<NodeClique> result;
	result.resize(4);
	for (int i = 0; i < left.size(); i++)
	{
		auto cliqueNum = (solution >> (i * 2)) % 4;
		result[cliqueNum].left.push_back(left[i]);
	}
	for (int i = 0; i < right.size(); i++)
	{
		auto cliqueNum = (solution >> (8 + i * 2)) % 4;
		result[cliqueNum].right.push_back(right[i]);
	}
	while (result.size() > 0 && result.back().left.size() == 0 && result.back().right.size() == 0) result.pop_back();
	assert(result.size() > 0);
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
	auto size = std::max(grouping.items.back().leftCliques[grouping.items.back().leftCliques.size()-1], grouping.items.back().rightCliques[grouping.items.back().rightCliques.size()-1])+1;
	result.originalNodeAndOffset.resize(size, {0, -1});
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
			std::pair<size_t, int> nodeAndOffset { grouping.items[i].originalNodeId, grouping.items[i].originalNodeOffset + 1 };
			result.originalNodeAndOffset[clique] = nodeAndOffset;
		}
		for (size_t j = 0; j < grouping.items[i].leftCliques.size(); j++)
		{
			auto clique = grouping.items[i].leftCliques[j];
			assert(clique < size);
			std::pair<size_t, int> nodeAndOffset { grouping.items[i].originalNodeId, grouping.items[i].originalNodeOffset };
			result.originalNodeAndOffset[clique] = nodeAndOffset;
		}
	}
#ifndef NDEBUG
	for (size_t i = 0; i < result.originalNodeAndOffset.size(); i++)
	{
		assert(result.originalNodeAndOffset[i].second != -1);
	}
#endif
	return result;
}

GroupingBluntify loadUnbluntGFA(std::string filename)
{
	size_t resultSize = 0;
	size_t numNodes = 0;
	std::unordered_map<int, size_t> originalNodeStarts;
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
				assert(fromstart == "+");
				assert(toend == "+");
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
				originalNodeStarts[id] = resultSize;
				resultSize += seq.size() - result.overlap;
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
				result.items[pos].originalNodeId = id;
				result.items[pos].originalNodeOffset = 0;
				assert(seq.size() > result.overlap);
				for (size_t i = 1; i < seq.size() - result.overlap; i++)
				{
					result.items[pos+i-1].leftEdges.push_back(pos+i);
					result.items[pos+i].rightEdges.push_back(pos+i-1);
					result.items[pos+i].originalNodeId = id;
					result.items[pos+i].originalNodeOffset = i;
				}
				pos += seq.size() - result.overlap;
			}
			if (linestr[0] == 'L')
			{
				std::string dummy, fromstart, toend, overlap;
				int fromid, toid;
				line >> dummy >> fromid >> fromstart >> toid >> toend >> overlap;
				assert(fromstart == "+");
				assert(toend == "+");
				auto start = originalNodeEnds[fromid];
				auto end = originalNodeStarts[toid];
				bool found = false;
				for (size_t i = 0; i < result.items[start].leftEdges.size(); i++)
				{
					if (result.items[start].leftEdges[i] == end) found = true;
				}
				if (!found) result.items[start].leftEdges.push_back(end);
				found = false;
				for (size_t i = 0; i < result.items[end].rightEdges.size(); i++)
				{
					if (result.items[end].rightEdges[i] == start) found = true;
				}
				if (!found) result.items[end].rightEdges.push_back(start);
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
			auto node = graph.items[i].originalNodeId;
			auto pos = graph.items[i].originalNodeOffset;
			assert(originalSequences.count(node) == 1);
			assert(originalSequences[node].size() > pos);
			file << "S\t" << i << "\t" << originalSequences[node][pos] << std::endl;
			for (size_t j = 0; j < graph.items[i].leftEdges.size(); j++)
			{
				auto to = graph.items[i].leftEdges[j];
				file << "L\t" << i << "\t+\t" << to << "\t+\t" << graph.overlap << "M" << std::endl;
			}
		}
	}
}

int main(int argc, char** argv)
{
	auto graph = loadUnbluntGFA(argv[1]);
	while (graph.overlap > 0)
	{
		std::cerr << "overlap " << graph.overlap << " size " << graph.items.size() << std::endl;
		fillEdgeCliques(graph);
		auto cliques = getBluntedCliques(graph);
		graph = getSmallerGraph(cliques, graph.overlap);
	}
	writeResultGfa(graph, argv[1], argv[2]);
}