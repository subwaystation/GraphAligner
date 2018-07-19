#include <algorithm>
#include <fstream>
#include <unordered_map>
#include "GfaGraph.h"
#include "CommonUtils.h"

namespace std {
	template <>
	struct hash<std::pair<int, bool>>
	{
		std::size_t operator()(const std::pair<int, bool>& k) const
		{
			return std::hash<size_t>()(k.first) ^ std::hash<bool>()(k.second);
		}
	};
}

struct PathSuffixNode
{
	int nodeId;
	bool reverse;
	size_t parent;
	std::unordered_map<std::pair<int, bool>, size_t> children;
	bool equivalenceFlag;
	size_t equivalentNode;
};

struct GraphNode
{
	std::vector<std::pair<int, bool>> originalNodes;
};

struct Graph
{
	std::vector<GraphNode> nodes;
	std::vector<std::pair<size_t, size_t>> edges;
};

std::vector<std::vector<std::pair<int, bool>>> getAlignments(std::string filename)
{
	auto alns = CommonUtils::LoadVGAlignments(filename);
	std::vector<std::vector<std::pair<int, bool>>> result;
	for (auto aln : alns)
	{
		result.emplace_back();
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			result.back().emplace_back(aln.path().mapping(i).position().node_id(), aln.path().mapping(i).position().is_reverse());
		}
		//also add reverse complement
		result.emplace_back();
		for (int i = aln.path().mapping_size()-1; i >= 0; i--)
		{
			result.back().emplace_back(aln.path().mapping(i).position().node_id(), !aln.path().mapping(i).position().is_reverse());
		}
	}
	return result;
}

void addSuffix(std::vector<PathSuffixNode>& tree, const std::vector<std::pair<int, bool>>& aln, size_t start)
{
	size_t alnpos = start;
	size_t treepos = 0;
	while (alnpos < aln.size() && tree[treepos].children.count(aln[alnpos]) == 1)
	{
		treepos = tree[treepos].children[aln[alnpos]];
		alnpos++;
	}
	for (;alnpos < aln.size(); alnpos++)
	{
		tree.emplace_back();
		tree.back().nodeId = aln[alnpos].first;
		tree.back().reverse = aln[alnpos].second;
		tree.back().parent = treepos;
		tree[treepos].children[aln[alnpos]] = tree.size()-1;
		tree.back().equivalenceFlag = false;
		tree.back().equivalentNode = 0;
		treepos = tree.size()-1;
	}
}

void verifyTreeTopology(const std::vector<PathSuffixNode>& tree)
{
	for (size_t i = 1; i < tree.size(); i++)
	{
		assert(tree[tree[i].parent].children.count(std::make_pair(tree[i].nodeId, tree[i].reverse)) == 1);
		assert(tree[tree[i].parent].children.at(std::make_pair(tree[i].nodeId, tree[i].reverse)) == i);
	}
	std::vector<bool> reachableFromRoot;
	reachableFromRoot.resize(tree.size(), false);
	reachableFromRoot[0] = true;
	for (size_t i = 0; i < tree.size(); i++)
	{
		assert(reachableFromRoot[i]);
		for (auto child : tree[i].children)
		{
			reachableFromRoot[child.second] = true;
		}
	}
}

void verifyTreeEquivalences(const std::vector<PathSuffixNode>& tree)
{
	for (size_t i = 1; i < tree.size(); i++)
	{
		if (tree[i].equivalentNode != 0)
		{
			assert(tree[i].children.size() == 0);
			size_t pos = i;
			std::vector<std::pair<int, bool>> stackHere;
			while (pos != 0)
			{
				stackHere.emplace_back(tree[pos].nodeId, tree[pos].reverse);
				pos = tree[pos].parent;
			}
			std::vector<std::pair<int, bool>> stackThere;
			pos = tree[i].equivalentNode;
			while (pos != 0)
			{
				stackThere.emplace_back(tree[pos].nodeId, tree[pos].reverse);
				pos = tree[pos].parent;
			}
			assert(stackThere.size() > 0);
			assert(stackHere.size() == stackThere.size() + 1);
			for (size_t i = 0; i < stackThere.size(); i++)
			{
				assert(stackHere[i] == stackThere[i]);
			}
		}
		else
		{
			assert(tree[i].parent == 0 || tree[i].children.size() > 0);
		}
	}
}

std::vector<PathSuffixNode> getPathSuffixTree(const std::vector<std::vector<std::pair<int, bool>>>& alns)
{
	std::vector<PathSuffixNode> result;
	result.emplace_back();
	result[0].nodeId = 0;
	result[0].reverse = false;
	result[0].parent = 0;
	result[0].equivalentNode = 0;
	result[0].equivalenceFlag = false;
	for (auto aln : alns)
	{
		for (size_t i = aln.size()-1; i < aln.size(); i--)
		{
			addSuffix(result, aln, i);
		}
	}
	verifyTreeTopology(result);
	for (size_t i = 0; i < result.size(); i++)
	{
		if (result[i].children.size() != 0) continue;
		std::vector<std::pair<size_t, bool>> stack;
		size_t treepos = i;
		while (treepos != 0)
		{
			stack.emplace_back(result[treepos].nodeId, result[treepos].reverse);
			treepos = result[treepos].parent;
		}
		stack.pop_back();
		std::reverse(stack.begin(), stack.end());
		size_t equivalentTreePos = 0;
		for (size_t i = 0; i < stack.size(); i++)
		{
			assert(result[equivalentTreePos].children.count(stack[i]) == 1);
			equivalentTreePos = result[equivalentTreePos].children[stack[i]];
		}
		if (equivalentTreePos != 0)
		{
			result[i].equivalentNode = equivalentTreePos;
		}
	}
	verifyTreeEquivalences(result);
	// for (size_t i = 0; i < result.size(); i++)
	// {
	// 	if (result[i].equivalentNode != 0)
	// 	{
	// 		size_t posHere = i;
	// 		size_t posThere = result[i].equivalentNode;
	// 		assert(posHere != 0);
	// 		assert(posThere != 0);
	// 		while (result[result[posHere].parent].children.size() == 1 && result[result[posThere].parent].children.size() == 1)
	// 		{
	// 			posHere = result[posHere].parent;
	// 			posThere = result[posThere].parent;
	// 		}
	// 		result[i].equivalentNode = 0;
	// 		result[posHere].equivalentNode = posThere;
	// 	}
	// }
	for (size_t i = 0; i < result.size(); i++)
	{
		result[result[i].equivalentNode].equivalenceFlag = true;
	}
	return result;
}

void addGraphNodesRec(const std::vector<PathSuffixNode>& tree, size_t start, size_t current, std::vector<size_t> stack, Graph& result, std::vector<std::vector<size_t>>& graphNodesFromTree)
{
	stack.push_back(current);
	if (tree[current].children.size() == 0)
	{
		if (graphNodesFromTree[stack.back()].size() > 0) return;
		graphNodesFromTree[start].push_back(result.nodes.size());
		result.nodes.emplace_back();
		result.nodes.back().originalNodes.emplace_back(tree[start].nodeId, tree[start].reverse);
		for (size_t i = 1; i < stack.size(); i++)
		{
			graphNodesFromTree[stack[i]].push_back(result.nodes.size());
			result.nodes.emplace_back();
			result.nodes.back().originalNodes.emplace_back(tree[stack[i]].nodeId, tree[stack[i]].reverse);
			assert(result.nodes.size() >= 2);
			result.edges.emplace_back(result.nodes.size()-2, result.nodes.size()-1);
		}
	}
	else
	{
		for (auto child : tree[current].children)
		{
			addGraphNodesRec(tree, start, child.second, stack, result, graphNodesFromTree);
		}
	}
}

void addPathsWithoutEquivalenceRec(size_t pos, std::vector<size_t> stack, const std::vector<PathSuffixNode>& tree, Graph& result, std::vector<std::vector<size_t>>& graphNodesFromTree)
{
	stack.push_back(pos);
	if (tree[pos].children.size() == 0)
	{
		if (tree[stack.back()].equivalenceFlag || graphNodesFromTree[stack.back()].size() != 0) return;
		result.nodes.emplace_back();
		result.nodes.back().originalNodes.emplace_back(tree[stack.back()].nodeId, tree[stack.back()].reverse);
		graphNodesFromTree[stack.back()].push_back(result.nodes.size()-1);
		for (size_t i = stack.size()-2; i < stack.size(); i--)
		{
			if (tree[stack[i]].equivalenceFlag || graphNodesFromTree[stack[i]].size() != 0)
			{
				for (auto node : graphNodesFromTree[stack[i]])
				{
					result.edges.emplace_back(node, result.nodes.size()-1);
				}
				return;
			}
			result.nodes.emplace_back();
			result.nodes.back().originalNodes.emplace_back(tree[stack[i]].nodeId, tree[stack[i]].reverse);
			//edge this way because nodes are added in the reverse order
			result.edges.emplace_back(result.nodes.size()-1, result.nodes.size()-2);
			graphNodesFromTree[stack[i]].push_back(result.nodes.size()-1);
		}
	}
	else
	{
		for (auto child : tree[pos].children)
		{
			addPathsWithoutEquivalenceRec(child.second, stack, tree, result, graphNodesFromTree);
		}
	}
}

void addPathsWithoutEquivalence(const std::vector<PathSuffixNode>& tree, Graph& result, std::vector<std::vector<size_t>>& graphNodesFromTree)
{
	for (auto child : tree[0].children)
	{
		std::vector<size_t> stack;
		addPathsWithoutEquivalenceRec(child.second, stack, tree, result, graphNodesFromTree);
	}
}

Graph getGraph(const std::vector<PathSuffixNode>& tree)
{
	std::set<size_t> expandEdges;
	std::set<size_t> expandGraphNodes;
	std::vector<std::vector<size_t>> graphNodesFromTree;
	graphNodesFromTree.resize(tree.size());
	Graph result;
	for (size_t i = 0; i < tree.size(); i++)
	{
		if (tree[i].children.size() == 0 && !tree[i].equivalenceFlag && tree[i].equivalentNode != 0)
		{
			expandEdges.insert(i);
			expandGraphNodes.insert(tree[i].equivalentNode);
		}
	}
	for (auto node : expandGraphNodes)
	{
		for (auto child : tree[node].children)
		{
			std::vector<size_t> stack;
			addGraphNodesRec(tree, child.second, child.second, stack, result, graphNodesFromTree);
		}
	}
	addPathsWithoutEquivalence(tree, result, graphNodesFromTree);
	for (auto node : expandEdges)
	{
		assert(graphNodesFromTree[node].size() == 1);
		size_t treepos = node;
		while (tree[treepos].children.size() == 0)
		{
			treepos = tree[treepos].equivalentNode;
		}
		if (treepos == 0) continue;
		for (auto child : tree[treepos].children)
		{
			for (auto source : graphNodesFromTree[node])
			{
				for (auto target : graphNodesFromTree[child.second])
				{
					result.edges.emplace_back(source, target);
				}
			}
		}
	}
	return result;
}

Graph mergeGraphUnitigs(const Graph& raw)
{
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<std::vector<size_t>> inNeighbors;
	outNeighbors.resize(raw.nodes.size());
	inNeighbors.resize(raw.nodes.size());
	for (auto edge : raw.edges)
	{
		outNeighbors[edge.first].push_back(edge.second);
		inNeighbors[edge.second].push_back(edge.first);
	}
	std::vector<size_t> unitigMapping;
	unitigMapping.resize(raw.nodes.size(), std::numeric_limits<size_t>::max());
	Graph result;
	for (size_t i = 0; i < raw.nodes.size(); i++)
	{
		if (inNeighbors[i].size() == 1 && outNeighbors[inNeighbors[i][0]].size() == 1) continue;
		size_t unitig = result.nodes.size();
		result.nodes.emplace_back();
		size_t pos = i;
		while (true)
		{
			assert(unitigMapping[pos] == std::numeric_limits<size_t>::max());
			unitigMapping[pos] = unitig;
			assert(raw.nodes[pos].originalNodes.size() == 1);
			result.nodes[unitig].originalNodes.push_back(raw.nodes[pos].originalNodes[0]);
			if (outNeighbors[pos].size() != 1) break;
			pos = outNeighbors[pos][0];
			if (inNeighbors[pos].size() != 1) break;
		}
	}
	for (size_t i = 0; i < raw.nodes.size(); i++)
	{
		assert(unitigMapping[i] != std::numeric_limits<size_t>::max());
	}
	for (auto edge : raw.edges)
	{
		assert(unitigMapping[edge.first] != std::numeric_limits<size_t>::max());
		assert(unitigMapping[edge.second] != std::numeric_limits<size_t>::max());
		if (unitigMapping[edge.first] == unitigMapping[edge.second]) continue;
		result.edges.emplace_back(unitigMapping[edge.first], unitigMapping[edge.second]);
	}
	return result;
}

void writeGraph(const Graph& graph, const GfaGraph& originalGfa, std::string outputFileName)
{
	std::ofstream output { outputFileName };
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		std::string nodeLabel;
		assert(graph.nodes[i].originalNodes.size() > 0);
		if (graph.nodes[i].originalNodes[0].second)
		{
			nodeLabel += CommonUtils::ReverseComplement(originalGfa.nodes.at(graph.nodes[i].originalNodes[0].first)).substr(0, originalGfa.edgeOverlap);
		}
		else
		{
			nodeLabel += originalGfa.nodes.at(graph.nodes[i].originalNodes[0].first).substr(0, originalGfa.edgeOverlap);
		}
		for (auto original : graph.nodes[i].originalNodes)
		{
			if (original.second)
			{
				nodeLabel += CommonUtils::ReverseComplement(originalGfa.nodes.at(original.first)).substr(originalGfa.edgeOverlap);
			}
			else
			{
				nodeLabel += originalGfa.nodes.at(original.first).substr(originalGfa.edgeOverlap);
			}
		}
		output << "S\t" << i << "\t" << nodeLabel << std::endl;
	}
	for (auto edge : graph.edges)
	{
		output << "L\t" << edge.first << "\t+\t" << edge.second << "\t+\t" << originalGfa.edgeOverlap << "M" << std::endl;
	}
}

int main(int argc, char** argv)
{
	std::string inputAlignmentFile { argv[1] };
	std::string inputGraphFile { argv[2] };
	std::string outputGraphFile { argv[3] };

	auto alns = getAlignments(inputAlignmentFile);
	auto tree = getPathSuffixTree(alns);
	auto graph = getGraph(tree);
	auto graphUnitigs = mergeGraphUnitigs(graph);
	auto gfaGraph = GfaGraph::LoadFromFile(inputGraphFile);
	// writeGraph(graph, gfaGraph, outputGraphFile);
	writeGraph(graphUnitigs, gfaGraph, outputGraphFile);
}
