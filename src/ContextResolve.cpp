#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <cassert>
#include "Gfagraph.h"

std::vector<NodePos> getReversePath(const std::vector<NodePos>& path)
{
	auto reversePath = path;
	std::reverse(reversePath.begin(), reversePath.end());
	for (size_t i = 0; i < reversePath.size(); i++) reversePath[i] = reversePath[i].Reverse();
	return reversePath;
}

std::pair<NodePos, NodePos> canon(NodePos left, NodePos right)
{
	if (left.id == right.id)
	{
		if (!left.end && !right.end) return std::make_pair(right.Reverse(), left.Reverse());
		return std::make_pair(left, right);
	}
	if (left < right) return std::make_pair(left, right);
	assert(right.Reverse() < left.Reverse());
	return std::make_pair(right.Reverse(), left.Reverse());
}

bool comparePath(const std::vector<NodePos>& left, const std::vector<NodePos>& right)
{
	for (size_t i = 0; i < left.size() && i < right.size(); i++)
	{
		if (left[i].id < right[i].id) return true;
		if (left[i].id > right[i].id) return false;
		if (!left[i].end && right[i].end) return true;
		if (left[i].end && !right[i].end) return false;
	}
	return false;
}

class PathTree
{
public:
	std::vector<size_t> parent;
	std::vector<std::vector<size_t>> children;
	std::vector<NodePos> node;
	std::vector<size_t> count;
	void addPath(const std::vector<NodePos>& path, size_t offset)
	{
		if (node.size() == 0)
		{
			node.emplace_back(path[offset]);
			parent.emplace_back(0);
			children.emplace_back();
			count.emplace_back(0);
		}
		assert(path[offset] == node[0]);
		size_t currentNode = 0;
		count[0] += 1;
		for (size_t i = offset+1; i < path.size(); i++)
		{
			size_t nextNode = node.size();
			for (auto child : children[currentNode])
			{
				if (node[child] == path[i])
				{
					nextNode = child;
				}
			}
			if (nextNode == node.size())
			{
				node.emplace_back(path[i]);
				children[currentNode].emplace_back(nextNode);
				parent.emplace_back(currentNode);
				count.emplace_back(0);
				children.emplace_back();
			}
			currentNode = nextNode;
			count[currentNode] += 1;
		}
	}
};

void recursePaths(const PathTree& tree, std::vector<std::vector<NodePos>>& result, double currentCoverage, size_t minCoverage, size_t currentNode)
{
	if (currentCoverage < minCoverage) return;
	if (tree.children[currentNode].size() == 0)
	{
		result.emplace_back();
		while (currentNode != 0)
		{
			result.back().emplace_back(tree.node[currentNode]);
			currentNode = tree.parent[currentNode];
		}
		result.back().emplace_back(tree.node[0]);
		std::reverse(result.back().begin(), result.back().end());
		return;
	}
	size_t coverageSum = 0;
	for (auto child : tree.children[currentNode])
	{
		coverageSum += tree.count[child];
	}
	bool recursedAny = false;
	for (auto child : tree.children[currentNode])
	{
		double coverageHere = currentCoverage * (double)tree.count[child] / (double)coverageSum;
		if (coverageHere >= minCoverage)
		{
			recursePaths(tree, result, coverageHere, minCoverage, child);
			recursedAny = true;
		}
	}
	if (!recursedAny)
	{
		result.emplace_back();
		while (currentNode != 0)
		{
			result.back().emplace_back(tree.node[currentNode]);
			currentNode = tree.parent[currentNode];
		}
		result.back().emplace_back(tree.node[0]);
		std::reverse(result.back().begin(), result.back().end());
		return;
	}
}

std::vector<std::vector<NodePos>> getPaths(const PathTree& tree, size_t minCoverage)
{
	std::vector<std::vector<NodePos>> result;
	if (tree.count.size() == 0) return result;
	double currentCoverage = tree.count[0];
	recursePaths(tree, result, currentCoverage, minCoverage, 0);
	return result;
}

std::vector<std::vector<NodePos>> getDistinguishingPrefixes(const std::vector<std::vector<NodePos>>& paths)
{
	std::vector<std::vector<NodePos>> result = paths;
	std::sort(result.begin(), result.end(), comparePath);
	std::vector<size_t> longestPrefixMatch;
	longestPrefixMatch.resize(result.size(), 0);
	for (size_t i = 1; i < result.size(); i++)
	{
		size_t match = std::numeric_limits<size_t>::max();
		assert(result[i-1][0] == result[i][0]);
		for (size_t j = 1; j < result[i-1].size() && j < result[i].size(); j++)
		{
			if (result[i-1][j] != result[i][j])
			{
				match = j-1;
				break;
			}
		}
		assert(match != std::numeric_limits<size_t>::max());
		longestPrefixMatch[i-1] = std::max(longestPrefixMatch[i-1], match);
		longestPrefixMatch[i] = std::max(longestPrefixMatch[i], match);
	}
	if (result.size() > 1)
	{
		for (size_t i = 0; i < result.size(); i++)
		{
			assert(longestPrefixMatch[i]+2 <= result[i].size());
			result[i].erase(result[i].begin()+longestPrefixMatch[i]+2, result[i].end());
		}
	}
	return result;
}

class NodeTree
{
public:
	std::vector<std::vector<NodePos>> distinguishingPrefixes;
	std::vector<int> prefixComponent;
	template <typename Iter>
	std::vector<int> getPossibleComponents(Iter start, Iter end) const
	{
		std::vector<int> result;
		for (size_t prefixI = 0; prefixI < distinguishingPrefixes.size(); prefixI++)
		{
			bool match = true;
			Iter pos = start;
			for (size_t i = 0; i < distinguishingPrefixes[prefixI].size() && pos != end; i++, ++pos)
			{
				if (*pos != distinguishingPrefixes[prefixI][i])
				{
					match = false;
					break;
				}
			}
			if (match) result.emplace_back(prefixComponent[prefixI]);
		}
		return result;
	}
	template <typename Iter>
	int getUniqueComponent(Iter start, Iter end) const
	{
		int current = -1;
		for (size_t prefixI = 0; prefixI < distinguishingPrefixes.size(); prefixI++)
		{
			bool match = true;
			Iter pos = start;
			for (size_t i = 0; i < distinguishingPrefixes[prefixI].size() && pos != end; i++, ++pos)
			{
				if (*pos != distinguishingPrefixes[prefixI][i])
				{
					match = false;
					break;
				}
			}
			if (match)
			{
				if (current == -1)
				{
					current = prefixComponent[prefixI];
				}
				else
				{
					if (prefixComponent[prefixI] != current)
					{
						return -1;
					}
				}
			}
		}
		return current;
	}
};

class NodeTwowayTree
{
public:
	NodePos originalNode;
	NodeTree right;
	NodeTree left;
	int getUniqueComponent(const std::vector<NodePos>& path, const std::vector<NodePos>& reversePath, size_t offset) const
	{
		assert(path[offset] == originalNode);
		assert(reversePath[path.size() - offset - 1] == originalNode.Reverse());
		auto leftResult = left.getUniqueComponent(reversePath.begin() + (path.size() - offset - 1), reversePath.end());
		auto rightResult = right.getUniqueComponent(path.begin() + offset, path.end());
		if (leftResult == -1) return rightResult;
		if (rightResult == -1) return leftResult;
		if (leftResult != rightResult) return -1;
		return leftResult;
	}
	std::pair<std::vector<int>, std::vector<int>> getPossibleComponents(const std::vector<NodePos>& path, const std::vector<NodePos>& reversePath, size_t offset) const
	{
		assert(path[offset] == originalNode);
		assert(reversePath[path.size() - offset - 1] == originalNode.Reverse());
		auto leftResult = left.getPossibleComponents(reversePath.begin() + (path.size() - offset - 1), reversePath.end());
		auto rightResult = right.getPossibleComponents(path.begin() + offset, path.end());
		return std::make_pair(leftResult, rightResult);
	}
};

std::vector<int> merge(std::vector<int> leftResult, std::vector<int> rightResult)
{
	std::sort(leftResult.begin(), leftResult.end());
	std::sort(rightResult.begin(), rightResult.end());
	std::vector<int> result;
	size_t leftPos = 0;
	size_t rightPos = 0;
	while (leftPos < leftResult.size() || rightPos < rightResult.size())
	{
		if (leftPos == leftResult.size())
		{
			result.emplace_back(rightResult[rightPos]);
			rightPos++;
			continue;
		}
		if (rightPos == rightResult.size())
		{
			result.emplace_back(leftResult[leftPos]);
			leftPos++;
			continue;
		}
		if (leftResult[leftPos] < rightResult[rightPos])
		{
			result.emplace_back(leftResult[leftPos]);
			leftPos++;
			continue;
		}
		if (rightResult[rightPos] < leftResult[leftPos])
		{
			result.emplace_back(rightResult[rightPos]);
			rightPos++;
			continue;
		}
		assert(leftResult[leftPos] == rightResult[rightPos]);
		result.emplace_back(rightResult[rightPos]);
		rightPos++;
		leftPos++;
		continue;
	}
	return result;
}

std::vector<std::vector<NodePos>> rejigPaths(const std::vector<std::vector<NodePos>>& paths, const std::vector<NodeTwowayTree>& trees)
{
	std::vector<std::vector<NodePos>> result;
	size_t solved = 0;
	size_t unsolved = 0;
	for (auto path : paths)
	{
		result.emplace_back();
		auto reversePath = getReversePath(path);
		for (size_t i = 0; i < path.size(); i++)
		{
			int nodeId = 0;
			if (!path[i].end)
			{
				nodeId = trees[path[i].id].getUniqueComponent(reversePath, path, path.size()-i-1);
			}
			else
			{
				nodeId = trees[path[i].id].getUniqueComponent(path, reversePath, i);
			}
			result.back().emplace_back(nodeId, path[i].end);
			if (nodeId != -1)
			{
				solved += 1;
			}
			else
			{
				unsolved += 1;
			}
		}
	}
	std::cout << "solved " << solved << ", unsolved " << unsolved << std::endl;
	return result;
}

std::vector<std::pair<NodePos, NodePos>> getComponentEdges(const std::vector<std::vector<NodePos>>& paths, const size_t minCoverage)
{
	std::unordered_map<std::pair<NodePos, NodePos>, size_t> edgeCoverage;
	for (auto path : paths)
	{
		for (size_t i = 1; i < path.size(); i++)
		{
			if (path[i].id == -1 || path[i-1].id == -1) continue;
			auto edge = canon(path[i-1], path[i]);
			edgeCoverage[edge] += 1;
		}
	}
	std::vector<std::pair<NodePos, NodePos>> result;
	for (auto pair : edgeCoverage)
	{
		if (pair.second < minCoverage) continue;
		result.push_back(pair.first);
	}
	std::cout << "edges " << result.size() << std::endl;
	return result;
}

void writeGraph(const std::string& fileName, const std::vector<NodeTwowayTree>& trees, const std::vector<std::pair<NodePos, NodePos>>& edges, const GfaGraph& oldGraph)
{
	std::ofstream file { fileName };
	std::unordered_set<size_t> outputted;
	for (size_t i = 0; i < trees.size(); i++)
	{
		assert(trees[i].originalNode.end);
		if (oldGraph.nodes.count(trees[i].originalNode.id) == 0) continue;
		std::string sequence = oldGraph.nodes.at(trees[i].originalNode.id);
		for (auto id : trees[i].left.prefixComponent)
		{
			if (outputted.count(id) == 1) continue;
			outputted.insert(id);
			file << "S\t" << id << "\t" << sequence << std::endl;
		}
		for (auto id : trees[i].right.prefixComponent)
		{
			if (outputted.count(id) == 1) continue;
			outputted.insert(id);
			file << "S\t" << id << "\t" << sequence << std::endl;
		}
	}
	std::cout << "nodes " << outputted.size() << std::endl;
	for (auto edge : edges)
	{
		file << "L\t" << edge.first.id << "\t" << (edge.first.end ? "+" : "-") << "\t" << edge.second.id << "\t" << (edge.second.end ? "+" : "-") << "\t0M" << std::endl;
	}
}

std::vector<std::vector<NodePos>> loadGaf(const std::string& gafFile)
{
	std::ifstream file { gafFile };
	std::vector<std::vector<NodePos>> result;
	while (file.good())
	{
		std::string str;
		std::getline(file, str);
		if (!file.good()) break;
		std::stringstream parseStr { str };
		std::string field;
		for (size_t i = 0; i < 6; i++)
		{
			std::getline(parseStr, field, '\t');
		}
		std::stringstream parseField { field };
		result.emplace_back();
		for (size_t i = 0; i < field.size(); i++)
		{
			if (field[i] == '>')
			{
				int nodeID = std::stoi(field.data()+i+1);
				result.back().emplace_back(nodeID, true);
			}
			else if (field[i] == '<')
			{
				int nodeID = std::stoi(field.data()+i+1);
				result.back().emplace_back(nodeID, false);
			}
		}
		assert(result.back().size() != 0);
	}
	std::cout << "paths " << result.size() << std::endl;
	return result;
}

std::vector<PathTree> getPathTrees(const std::vector<std::vector<NodePos>>& paths)
{
	std::vector<PathTree> pathTrees;
	int maxNode = 0;
	for (auto path : paths)
	{
		for (auto node : path)
		{
			maxNode = std::max(maxNode, node.id);
		}
	}
	pathTrees.resize(maxNode*2+2);
	for (auto path : paths)
	{
		auto reversePath = getReversePath(path);
		for (size_t j = 0; j < path.size(); j++)
		{
			pathTrees[path[j].id * 2 + (path[j].end ? 1 : 0)].addPath(path, j);
			pathTrees[reversePath[j].id * 2 + (reversePath[j].end ? 1 : 0)].addPath(reversePath, j);
		}
	}
	return pathTrees;
}

std::vector<NodeTwowayTree> getNodeTrees(const std::vector<PathTree>& pathTrees, size_t minCoverage)
{
	std::vector<NodeTwowayTree> result;
	result.resize(pathTrees.size() / 2);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].originalNode.id = i;
		result[i].originalNode.end = true;
		result[i].right.distinguishingPrefixes = getDistinguishingPrefixes(getPaths(pathTrees[i*2+1], minCoverage));
		result[i].left.distinguishingPrefixes = getDistinguishingPrefixes(getPaths(pathTrees[i*2], minCoverage));
		assert(result[i].right.distinguishingPrefixes.size() > 0 == result[i].left.distinguishingPrefixes.size() > 0);
	}
	return result;
}

size_t find(std::vector<size_t>& parent, size_t item)
{
	if (parent[item] == item) return item;
	auto result = find(parent, parent[item]);
	parent[item] = result;
	return result;
}

void merge(std::vector<size_t>& parent, std::vector<size_t>& rank, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	if (rank[left] < rank[right]) std::swap(left, right);
	parent[right] = left;
	if (rank[right] == rank[left]) rank[left] += 1;
}

void addComponents(std::vector<NodeTwowayTree>& trees, const std::vector<std::vector<NodePos>>& paths, size_t minCoverage)
{
	int currentComponent = 0;
	for (auto& tree : trees)
	{
		tree.left.prefixComponent.resize(tree.left.distinguishingPrefixes.size());
		for (size_t i = 0; i < tree.left.prefixComponent.size(); i++)
		{
			tree.left.prefixComponent[i] = currentComponent;
			currentComponent++;
		}
		tree.right.prefixComponent.resize(tree.right.distinguishingPrefixes.size());
		for (size_t i = 0; i < tree.right.prefixComponent.size(); i++)
		{
			tree.right.prefixComponent[i] = currentComponent;
			currentComponent++;
		}
	}
	std::vector<std::vector<std::vector<int>>> hyperEdges;
	hyperEdges.resize(trees.size());
	for (auto path : paths)
	{
		auto reversePath = getReversePath(path);
		for (size_t i = 0; i < path.size(); i++)
		{
			std::pair<std::vector<int>, std::vector<int>> possibles;
			if (trees[path[i].id].left.prefixComponent.size() == 0)
			{
				assert(trees[path[i].id].right.prefixComponent.size() == 0);
				continue;
			}
			if (!path[i].end)
			{
				possibles = trees[path[i].id].getPossibleComponents(reversePath, path, path.size()-i-1);
			}
			else
			{
				possibles = trees[path[i].id].getPossibleComponents(path, reversePath, i);
			}
			if (possibles.first.size() == 0 || possibles.second.size() == 0) continue;
			hyperEdges[path[i].id].push_back(merge(possibles.first, possibles.second));
			std::sort(hyperEdges[path[i].id].back().begin(), hyperEdges[path[i].id].back().end());
			assert(hyperEdges[path[i].id].back().size() >= 2);
		}
	}
	std::vector<size_t> parent;
	std::vector<size_t> rank;
	parent.resize(currentComponent);
	rank.resize(currentComponent);
	for (size_t i = 0; i < currentComponent; i++)
	{
		parent[i] = i;
		rank[i] = 0;
	}
	for (size_t i = 0; i < hyperEdges.size(); i++)
	{
		std::sort(hyperEdges[i].begin(), hyperEdges[i].end(), [](const std::vector<int>& left, const std::vector<int>& right) 
		{
			for (size_t i = 0; i < left.size() && i < right.size(); i++)
			{
				if (left[i] < right[i]) return true;
				if (right[i] < left[i]) return false;
			}
			return false;
		});
		std::stable_sort(hyperEdges[i].begin(), hyperEdges[i].end(), [](const std::vector<int>& left, const std::vector<int>& right) { return left.size() < right.size(); });
		std::unordered_set<int> hasEdge;
		std::unordered_set<int> getsEdge;
		size_t currentCount = 1;
		for (size_t j = 1; j < hyperEdges[i].size(); j++)
		{
			if (hyperEdges[i][j] != hyperEdges[i][j-1])
			{
				if (currentCount >= minCoverage)
				{
					getsEdge.insert(hyperEdges[i][j-1].begin(), hyperEdges[i][j-1].end());
					std::vector<int> canMerge;
					for (auto k : hyperEdges[i][j-1])
					{
						if (hasEdge.count(k) == 0) canMerge.push_back(k);
					}
					for (size_t k = 1; k < canMerge.size(); k++)
					{
						merge(parent, rank, canMerge[0], canMerge[k]);
					}
				}
				currentCount = 0;
			}
			if (hyperEdges[i][j].size() != hyperEdges[i][j-1].size())
			{
				hasEdge.insert(getsEdge.begin(), getsEdge.end());
				getsEdge.clear();
			}
			currentCount += 1;
		}
		if (hyperEdges[i].size() > 0 && currentCount >= minCoverage)
		{
			std::vector<int> canMerge;
			for (auto k : hyperEdges[i].back())
			{
				if (hasEdge.count(k) == 0) canMerge.push_back(k);
			}
			for (size_t k = 1; k < canMerge.size(); k++)
			{
				merge(parent, rank, canMerge[0], canMerge[k]);
			}
		}
	}
	for (auto& tree : trees)
	{
		for (size_t i = 0; i < tree.left.prefixComponent.size(); i++)
		{
			tree.left.prefixComponent[i] = find(parent, tree.left.prefixComponent[i]);
		}
		for (size_t i = 0; i < tree.right.prefixComponent.size(); i++)
		{
			tree.right.prefixComponent[i] = find(parent, tree.right.prefixComponent[i]);
		}
	}
}

std::vector<NodeTwowayTree> solveTrees(const std::vector<std::vector<NodePos>>& paths, size_t minCoverage)
{
	auto pathTrees = getPathTrees(paths);
	auto nodeTrees = getNodeTrees(pathTrees, minCoverage);
	// addComponents(nodeTrees, paths, 1);
	addComponents(nodeTrees, paths, minCoverage);
	return nodeTrees;
}

int main(int argc, char** argv)
{
	std::string graphFile = argv[1];
	std::string gafFile = argv[2];
	std::string outputFile = argv[3];

	auto paths = loadGaf(gafFile);
	auto trees = solveTrees(paths, 2);
	auto resolvedPaths = rejigPaths(paths, trees);
	auto componentEdges = getComponentEdges(resolvedPaths, 2);
	auto graph = GfaGraph::LoadFromFile(graphFile);
	writeGraph(outputFile, trees, componentEdges, graph);
}
