#include <fstream>
#include "stream.hpp"
#include "GfaGraph.h"
#include "CommonUtils.h"

namespace std {
	template <> struct hash<std::pair<int, int>>
	{
		size_t operator()(const std::pair<int, int>& pair) const
		{
			return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
		}
	};
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

std::unordered_set<int> getSafeChains(const GfaGraph& graph, size_t safeChainSize, const std::unordered_map<size_t, size_t>& coveragePerNode, double averageCoverage)
{
	std::unordered_map<int, size_t> chainSize;
	std::unordered_map<int, size_t> chainKmers;
	size_t totalKmers = 0;
	size_t totalSize = 0;
	for (auto tags : graph.tags)
	{
		size_t size = 0;
		int chain = 0;
		size_t kmers = 0;
		bool hasKmers = true;
		bool hasSize = false;
		bool hasChain = false;
		std::stringstream sstr { tags.second };
		if (coveragePerNode.count(tags.first) == 1)
		{
			hasKmers = true;
			kmers = coveragePerNode.at(tags.first);
		}		
		while (sstr.good())
		{
			std::string tag;
			sstr >> tag;
			// if (tag.substr(0, 5) == "RC:i:")
			// {
			// 	assert(!hasKmers);
			// 	hasKmers = true;
			// 	kmers = std::stol(tag.substr(5));
			// }
			if (tag.substr(0, 5) == "LN:i:")
			{
				assert(!hasSize);
				hasSize = true;
				size = std::stol(tag.substr(5));
			}
			if (tag.substr(0, 5) == "bc:Z:")
			{
				assert(!hasChain);
				hasChain = true;
				chain = std::stoi(tag.substr(5));
			}
		}
		if (!hasSize || !hasChain) continue;
		chainSize[chain] += size;
		chainKmers[chain] += kmers;
		totalKmers += kmers;
		totalSize += size;
	}
	std::cerr << "average coverage " << averageCoverage << std::endl;
	std::unordered_set<int> result;
	for (auto pair : chainSize)
	{
		if (pair.second < safeChainSize) continue;
		assert(chainKmers.count(pair.first) == 1);
		double coverage = (double)chainKmers.at(pair.first) / (double)pair.second;
		if (coverage < averageCoverage * 0.6) continue;
		if (coverage > averageCoverage * 1.4) continue;
		result.insert(pair.first);
	}
	return result;
}

std::unordered_map<int, int> getChainBelongers(const GfaGraph& graph)
{
	std::unordered_map<int, int> result;
	for (auto tags : graph.tags)
	{
		std::stringstream sstr { tags.second };
		while (sstr.good())
		{
			std::string tag;
			sstr >> tag;
			if (tag.substr(0, 5) == "bc:Z:")
			{
				int chain = std::stoi(tag.substr(5));
				result[tags.first] = chain;
			}
		}
	}
	return result;
}

struct ResolvableComponent
{
	std::unordered_set<int> nodeIDs;
	std::unordered_set<std::pair<NodePos, NodePos>> edges;
	std::unordered_set<int> removedNodes;
	std::unordered_set<std::pair<NodePos, NodePos>> removedEdges;
	std::unordered_map<int, NodePos> newNodes;
	std::vector<std::pair<NodePos, NodePos>> newEdges;
};

struct Path
{
	std::string name;
	std::vector<NodePos> path;
};

struct Subpath
{
	size_t index;
	size_t start;
	std::vector<NodePos> path;
	std::vector<NodePos> replacement;
};

template <typename T>
T find(std::unordered_map<T, T>& parent, T k)
{
	std::vector<T> fix;
	if (parent.count(k) == 0)
	{
		parent[k] = k;
	}
	while (parent.at(k) != k)
	{
		fix.push_back(k);
		k = parent.at(k);
	}
	for (auto i : fix)
	{
		parent[i] = k;
	}
	return k;
}

template <typename T>
void merge(std::unordered_map<T, T>& parent, T k1, T k2)
{
	parent[find(parent, k2)] = find(parent, k1);
}

std::vector<ResolvableComponent> getComponents(const GfaGraph& graph, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	std::unordered_map<std::pair<NodePos, NodePos>, std::pair<NodePos, NodePos>> parent;
	for (auto edge : graph.edges)
	{
		auto firstKey = std::make_pair(edge.first, edge.second[0]);
		for (auto target : edge.second)
		{
			// either all out-edges lead to the same chain or all of them lead to a different one
			if (belongsToChain.at(edge.first.id) == belongsToChain.at(target.id) && safeChains.count(belongsToChain.at(edge.first.id)) == 1) break;
			auto keyHere = std::make_pair(edge.first, target);
			auto reverseKeyHere = std::make_pair(target.Reverse(), edge.first.Reverse());
			merge(parent, firstKey, keyHere);
			merge(parent, keyHere, reverseKeyHere);
		}
		if (safeChains.count(belongsToChain.at(edge.first.id)) == 0)
		{
			if (graph.edges.count(edge.first.Reverse()) == 1)
			{
				for (auto target : graph.edges.at(edge.first.Reverse()))
				{
					// either all out-edges lead to the same chain or all of them lead to a different one
					if (belongsToChain.at(edge.first.id) == belongsToChain.at(target.id) && safeChains.count(belongsToChain.at(edge.first.id)) == 1) break;
					auto keyHere = std::make_pair(edge.first.Reverse(), target);
					auto reverseKeyHere = std::make_pair(target.Reverse(), edge.first);
					merge(parent, firstKey, keyHere);
					merge(parent, keyHere, reverseKeyHere);
				}
			}
		}
	}
	std::vector<std::pair<NodePos, NodePos>> keys;
	for (auto pair : parent)
	{
		keys.push_back(pair.second);
	}
	for (auto key : keys)
	{
		find(parent, key);
	}
	std::unordered_map<std::pair<NodePos, NodePos>, size_t> componentNumber;
	int nextComponent = 1;
	for (auto pair : parent)
	{
		if (componentNumber.count(pair.second) == 0)
		{
			componentNumber[pair.second] = nextComponent;
			nextComponent += 1;
		}
	}
	std::vector<ResolvableComponent> result;
	result.resize(nextComponent);
	for (auto edge : graph.edges)
	{
		for (auto target : edge.second)
		{
			auto keyHere = std::make_pair(edge.first, target);
			auto componentNum = componentNumber[find(parent, keyHere)];
			if (componentNum == 0) continue;
			componentNum -= 1;
			result[componentNum].edges.insert(keyHere);
			if (safeChains.count(belongsToChain.at(edge.first.id)) == 0) result[componentNum].nodeIDs.insert(edge.first.id);
			if (safeChains.count(belongsToChain.at(target.id)) == 0) result[componentNum].nodeIDs.insert(target.id);
		}
	}
	while (result.size() > 0 && result[0].edges.size() == 0)
	{
		std::swap(result[0], result.back());
		result.pop_back();
	}
	for (size_t i = 1; i < result.size(); i++)
	{
		if (result[i].edges.size() == 0)
		{
			std::swap(result[i], result.back());
			result.pop_back();
			i--;
		}
	}
	return result;
}

std::tuple<std::vector<Path>, std::unordered_map<size_t, size_t>, size_t> loadAlignmentsAsPaths(std::string filename)
{
	std::vector<Path> result;
	std::ifstream file { filename, std::ios::in | std::ios::binary };
	size_t bpCount = 0;
	std::unordered_map<size_t, size_t> coveragePerNode;
	std::function<void(vg::Alignment&)> lambda = [&result, &bpCount, &coveragePerNode](vg::Alignment& g) {
		result.emplace_back();
		result.back().name = g.name();
		bpCount += g.sequence().size();
		for (int i = 0; i < g.path().mapping_size(); i++)
		{
			result.back().path.emplace_back(g.path().mapping(i).position().node_id(), !g.path().mapping(i).position().is_reverse());
			for (size_t j = 0; j < g.path().mapping(i).edit_size(); j++)
			{
				coveragePerNode[g.path().mapping(i).position().node_id()] += g.path().mapping(i).edit(j).to_length();
			}
		}
	};
	stream::for_each(file, lambda);
	return std::make_tuple(result, coveragePerNode, bpCount);
}

std::vector<std::vector<Subpath>> splitPathsPerComponent(const std::vector<Path>& paths, const std::vector<ResolvableComponent>& components, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	std::vector<std::vector<Subpath>> result;
	result.resize(components.size());
	std::unordered_map<std::pair<NodePos, NodePos>, size_t> edgeBelongsToComponent;
	for (size_t i = 0; i < components.size(); i++)
	{
		for (auto edge : components[i].edges)
		{
			edgeBelongsToComponent[edge] = i;
		}
	}
	for (size_t path = 0; path < paths.size(); path++)
	{
		if (paths[path].path.size() < 2) continue;
		auto currentEdge = std::make_pair(paths[path].path[0], paths[path].path[1]);
		size_t previousBelongedToComponent = components.size();
		if (edgeBelongsToComponent.count(currentEdge) == 1) previousBelongedToComponent = edgeBelongsToComponent.at(currentEdge);
		Subpath currentPath;
		currentPath.index = path;
		currentPath.start = 0;
		currentPath.path.push_back(paths[path].path[0]);
		for (size_t i = 1; i < paths[path].path.size(); i++)
		{
			currentEdge = std::make_pair(paths[path].path[i-1], paths[path].path[i]);
			bool split = false;
			if (edgeBelongsToComponent.count(currentEdge) == 0 || edgeBelongsToComponent.at(currentEdge) != previousBelongedToComponent) split = true;
			if (safeChains.count(belongsToChain.at(paths[path].path[i].id)) == 1)
			{
				currentPath.path.push_back(paths[path].path[i]);
				split = true;
			}
			if (split)
			{
				if (previousBelongedToComponent < components.size())
				{
					assert(currentPath.path.size() >= 2);
					for (size_t j = 1; j < currentPath.path.size()-1; j++)
					{
						assert(safeChains.count(belongsToChain.at(currentPath.path[j].id)) == 0);
					}
					result[previousBelongedToComponent].push_back(currentPath);
				}
				currentPath.index = path;
				currentPath.start = i;
				currentPath.path.clear();
				previousBelongedToComponent = components.size();
				if (edgeBelongsToComponent.count(currentEdge) == 1 && safeChains.count(belongsToChain.at(paths[path].path[i].id)) == 0)
				{
					currentPath.path.push_back(paths[path].path[i-1]);
					previousBelongedToComponent = edgeBelongsToComponent.at(currentEdge);
				}
			}
			currentPath.path.push_back(paths[path].path[i]);
		}
		if (previousBelongedToComponent < components.size())
		{
			assert(currentPath.path.size() >= 2);
			for (size_t j = 1; j < currentPath.path.size()-1; j++)
			{
				assert(safeChains.count(belongsToChain.at(currentPath.path[j].id)) == 0);
			}
			result[previousBelongedToComponent].push_back(currentPath);
		}
	}
	return result;
}

bool canResolve(const std::vector<Subpath>& pathsPerComponent, const ResolvableComponent& component, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	size_t totalSafeCrossing = 0;
	std::unordered_map<int, size_t> pathsCrossingPerSafe;
	std::unordered_map<int, size_t> safeCrossingPerSafe;
	for (auto path : pathsPerComponent)
	{
		if (safeChains.count(belongsToChain.at(path.path[0].id)) == 1) pathsCrossingPerSafe[path.path[0].id] += 1;
		if (safeChains.count(belongsToChain.at(path.path.back().id)) == 1) pathsCrossingPerSafe[path.path.back().id] += 1;
		if (safeChains.count(belongsToChain.at(path.path[0].id)) == 0 || safeChains.count(belongsToChain.at(path.path.back().id)) == 0) continue;
		totalSafeCrossing += 1;
		safeCrossingPerSafe[path.path[0].id] += 1;
		safeCrossingPerSafe[path.path.back().id] += 1;
	}
	std::unordered_set<int> nodeids = component.nodeIDs;
	for (auto edge : component.edges)
	{
		nodeids.insert(edge.first.id);
		nodeids.insert(edge.second.id);
	}
	for (auto node : nodeids)
	{
		if (safeChains.count(belongsToChain.at(node)) == 1)
		{
			std::cerr << "(" << node << ") ";
		}
		else
		{
			std::cerr << node << " ";
		}
	}
	if (pathsPerComponent.size() == 0)
	{
		std::cerr << "cannot resolve, zero paths" << std::endl;
		return false;
	}
	if (totalSafeCrossing < pathsPerComponent.size() * 0.5)
	{
		std::cerr << "cannot resolve, too few total safe crossers (" << totalSafeCrossing << ", " << pathsPerComponent.size() << ")" << std::endl;
		return false;
	}
	for (auto edge : component.edges)
	{
		if (safeChains.count(belongsToChain.at(edge.first.id)) == 1 && (safeCrossingPerSafe.count(edge.first.id) == 0 || safeCrossingPerSafe.at(edge.first.id) == 0))
		{
			std::cerr << "cannot resolve, zero safe crossers for node " << edge.first.id;
			return false;
		}
		if (safeChains.count(belongsToChain.at(edge.second.id)) == 1 && (safeCrossingPerSafe.count(edge.second.id) == 0 || safeCrossingPerSafe.at(edge.second.id) == 0))
		{
			std::cerr << "cannot resolve, zero safe crossers for node " << edge.second.id;
			return false;
		}
	}
	for (auto pair : pathsCrossingPerSafe)
	{
		if (safeCrossingPerSafe.count(pair.first) == 0 || safeCrossingPerSafe.at(pair.first) == 0)
		{
			std::cerr << "cannot resolve, zero safe crossers for node " << pair.first;
			return false;
		}
		if (safeCrossingPerSafe.at(pair.first) < pair.second * 0.5)
		{
			std::cerr << "cannot resolve, too few safe crossers for node " << pair.first << " (" << totalSafeCrossing << ", " << pathsPerComponent.size() << ")" << std::endl;
			return false;
		}
	}
	std::cerr << "resolve" << std::endl;
	return true;
}

bool canPartiallyResolve(const std::vector<Subpath>& pathsPerComponent, const ResolvableComponent& component, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	size_t totalSafeCrossing = 0;
	std::unordered_map<int, std::unordered_map<int, size_t>> crossers;
	for (auto path : pathsPerComponent)
	{
		if (safeChains.count(belongsToChain.at(path.path[0].id)) == 0 || safeChains.count(belongsToChain.at(path.path.back().id)) == 0) continue;
		crossers[path.path[0].id][path.path.back().id] += 1;
		crossers[path.path.back().id][path.path[0].id] += 1;
	}
	std::unordered_map<int, size_t> totalCrossers;
	for (auto pair : crossers)
	{
		for (auto pair2 : pair.second)
		{
			totalCrossers[pair.first] += pair2.second;
		}
	}
	std::unordered_set<std::pair<int, int>> potentiallyResolvable;
	for (auto pair : crossers)
	{
		size_t total = totalCrossers.at(pair.first);
		assert(total > 0);
		for (auto other : pair.second)
		{
			if (other.second > total * 0.7) potentiallyResolvable.emplace(pair.first, other.first);
		}
	}
	std::vector<std::pair<int, int>> resolvable;
	for (auto pair : potentiallyResolvable)
	{
		if (pair.second < pair.first) continue;
		if (potentiallyResolvable.count(std::make_pair(pair.second, pair.first)) == 1)
		{
			resolvable.push_back(pair);
		}
	}
	if (resolvable.size() > 0)
	{
		std::cerr << "can partially resolve";
		for (auto pair : resolvable)
		{
			std::cerr << " " << pair.first << "-" << pair.second;
		}
		std::cerr << std::endl;
		return true;
	}
	std::cerr << "cannot partially resolve" << std::endl;
	return false;
}

std::vector<std::pair<size_t, size_t>> align(const std::vector<std::pair<int, NodePos>>& nodes, const std::list<size_t>& order, const std::unordered_map<size_t, std::vector<size_t>>& inNeighbors, const std::vector<NodePos>& path, const std::unordered_map<int, size_t>& nodeSizes)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> backtrace;
	std::vector<std::vector<double>> score;
	score.resize(path.size());
	backtrace.resize(path.size());
	for (size_t i = 0; i < score.size(); i++)
	{
		score[i].resize(nodes.size(), 0);
		backtrace[i].resize(nodes.size());
	}
	for (size_t j : order)
	{
		for (size_t i = 0; i < path.size(); i++)
		{
			if (i == 0 && j == 0)
			{
				assert(inNeighbors.count(j) == 0);
				assert(nodes[0].second == path[0]);
				score[i][j] = nodeSizes.at(path[0].id);
				backtrace[i][j] = std::make_pair(0, 0);
			}
			else if (i == 0)
			{
				assert(inNeighbors.count(j) == 1);
				score[i][j] = std::numeric_limits<double>::min();
				size_t nodeSize = nodeSizes.at(nodes[j].second.id);
				for (auto neighbor : inNeighbors.at(j))
				{
					if (score[i][neighbor] - nodeSize > score[i][j])
					{
						score[i][j] = score[i][neighbor] - nodeSize;
						backtrace[i][j] = std::make_pair(i, neighbor);
					}
				}
			}
			else if (j == 0)
			{
				assert(inNeighbors.count(j) == 0);
				score[i][j] = score[i-1][j] - nodeSizes.at(path[i].id);
				backtrace[i][j] = std::make_pair(i-1, j);
			}
			else
			{
				assert(inNeighbors.count(j) == 1);
				bool match = path[i] == nodes[j].second;
				score[i][j] = score[i-1][j] - nodeSizes.at(path[i].id);
				backtrace[i][j] = std::make_pair(i-1, j);
				double deletionSize = nodeSizes.at(nodes[j].second.id);
				double matchSize = std::max(nodeSizes.at(nodes[j].second.id), nodeSizes.at(path[i].id));
				for (auto neighbor : inNeighbors.at(j))
				{
					if (match)
					{
						if (score[i-1][neighbor] + matchSize > score[i][j])
						{
							score[i][j] = score[i-1][neighbor] + matchSize;
							backtrace[i][j] = std::make_pair(i-1, neighbor);
						}
					}
					else
					{
						if (score[i-1][neighbor] - matchSize > score[i][j])
						{
							score[i][j] = score[i-1][neighbor] - matchSize;
							backtrace[i][j] = std::make_pair(i-1, neighbor);
						}
					}
					if (score[i][neighbor] - deletionSize > score[i][j])
					{
						score[i][j] = score[i][neighbor] - deletionSize;
						backtrace[i][j] = std::make_pair(i, neighbor);
					}
				}
			}
		}
	}
	std::vector<std::pair<size_t, size_t>> matches;
	size_t i = path.size()-1;
	size_t j = *(--order.end());
	while (i != 0 || j != 0)
	{
		bool match = path[i] == nodes[j].second;
		if (match && backtrace[i][j].first != i && backtrace[i][j].second != j) matches.emplace_back(i, j);
		std::tie(i, j) = backtrace[i][j];
	}
	matches.emplace_back(0, 0);
	std::reverse(matches.begin(), matches.end());
	return matches;
}

void resolve(int& nextNodeId, const std::unordered_map<int, size_t>& nodeSizes, std::vector<Subpath>& pathsPerComponent, ResolvableComponent& component, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	std::unordered_map<std::pair<NodePos, NodePos>, std::vector<size_t>> subpathsPerConnection;
	for (size_t i = 0; i < pathsPerComponent.size(); i++)
	{
		if (safeChains.count(belongsToChain.at(pathsPerComponent[i].path[0].id)) == 0) continue;
		if (safeChains.count(belongsToChain.at(pathsPerComponent[i].path.back().id)) == 0) continue;
		for (size_t j = 1; j < pathsPerComponent[i].path.size()-1; j++)
		{
			assert(safeChains.count(belongsToChain.at(pathsPerComponent[i].path[j].id)) == 0);
		}
		subpathsPerConnection[canon(pathsPerComponent[i].path[0], pathsPerComponent[i].path.back())].push_back(i);
	}
	for (auto pair : subpathsPerConnection)
	{
		assert(pair.second.size() > 0);
		for (auto i : pair.second)
		{
			if (pathsPerComponent[i].path[0] != pair.first.first)
			{
				std::reverse(pathsPerComponent[i].path.begin(), pathsPerComponent[i].path.end());
				for (size_t j = 0; j < pathsPerComponent[i].path.size(); j++)
				{
					pathsPerComponent[i].path[j] = pathsPerComponent[i].path[j].Reverse();
				}
			}
			assert(pathsPerComponent[i].path[0] == pair.first.first);
			assert(pathsPerComponent[i].path.back() == pair.first.second);
		}
		std::vector<std::pair<int, NodePos>> nodes;
		std::list<size_t> order;
		std::unordered_map<size_t, std::vector<size_t>> inNeighbors;
		assert(pathsPerComponent[pair.second[0]].path.size() >= 2);
		for (size_t i = 0; i < pathsPerComponent[pair.second[0]].path.size(); i++)
		{
			nodes.emplace_back(nextNodeId, pathsPerComponent[pair.second[0]].path[i]);
			order.push_back(i);
			if (i > 0)
			{
				inNeighbors[i].push_back(i-1);
			}
			nextNodeId++;
		}
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			auto matches = align(nodes, order, inNeighbors, pathsPerComponent[pair.second[i]].path, nodeSizes);
			assert(matches.size() >= 2);
			assert(matches[0].first == 0);
			assert(matches[0].second == 0);
			assert(matches.back().first == pathsPerComponent[pair.second[i]].path.size()-1);
			assert(matches.back().second == *(--order.end()));
			auto pos = order.begin();
			for (size_t j = 1; j < matches.size(); j++)
			{
				while (*pos != matches[j].second)
				{
					++pos;
					assert(pos != order.end());
				}
				size_t lastNode = matches[j-1].second;
				for (size_t k = matches[j-1].first+1; k < matches[j].first; k++)
				{
					order.insert(pos, nodes.size());
					nodes.emplace_back(nextNodeId, pathsPerComponent[pair.second[i]].path[k]);
					if (k == matches[j-1].first+1)
					{
						inNeighbors[nodes.size()-1].push_back(matches[j-1].second);
					}
					else
					{
						inNeighbors[nodes.size()-1].push_back(nodes.size()-2);
					}
					lastNode = nodes.size()-1;
					nextNodeId += 1;
				}
				inNeighbors[matches[j].second].push_back(lastNode);
			}
		}
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (i == 0 || i == *(--order.end()))
			{
				assert(safeChains.count(belongsToChain.at(nodes[i].second.id)) == 1);
				continue;
			}
			assert(safeChains.count(belongsToChain.at(nodes[i].second.id)) == 0);
			assert(component.newNodes.count(nodes[i].first) == 0);
			component.newNodes[nodes[i].first] = nodes[i].second;
		}
		for (auto pair : inNeighbors)
		{
			size_t to = pair.first;
			for (auto from : pair.second)
			{
				NodePos fromPos;
				NodePos toPos;
				fromPos.id = nodes[from].first;
				toPos.id = nodes[to].first;
				fromPos.end = true;
				toPos.end = true;
				if (from == 0)
				{
					fromPos = nodes[0].second;
				}
				if (to == *(--order.end()))
				{
					toPos = nodes[to].second;
				}
				component.newEdges.emplace_back(fromPos, toPos);
			}
		}
	}
	component.removedNodes = component.nodeIDs;
}

void resolvePartially(int& nextNodeId, const std::unordered_map<int, size_t>& nodeSizes, std::vector<Subpath>& pathsPerComponent, ResolvableComponent& component, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain, const GfaGraph& graph)
{
	size_t totalSafeCrossing = 0;
	std::unordered_map<int, std::unordered_map<int, size_t>> crossers;
	for (auto path : pathsPerComponent)
	{
		if (safeChains.count(belongsToChain.at(path.path[0].id)) == 0 || safeChains.count(belongsToChain.at(path.path.back().id)) == 0) continue;
		crossers[path.path[0].id][path.path.back().id] += 1;
		crossers[path.path.back().id][path.path[0].id] += 1;
	}
	std::unordered_map<int, size_t> totalCrossers;
	for (auto pair : crossers)
	{
		for (auto pair2 : pair.second)
		{
			totalCrossers[pair.first] += pair2.second;
		}
	}
	std::unordered_set<std::pair<int, int>> potentiallyResolvable;
	for (auto pair : crossers)
	{
		size_t total = totalCrossers.at(pair.first);
		assert(total > 0);
		for (auto other : pair.second)
		{
			if (other.second > total * 0.7) potentiallyResolvable.emplace(pair.first, other.first);
		}
	}
	std::unordered_set<std::pair<int, int>> resolvable;
	for (auto pair : potentiallyResolvable)
	{
		if (pair.second < pair.first) continue;
		if (potentiallyResolvable.count(std::make_pair(pair.second, pair.first)) == 1)
		{
			resolvable.insert(pair);
		}
	}
	assert(resolvable.size() > 0);
	std::unordered_map<std::pair<NodePos, NodePos>, std::vector<size_t>> subpathsPerConnection;
	for (size_t i = 0; i < pathsPerComponent.size(); i++)
	{
		if (safeChains.count(belongsToChain.at(pathsPerComponent[i].path[0].id)) == 0) continue;
		if (safeChains.count(belongsToChain.at(pathsPerComponent[i].path.back().id)) == 0) continue;
		auto key = canon(pathsPerComponent[i].path[0], pathsPerComponent[i].path.back());
		std::pair<int, int> check { key.first.id, key.second.id };
		if (resolvable.count(check) == 0) continue;
		subpathsPerConnection[key].push_back(i);
	}
	assert(subpathsPerConnection.size() == resolvable.size());
	for (auto pair : subpathsPerConnection)
	{
		assert(graph.edges.count(pair.first.first) == 1);
		assert(graph.edges.count(pair.first.second.Reverse()) == 1);
		for (auto edge : graph.edges.at(pair.first.first))
		{
			component.removedEdges.emplace(pair.first.first, edge);
		}
		for (auto edge : graph.edges.at(pair.first.second.Reverse()))
		{
			component.removedEdges.emplace(pair.first.second.Reverse(), edge);
		}
		assert(pair.second.size() > 0);
		for (auto i : pair.second)
		{
			if (pathsPerComponent[i].path[0] != pair.first.first)
			{
				std::reverse(pathsPerComponent[i].path.begin(), pathsPerComponent[i].path.end());
				for (size_t j = 0; j < pathsPerComponent[i].path.size(); j++)
				{
					pathsPerComponent[i].path[j] = pathsPerComponent[i].path[j].Reverse();
				}
			}
			assert(pathsPerComponent[i].path[0] == pair.first.first);
			assert(pathsPerComponent[i].path.back() == pair.first.second);
		}
		std::vector<std::pair<int, NodePos>> nodes;
		std::list<size_t> order;
		std::unordered_map<size_t, std::vector<size_t>> inNeighbors;
		assert(pathsPerComponent[pair.second[0]].path.size() >= 2);
		for (size_t i = 0; i < pathsPerComponent[pair.second[0]].path.size(); i++)
		{
			nodes.emplace_back(nextNodeId, pathsPerComponent[pair.second[0]].path[i]);
			order.push_back(i);
			if (i > 0)
			{
				inNeighbors[i].push_back(i-1);
			}
			nextNodeId++;
		}
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			auto matches = align(nodes, order, inNeighbors, pathsPerComponent[pair.second[i]].path, nodeSizes);
			assert(matches.size() >= 2);
			assert(matches[0].first == 0);
			assert(matches[0].second == 0);
			assert(matches.back().first == pathsPerComponent[pair.second[i]].path.size()-1);
			assert(matches.back().second == *(--order.end()));
			auto pos = order.begin();
			for (size_t j = 1; j < matches.size(); j++)
			{
				while (*pos != matches[j].second)
				{
					++pos;
					assert(pos != order.end());
				}
				size_t lastNode = matches[j-1].second;
				for (size_t k = matches[j-1].first+1; k < matches[j].first; k++)
				{
					order.insert(pos, nodes.size());
					nodes.emplace_back(nextNodeId, pathsPerComponent[pair.second[i]].path[k]);
					if (k == matches[j-1].first+1)
					{
						inNeighbors[nodes.size()-1].push_back(matches[j-1].second);
					}
					else
					{
						inNeighbors[nodes.size()-1].push_back(nodes.size()-2);
					}
					lastNode = nodes.size()-1;
					nextNodeId += 1;
				}
				inNeighbors[matches[j].second].push_back(lastNode);
			}
		}
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (i == 0 || i == *(--order.end()))
			{
				assert(safeChains.count(belongsToChain.at(nodes[i].second.id)) == 1);
				continue;
			}
			assert(safeChains.count(belongsToChain.at(nodes[i].second.id)) == 0);
			assert(component.newNodes.count(nodes[i].first) == 0);
			component.newNodes[nodes[i].first] = nodes[i].second;
		}
		for (auto pair : inNeighbors)
		{
			size_t to = pair.first;
			for (auto from : pair.second)
			{
				NodePos fromPos;
				NodePos toPos;
				fromPos.id = nodes[from].first;
				toPos.id = nodes[to].first;
				fromPos.end = true;
				toPos.end = true;
				if (from == 0)
				{
					fromPos = nodes[0].second;
				}
				if (to == *(--order.end()))
				{
					toPos = nodes[to].second;
				}
				component.newEdges.emplace_back(fromPos, toPos);
			}
		}
	}
}

void updateGraph(GfaGraph& graph, const ResolvableComponent& component, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain)
{
	for (auto node : component.newNodes)
	{
		assert(safeChains.count(belongsToChain.at(node.second.id)) == 0);
		assert(graph.nodes.count(node.first) == 0);
		assert(graph.nodes.count(node.second.id) == 1);
		std::string seq = graph.nodes.at(node.second.id);
		if (!node.second.end) seq = CommonUtils::ReverseComplement(seq);
		graph.nodes[node.first] = seq;
	}
	for (auto edge : component.newEdges)
	{
		assert(graph.nodes.count(edge.first.id) == 1);
		assert(graph.nodes.count(edge.second.id) == 1);
		graph.edges[edge.first].push_back(edge.second);
		NodePos oldFrom = edge.first;
		if (belongsToChain.count(oldFrom.id) == 0 || safeChains.count(belongsToChain.at(oldFrom.id)) == 0) oldFrom = component.newNodes.at(oldFrom.id);
		if (!edge.first.end) oldFrom = oldFrom.Reverse();
		NodePos oldTo = edge.second;
		if (belongsToChain.count(oldTo.id) == 0 || safeChains.count(belongsToChain.at(oldTo.id)) == 0) oldTo = component.newNodes.at(oldTo.id);
		if (!edge.second.end) oldTo = oldTo.Reverse();
		if (graph.varyingOverlaps.count(std::make_pair(oldFrom, oldTo)) == 1) graph.varyingOverlaps[edge] = graph.varyingOverlaps.at(std::make_pair(oldFrom, oldTo));
	}
	for (auto node : component.removedNodes)
	{
		assert(safeChains.count(belongsToChain.at(node)) == 0);
		graph.nodes.erase(node);
	}
	for (auto edge : component.removedEdges)
	{
		if (graph.edges.count(edge.first) == 1)
		{
			for (size_t i = 0; i < graph.edges.at(edge.first).size(); i++)
			{
				if (graph.edges.at(edge.first)[i] == edge.second)
				{
					std::swap(graph.edges.at(edge.first)[i], graph.edges.at(edge.first).back());
					graph.edges.at(edge.first).pop_back();
					break;
				}
			}
		}
		if (graph.edges.count(edge.second.Reverse()) == 1)
		{
			for (size_t i = 0; i < graph.edges.at(edge.second.Reverse()).size(); i++)
			{
				if (graph.edges.at(edge.second.Reverse())[i] == edge.first.Reverse())
				{
					std::swap(graph.edges.at(edge.second.Reverse())[i], graph.edges.at(edge.second.Reverse()).back());
					graph.edges.at(edge.second.Reverse()).pop_back();
					break;
				}
			}
		}
	}
}

void resolveComponentsAndReplacePaths(GfaGraph& graph, const std::unordered_set<int>& safeChains, const std::unordered_map<int, int>& belongsToChain, std::vector<ResolvableComponent>& components, std::vector<Path>& paths)
{
	auto pathsPerComponent = splitPathsPerComponent(paths, components, safeChains, belongsToChain);
	std::unordered_map<int, size_t> nodeSizes;
	int nextNodeId = 0;
	for (auto node : graph.nodes)
	{
		nextNodeId = std::max(nextNodeId, node.first);
		if (node.second.size() > graph.edgeOverlap)
		{
			nodeSizes[node.first] = node.second.size() - graph.edgeOverlap;
		}
		else
		{
			nodeSizes[node.first] = 1;
		}
	}
	nextNodeId += 1;
	size_t unresolvableComponents = 0;
	size_t resolvedComponents = 0;
	size_t partiallyResolvedComponents = 0;
	size_t tooBigs = 0;
	for (size_t i = 0; i < components.size(); i++)
	{
		if (components[i].edges.size() == 0)
		{
			continue;
		}
		size_t size = 0;
		for (auto node : components[i].nodeIDs)
		{
			size_t nodeSize = graph.nodes.at(node).size();
			if (nodeSize > graph.edgeOverlap) nodeSize -= graph.edgeOverlap;
			size += nodeSize;
		}
		// if (canResolve(pathsPerComponent[i], components[i], safeChains, belongsToChain))
		// {
		// 	resolve(nextNodeId, nodeSizes, pathsPerComponent[i], components[i], safeChains, belongsToChain);
		// 	resolvedComponents += 1;
		// 	updateGraph(graph, components[i], safeChains, belongsToChain);
		// }
		if (canPartiallyResolve(pathsPerComponent[i], components[i], safeChains, belongsToChain))
		{
			resolvePartially(nextNodeId, nodeSizes, pathsPerComponent[i], components[i], safeChains, belongsToChain, graph);
			partiallyResolvedComponents += 1;
			updateGraph(graph, components[i], safeChains, belongsToChain);
		}
		else
		{
			unresolvableComponents += 1;
		}
	}
	std::cerr << resolvedComponents << " components resolved, " << partiallyResolvedComponents << " partially resolved, " << unresolvableComponents << " unresolved, " << tooBigs << " too big components" << std::endl;
}

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string alignmentFile { argv[2] };
	size_t safeChainSize = std::stol(argv[3]);
	size_t genomeSize = std::stol(argv[4]);
	std::string outputGraphFile { argv[5] };
	std::string translatedAlignmentFile { argv[6] };

	auto graph = GfaGraph::LoadFromFile(graphFile, true);
	graph.confirmDoublesidedEdges();
	auto alninfo = loadAlignmentsAsPaths(alignmentFile);
	auto paths = std::get<0>(alninfo);
	auto coveragePerNode = std::get<1>(alninfo);
	size_t bpCount = std::get<2>(alninfo);
	auto safeChains = getSafeChains(graph, safeChainSize, coveragePerNode, (double)bpCount/(double)genomeSize);
	std::cerr << safeChains.size() << " safe chains" << std::endl;
	auto belongsToChain = getChainBelongers(graph);
	auto components = getComponents(graph, safeChains, belongsToChain);
	std::cerr << components.size() << " components" << std::endl;
	resolveComponentsAndReplacePaths(graph, safeChains, belongsToChain, components, paths);
	graph.SaveToFile(outputGraphFile);
	// writePaths(translatedAlignmentFile);
}
