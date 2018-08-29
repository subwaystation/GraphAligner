#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "CommonUtils.h"

std::pair<std::vector<std::vector<NodePos>>, std::unordered_map<NodePos, NodePos>> getUnitigs(const GfaGraph& graph)
{
	std::vector<std::vector<NodePos>> result;
	std::unordered_map<NodePos, NodePos> nodeToUnitigMapping;
	for (auto node : graph.nodes)
	{
		NodePos fw, bw;
		fw.id = node.first;
		fw.end = true;
		bw = fw.Reverse();
		if (nodeToUnitigMapping.count(fw) == 1)
		{
			assert(nodeToUnitigMapping.count(bw) == 1);
			continue;
		}
		result.emplace_back();
		nodeToUnitigMapping[fw] = NodePos { result.size()-1, true };
		nodeToUnitigMapping[fw.Reverse()] = NodePos { result.size()-1, false };
		while (graph.edges.count(bw) == 1 && graph.edges.at(bw).size() == 1)
		{
			NodePos next = graph.edges.at(bw)[0];
			assert(graph.edges.count(next.Reverse()) == 1);
			if (graph.edges.at(next.Reverse()).size() != 1) break;
			if (next.id == node.first) break;
			bw = next;
			result.back().push_back(bw.Reverse());
			assert(nodeToUnitigMapping.count(bw) == 0);
			assert(nodeToUnitigMapping.count(bw.Reverse()) == 0);
			nodeToUnitigMapping[bw] = NodePos { result.size()-1, false };
			nodeToUnitigMapping[bw.Reverse()] = NodePos { result.size()-1, true };
		}
		std::reverse(result.back().begin(), result.back().end());
		result.back().push_back(fw);
		while (graph.edges.count(fw) == 1 && graph.edges.at(fw).size() == 1)
		{
			NodePos next = graph.edges.at(fw)[0];
			assert(graph.edges.count(next.Reverse()) == 1);
			if (graph.edges.at(next.Reverse()).size() != 1) break;
			if (next.id == node.first) break;
			fw = next;
			result.back().push_back(fw);
			assert(nodeToUnitigMapping.count(fw) == 0);
			assert(nodeToUnitigMapping.count(fw.Reverse()) == 0);
			nodeToUnitigMapping[fw] = NodePos { result.size()-1, true };
			nodeToUnitigMapping[fw.Reverse()] = NodePos { result.size()-1, false };
		}
	}
	return std::make_pair(result, nodeToUnitigMapping);
}

std::vector<std::vector<NodePos>> getAlnUnitigPaths(const std::vector<vg::Alignment>& alns, const std::unordered_map<NodePos, NodePos>& nodeToUnitigMapping)
{
	std::vector<std::vector<NodePos>> result;
	for (auto aln : alns)
	{
		result.emplace_back();
		NodePos lastUnitig = nodeToUnitigMapping.at(NodePos { aln.path().mapping(0).position().node_id(), !aln.path().mapping(0).position().is_reverse() });
		for (int i = 1; i < aln.path().mapping_size(); i++)
		{
			NodePos thisUnitig = nodeToUnitigMapping.at(NodePos { aln.path().mapping(i).position().node_id(), !aln.path().mapping(i).position().is_reverse() });
			if (thisUnitig != lastUnitig)
			{
				result.back().push_back(lastUnitig);
			}
			lastUnitig = thisUnitig;
		}
		result.back().push_back(lastUnitig);
	}
	return result;
}

std::unordered_set<int> getRepeatUnitigs(const GfaGraph& graph, size_t safeNonrepeatLength, const std::vector<std::vector<NodePos>>& unitigs, const std::unordered_map<NodePos, NodePos>& nodeToUnitigMapping)
{
	std::unordered_set<int> result;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		size_t unitigLength = 0;
		for (auto node : unitigs[i])
		{
			unitigLength += graph.nodes.at(node.id).size() - graph.edgeOverlap;
		}
		if (unitigLength >= safeNonrepeatLength) continue;
		NodePos leftEnd = unitigs[i][0].Reverse();
		if (graph.edges.count(leftEnd) == 1 && graph.edges.at(leftEnd).size() > 1)
		{
			result.insert(i);
			continue;
		}
		NodePos rightEnd = unitigs[i].back();
		if (graph.edges.count(rightEnd) == 1 && graph.edges.at(rightEnd).size() > 1)
		{
			result.insert(i);
			continue;
		}
	}
	return result;
}

std::vector<std::vector<NodePos>> getRepeatCrossingPaths(const std::unordered_set<int>& repeatUnitigs, const std::vector<std::vector<NodePos>>& paths)
{
	std::vector<std::vector<NodePos>> result;
	for (auto path : paths)
	{
		size_t lastNonrepetitive = 0;
		bool currentlyRepeat = repeatUnitigs.count(path[0].id) == 1;
		for (size_t i = 1; i < path.size(); i++)
		{
			if (repeatUnitigs.count(path[i].id) == 0)
			{
				if (currentlyRepeat)
				{
					result.emplace_back(path.begin() + lastNonrepetitive, path.begin() + i + 1);
					currentlyRepeat = false;
				}
				lastNonrepetitive = i;
				currentlyRepeat = false;
			}
			else
			{
				currentlyRepeat = true;
			}
		}
	}
	return result;
}

std::vector<std::vector<NodePos>> getRepeatPassthroughPaths(const std::unordered_set<int>& repeatUnitigs, const std::vector<std::vector<NodePos>>& paths)
{
	std::vector<std::vector<NodePos>> result;
	for (auto path : paths)
	{
		assert(path.size() >= 1);
		assert(path.size() >= 2 || repeatUnitigs.count(path[0].id) == 1);
		assert(path.size() >= 3 || (repeatUnitigs.count(path[0].id) == 1 || repeatUnitigs.count(path[1].id) == 1));
#ifndef NDEBUG
		for (size_t i = 1; i < path.size()-1; i++)
		{
			assert(repeatUnitigs.count(path[i].id) == 1);
		}
#endif
		if (repeatUnitigs.count(path[0].id) == 0 && repeatUnitigs.count(path.back().id) == 0) result.emplace_back(path);
	}
	return result;
}

std::vector<std::vector<NodePos>> canonizePaths(const std::vector<std::vector<NodePos>>& paths)
{
	std::vector<std::vector<NodePos>> result;
	for (auto path : paths)
	{
		assert(path.size() >= 3);
		result.push_back(path);
		if (path.back() < path[0])
		{
			std::reverse(result.back().begin(), result.back().end());
			for (auto& node : result.back())
			{
				node = node.Reverse();
			}
		}
	}
	return result;
}

std::vector<std::vector<NodePos>> getUniquePaths(const std::vector<std::vector<NodePos>>& paths)
{
	std::set<std::vector<NodePos>> uniques { paths.begin(), paths.end() };
	std::vector<std::vector<NodePos>> result { uniques.begin(), uniques.end() };
	return result;
}

bool pathsAreCompatible(const std::vector<NodePos>& subpath, const std::vector<NodePos>& superpath)
{
	if (subpath.size() > superpath.size()) return false;
	for (size_t start = 0; start < superpath.size() - subpath.size() + 1; start++)
	{
		bool compatible = true;
		for (size_t i = 0; i < subpath.size(); i++)
		{
			if (subpath[i] != superpath[start+i])
			{
				compatible = false;
				break;
			}
		}
		if (compatible) return true;
	}
	for (size_t start = 0; start < superpath.size() - subpath.size() + 1; start++)
	{
		bool compatible = true;
		for (size_t i = 0; i < subpath.size(); i++)
		{
			if (subpath[subpath.size() - 1 - i].Reverse() != superpath[start+i])
			{
				compatible = false;
				break;
			}
		}
		if (compatible) return true;
	}
	return false;
}

std::vector<std::vector<NodePos>> resolvePaths(const std::vector<std::vector<NodePos>>& throughPaths, const std::vector<std::vector<NodePos>>& paths)
{
	std::unordered_map<int, std::vector<size_t>> nodeCrossingThroughPaths;
	for (size_t i = 0; i < throughPaths.size(); i++)
	{
		for (auto node : throughPaths[i])
		{
			nodeCrossingThroughPaths[node.id].push_back(i);
		}
	}
	for (auto path : paths)
	{
		std::set<size_t> possibleCompatibles;
		for (auto node : path)
		{
			possibleCompatibles.insert(nodeCrossingThroughPaths[node.id].begin(), nodeCrossingThroughPaths[node.id].end());
		}
		bool compatible = false;
		for (auto possibleCompatible : possibleCompatibles)
		{
			if (pathsAreCompatible(path, throughPaths[possibleCompatible]))
			{
				compatible = true;
				break;
			}
		}
		if (compatible) continue;
		// ??? not compatible, what do ???
		assert(false);
	}
	return throughPaths;
}

std::vector<NodePos> getUnitigPathNodes(const std::vector<NodePos>& unitigPath, const std::vector<std::vector<NodePos>>& unitigs)
{
	assert(unitigPath.size() >= 3);
	std::vector<NodePos> nodes;
	if (!unitigPath[0].end)
	{
		nodes.push_back(unitigs[unitigPath[0].id].back());
	}
	else
	{
		nodes.push_back(unitigs[unitigPath[0].id][0].Reverse());
	}
	for (size_t i = 1; i < unitigPath.size()-1; i++)
	{
		auto u = unitigPath[i];
		if (!u.end)
		{
			for (auto n : unitigs[u.id])
			{
				nodes.push_back(n);
			}
		}
		else
		{
			for (size_t i = unitigs[u.id].size()-1; i < unitigs[u.id].size(); i--)
			{
				nodes.push_back(unitigs[u.id][i].Reverse());
			}
		}
	}
	if (!unitigPath.back().end)
	{
		nodes.push_back(unitigs[unitigPath.back().id][0]);
	}
	else
	{
		nodes.push_back(unitigs[unitigPath.back().id].back().Reverse());
	}
	return nodes;
}

void writeGraph(const GfaGraph& graph, const std::vector<std::vector<NodePos>>& unitigResolvedPaths, const std::vector<std::vector<NodePos>>& unitigs, std::string filename)
{
	std::ofstream file { filename };
	std::unordered_set<int> resolvedNodes;
	for (auto p : unitigResolvedPaths)
	{
		std::vector<NodePos> nodes = getUnitigPathNodes(p, unitigs);
		assert(nodes.size() >= 3);
		for (size_t i = 1; i < nodes.size()-1; i++)
		{
			resolvedNodes.insert(nodes[i].id);
		}
	}
	int newId = 0;
	for (auto node : graph.nodes)
	{
		if (resolvedNodes.count(node.first) == 1) continue;
		file << "S\t" << node.first << "\t" << node.second << std::endl;
		newId = std::max(newId, node.first);
	}
	for (auto edge : graph.edges)
	{
		NodePos source = edge.first;
		if (resolvedNodes.count(source.id) == 1) continue;
		for (auto target : edge.second)
		{
			if (resolvedNodes.count(target.id) == 1) continue;
			file << "L\t" << source.id << "\t" << (source.end ? "+" : "-") << "\t" << target.id << "\t" << (target.end ? "+" : "-") << "\t" << graph.edgeOverlap << "M" << std::endl;
		}
	}
	newId++;
	for (auto p : unitigResolvedPaths)
	{
		std::vector<NodePos> nodes = getUnitigPathNodes(p, unitigs);
		assert(nodes.size() >= 3);
		assert(resolvedNodes.count(nodes[0].id) == 0);
		assert(resolvedNodes.count(nodes.back().id) == 0);
		for (size_t i = 1; i < nodes.size() - 1; i++)
		{
			assert(resolvedNodes.count(nodes[i].id) == 1);
			std::string seq = graph.nodes.at(nodes[i].id);
			if (!nodes[i].end) seq = CommonUtils::ReverseComplement(seq);
			file << "S\t" << newId + i - 1 << "\t" << seq << std::endl;
		}
		file << "L\t" << nodes[0].id << "\t" << (nodes[0].end ? "+" : "-") << "\t" << newId << "\t+\t" << graph.edgeOverlap << "M" << std::endl; 
		for (size_t i = 1; i < nodes.size() - 2; i++)
		{
			file << "L\t" << (newId + i - 1) << "\t+\t" << (newId + i) << "\t+\t" << graph.edgeOverlap << "M" << std::endl; 
		}
		file << "L\t" << (newId + nodes.size() - 3) << "\t+\t" << nodes.back().id << "\t" << (nodes.back().end ? "+" : "-") << "\t" << graph.edgeOverlap << "M" << std::endl;
		newId += nodes.size() - 2;
	}
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	int safeNonrepeatLength = std::stoi(argv[3]);
	std::string outputGraph { argv[4] };

	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	auto alns = CommonUtils::LoadVGAlignments(inputAlns);
	auto unitigs = getUnitigs(graph);
	auto alnUnitigPaths = getAlnUnitigPaths(alns, unitigs.second);
	auto repeatUnitigs = getRepeatUnitigs(graph, safeNonrepeatLength, unitigs.first, unitigs.second);
	auto paths = getRepeatCrossingPaths(repeatUnitigs, alnUnitigPaths);
	auto uniquePaths = getUniquePaths(paths);
	auto throughPaths = getRepeatPassthroughPaths(repeatUnitigs, uniquePaths);
	auto canonThroughs = canonizePaths(throughPaths);
	auto uniqueThroughs = getUniquePaths(canonThroughs);
	auto resolved = resolvePaths(uniqueThroughs, uniquePaths);
	writeGraph(graph, resolved, unitigs.first, outputGraph);
}