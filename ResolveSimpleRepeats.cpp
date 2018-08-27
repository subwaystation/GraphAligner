#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "CommonUtils.h"

std::unordered_set<int> getSolvableNodes(const GfaGraph& graph)
{
	std::unordered_set<int> result;
	for (auto node : graph.nodes)
	{
		NodePos pos;
		pos.id = node.first;
		pos.end = true;
		if (graph.edges.count(pos) == 0) continue;
		if (graph.edges.at(pos).size() <= 1) continue;
		pos.end = false;
		if (graph.edges.count(pos) == 0) continue;
		if (graph.edges.at(pos).size() <= 1) continue;
		result.insert(node.first);
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, NodePos>> getTriplets(const std::unordered_set<int>& solvableNodes, const std::vector<vg::Alignment>& alns)
{
	std::vector<std::tuple<NodePos, NodePos, NodePos>> result;
	for (auto aln : alns)
	{
		for (int i = 1; i < aln.path().mapping_size()-1; i++)
		{
			if (solvableNodes.count(aln.path().mapping(i).position().node_id()) == 1)
			{
				result.emplace_back();
				std::get<0>(result.back()).id = aln.path().mapping(i-1).position().node_id();
				std::get<0>(result.back()).end = !aln.path().mapping(i-1).position().is_reverse();
				std::get<1>(result.back()).id = aln.path().mapping(i).position().node_id();
				std::get<1>(result.back()).end = !aln.path().mapping(i).position().is_reverse();
				std::get<2>(result.back()).id = aln.path().mapping(i+1).position().node_id();
				std::get<2>(result.back()).end = !aln.path().mapping(i+1).position().is_reverse();
			}
		}
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, NodePos>> getUniqueTriplets(const std::vector<std::tuple<NodePos, NodePos, NodePos>>& triplets)
{
	std::set<std::tuple<NodePos, NodePos, NodePos>> uniques { triplets.begin(), triplets.end() };
	std::vector<std::tuple<NodePos, NodePos, NodePos>> result { uniques.begin(), uniques.end() };
	return result;
}

std::vector<std::tuple<NodePos, std::unordered_set<NodePos>, std::unordered_set<NodePos>>> resolveTriplets(const std::vector<std::tuple<NodePos, NodePos, NodePos>>& triplets)
{
	std::vector<std::tuple<NodePos, std::unordered_set<NodePos>, std::unordered_set<NodePos>>> result;
	std::unordered_map<NodePos, std::vector<std::pair<NodePos, NodePos>>> tripletsPerNode;
	for (auto tr : triplets)
	{
		tripletsPerNode[std::get<1>(tr)].emplace_back(std::get<0>(tr), std::get<2>(tr));
	}
	for (auto node : tripletsPerNode)
	{
		std::unordered_map<NodePos, size_t> component;
		for (auto pair : node.second)
		{
			size_t componentNum = component.size();
			component[pair.first] = componentNum;
			componentNum = component.size();
			component[pair.second] = componentNum;
		}
		bool repeat = true;
		while (repeat)
		{
			repeat = false;
			for (auto pair : node.second)
			{
				if (component[pair.first] != component[pair.second])
				{
					size_t replaceThis = component[pair.second];
					size_t replaceWith = component[pair.first]; 
					for (auto& p : component)
					{
						if (p.second == replaceThis) p.second = replaceWith;
					}
					repeat = true;
				}
			}
		}
		size_t numComponents = 0;
		std::unordered_map<size_t, size_t> componentMap;
		for (auto pair : component)
		{
			if (componentMap.count(pair.second) == 0)
			{
				componentMap[pair.second] = numComponents;
				numComponents++;
			}
		}
		assert(numComponents != -1);
		std::vector<std::tuple<NodePos, std::unordered_set<NodePos>, std::unordered_set<NodePos>>> partialResult;
		partialResult.resize(numComponents);
		for (auto pair : componentMap)
		{
			assert(std::get<0>(partialResult[pair.second]).id == 0 || std::get<0>(partialResult[pair.second]) == node.first);
			std::get<0>(partialResult[pair.second]) = node.first;
			for (auto c : node.second)
			{
				if (componentMap[component[c.first]] == pair.second)
				{
					std::get<1>(partialResult[pair.second]).insert(c.first);
				}
				if (componentMap[component[c.second]] == pair.second)
				{
					std::get<2>(partialResult[pair.second]).insert(c.second);
				}
			}
		}
		for (size_t i = 0; i < numComponents; i++)
		{
			assert(std::get<0>(partialResult[i]).id != 0);
		}
		result.insert(result.end(), partialResult.begin(), partialResult.end());
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, NodePos>> canonizeTriplets(const std::vector<std::tuple<NodePos, NodePos, NodePos>>& triplets)
{
	std::vector<std::tuple<NodePos, NodePos, NodePos>> result;
	for (auto t : triplets)
	{
		if (std::get<1>(t).end)
		{
			result.emplace_back(std::get<2>(t).Reverse(), std::get<1>(t).Reverse(), std::get<0>(t).Reverse());
		}
		else
		{
			result.push_back(t);
		}
	}
	return result;
}

void writeGraph(const GfaGraph& graph, const std::vector<std::tuple<NodePos, std::unordered_set<NodePos>, std::unordered_set<NodePos>>>& resolved, std::string filename)
{
	std::ofstream file { filename };
	std::unordered_set<int> resolvedNodes;
	for (auto r : resolved)
	{
		resolvedNodes.insert(std::get<0>(r).id);
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
	for (auto r : resolved)
	{
		assert(graph.nodes.count(std::get<0>(r).id) == 1);
		file << "S\t" << newId << "\t" << graph.nodes.at(std::get<0>(r).id) << std::endl;
		for (auto inNeighbor : std::get<1>(r))
		{
			file << "L\t" << inNeighbor.id << "\t" << (inNeighbor.end ? "-" : "+") << "\t" << newId << "\t+\t" << graph.edgeOverlap << "M" << std::endl;
		}
		for (auto outNeighbor : std::get<2>(r))
		{
			file << "L\t" << newId << "\t+\t" << outNeighbor.id << "\t" << (outNeighbor.end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
		}
		newId++;
	}
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	std::string outputGraph { argv[3] };

	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	auto alns = CommonUtils::LoadVGAlignments(inputAlns);
	auto solvableNodes = getSolvableNodes(graph);
	auto triplets = getTriplets(solvableNodes, alns);
	auto canons = canonizeTriplets(triplets);
	auto uniques = getUniqueTriplets(canons);
	auto resolved = resolveTriplets(uniques);
	writeGraph(graph, resolved, outputGraph);
}