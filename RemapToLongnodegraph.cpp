#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include "GfaGraph.h"
#include "CommonUtils.h"
#include "stream.hpp"

std::unordered_set<int> getExistingNodes(const GfaGraph& graph, const std::vector<vg::Alignment>& alns)
{
	std::unordered_set<int> existsInAlns;
	for (auto aln : alns)
	{
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			existsInAlns.insert(aln.path().mapping(i).position().node_id());
		}
	}
	std::unordered_set<int> result;
	for (auto node : graph.nodes)
	{
		if (existsInAlns.count(node.first) == 1)
		{
			result.insert(node.first);
		}
	}
	std::cout << result.size() << " solid nodes" << std::endl;
	return result;
}

std::unordered_map<std::pair<NodePos, NodePos>, std::pair<NodePos, bool>> getConnectors(const GfaGraph& graph, const std::unordered_set<int>& existingNodes)
{
	std::unordered_map<std::pair<NodePos, NodePos>, std::pair<NodePos, bool>> result;
	for (auto node : graph.nodes)
	{
		if (existingNodes.count(node.first) == 0) continue;
		NodePos start;
		start.id = node.first;
		start.end = false;
		if (graph.edges.count(start) == 1)
		{
			for (auto target : graph.edges.at(start))
			{
				if (existingNodes.count(target.id) == 1)
				{
					result[std::make_pair(start, target)] = std::make_pair(target, true);
					result[std::make_pair(target.Reverse(), start.Reverse())] = std::make_pair(target.Reverse(), true);
					continue;
				}
				assert(graph.edges.count(target) == 1);
				for (auto target2 : graph.edges.at(target))
				{
					assert(existingNodes.count(target2.id) == 1);
					result[std::make_pair(start, target2)] = std::make_pair(target, false);
					result[std::make_pair(target2.Reverse(), start.Reverse())] = std::make_pair(target.Reverse(), false);
				}
			}
		}
		start.id = node.first;
		start.end = true;
		if (graph.edges.count(start) == 1)
		{
			for (auto target : graph.edges.at(start))
			{
				if (existingNodes.count(target.id) == 1)
				{
					result[std::make_pair(start, target)] = std::make_pair(target, true);
					result[std::make_pair(target.Reverse(), start.Reverse())] = std::make_pair(target.Reverse(), true);
					continue;
				}
				assert(graph.edges.count(target) == 1);
				for (auto target2 : graph.edges.at(target))
				{
					assert(existingNodes.count(target2.id) == 1);
					result[std::make_pair(start, target2)] = std::make_pair(target, false);
					result[std::make_pair(target2.Reverse(), start.Reverse())] = std::make_pair(target.Reverse(), false);
				}
			}
		}
	}
	return result;
}

std::vector<vg::Alignment> remapAlignments(const std::vector<vg::Alignment>& oldAlns, const std::unordered_set<int>& existingNodes, const std::unordered_map<std::pair<NodePos, NodePos>, std::pair<NodePos, bool>>& connectors)
{
	std::vector<vg::Alignment> result;
	for (auto aln : oldAlns)
	{
		std::vector<NodePos> solids;
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			if (existingNodes.count(aln.path().mapping(i).position().node_id()) == 1)
			{
				solids.emplace_back(aln.path().mapping(i).position().node_id(), !aln.path().mapping(i).position().is_reverse());
			}
		}
		// for (size_t i = 0; i < solids.size(); i++)
		// {
		// 	std::cout << solids[i].id << (solids[i].end ? "-" : "+") << " ";
		// }
		// std::cout << std::endl;
		vg::Alignment oneResult;
		for (size_t i = 1; i < solids.size(); i++)
		{
			auto pair = std::make_pair(solids[i-1], solids[i]);
			if (connectors.count(pair) == 0)
			{
				assert(connectors.count(std::make_pair(solids[i].Reverse(), solids[i-1].Reverse())) == 0);
				if (oneResult.path().mapping_size() > 0) result.push_back(oneResult);
				oneResult = vg::Alignment {};
				continue;
			}
			auto connector = connectors.at(pair);
			if (oneResult.path().mapping_size() == 0)
			{
				auto mapping = oneResult.mutable_path()->add_mapping();
				mapping->mutable_position()->set_node_id(pair.first.id);
				mapping->mutable_position()->set_is_reverse(pair.first.end);
			}
			if (!connector.second)
			{
				auto mapping = oneResult.mutable_path()->add_mapping();
				mapping->mutable_position()->set_node_id(connector.first.id);
				mapping->mutable_position()->set_is_reverse(connector.first.end);
			}
			auto mapping = oneResult.mutable_path()->add_mapping();
			mapping->mutable_position()->set_node_id(pair.second.id);
			mapping->mutable_position()->set_is_reverse(pair.second.end);
		}
		if (oneResult.path().mapping_size() > 0) result.push_back(oneResult);
	}
	std::cout << result.size() << " parts" << std::endl;
	return result;
}

std::vector<vg::Alignment> connectorsToAlns(const std::unordered_map<std::pair<NodePos, NodePos>, std::pair<NodePos, bool>>& connectors)
{
	std::vector<vg::Alignment> result;
	for (auto pair : connectors)
	{
		vg::Alignment oneResult;
		auto mapping = oneResult.mutable_path()->add_mapping();
		mapping->mutable_position()->set_node_id(pair.first.first.id);
		mapping->mutable_position()->set_is_reverse(pair.first.first.end);
		if (!pair.second.second)
		{
			mapping = oneResult.mutable_path()->add_mapping();
			mapping->mutable_position()->set_node_id(pair.second.first.id);
			mapping->mutable_position()->set_is_reverse(pair.second.first.end);
		}
		mapping = oneResult.mutable_path()->add_mapping();
		mapping->mutable_position()->set_node_id(pair.first.second.id);
		mapping->mutable_position()->set_is_reverse(pair.first.second.end);
		result.push_back(oneResult);
	}
	std::cout << result.size() << " connectors" << std::endl;
	return result;
}

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string inputAlns { argv[2] };
	std::string outputAlns { argv[3] };

	auto graph = GfaGraph::LoadFromFile(graphFile);
	auto alns = CommonUtils::LoadVGAlignments(inputAlns);
	auto existingNodes = getExistingNodes(graph, alns);
	auto connectors = getConnectors(graph, existingNodes);
	auto newAlns = remapAlignments(alns, existingNodes, connectors);
	auto connectorAlns = connectorsToAlns(connectors);

	std::vector<vg::Alignment> result { newAlns.begin(), newAlns.end() };
	result.insert(result.end(), connectorAlns.begin(), connectorAlns.end());
	std::ofstream alignmentOut { outputAlns, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, result, 0);
}