#include <iostream>
#include <algorithm>
#include <fstream>
#include "fastqloader.h"
#include "GfaGraph.h"
#include "CommonUtils.h"

bool operator>(const NodePos& left, const NodePos& right)
{
	return left.id > right.id || (left.id == right.id && left.end && !right.end);
}

bool operator<(const NodePos& left, const NodePos& right)
{
	return left.id < right.id || (left.id == right.id && !left.end && right.end);
}

// std::vector<NodePos> topoSort(const GfaGraph& graph)
// {
// 	std::set<NodePos> possibleStarts;
// 	for (auto node : graph)
// 	{
// 		if (graph.edges.count(NodePos{node.first, false}) == 0) edgelessEnds.emplace(node.first, true);
// 		if (graph.edges.count(NodePos{node.first, true}) == 0) edgelessEnds.emplace(node.first, false);
// 	}
// 	std::vector<NodePos> result;
// 	std::unordered_map<int, bool> setOrder;
// 	while (true)
// 	{
// 		assert(possibleStarts.size() > 0);
// 		auto start = *possibleStarts.begin();
// 		possibleStarts.erase(start);
// 		if (setOrder.count(start.id) == 1) continue;
// 		std::queue<NodePos> expandables;
// 		expandables.push(start);
// 		while (expandables.size() > 0)
// 		{
// 			auto top = expandables.front();
// 			expandables.pop();
// 			assert(setOrder.count(top.id) == 0);
// 			result.push_back(top);
// 			setOrder[top.id] = top.end;
// 			for (auto neighbor : graph.edges[top])
// 			{
// 				expandables.push(neighbor);
// 			}
// 		}
// 	}
// 	assert(setOrder.size() == result.size());
// 	assert(setOrder.size() == graph.nodes.size());
// 	return result;
// }

// std::unordered_map<int, size_t> getOrderIndex(const std::vector<NodePos>& order)
// {
// 	std::unordered_map<int, size_t> result;
// 	for (size_t i = 0; i < order.size(); i++)
// 	{
// 		result[order[i].id] = i;
// 	}
// 	return result;
// }

void writeUnitigs(const GfaGraph& graph, std::string filename)
{
	std::unordered_map<NodePos, NodePos> nodeToUnitig;
	std::vector<std::string> unitigSequences;
	for (auto node : graph.nodes)
	{
		auto leftEnd = NodePos { node.first, false };
		auto rightEnd = NodePos { node.first, true };
		if (nodeToUnitig.count(leftEnd) == 1) continue;
		if (nodeToUnitig.count(rightEnd) == 1) continue;
		NodePos pos;
		if (graph.edges.count(leftEnd) == 0 || graph.edges.at(leftEnd).size() != 1 || graph.edges.at(graph.edges.at(leftEnd)[0].Reverse()).size() != 1)
		{
			pos = rightEnd;
		}
		else if (graph.edges.count(rightEnd) == 0 || graph.edges.at(rightEnd).size() != 1 || graph.edges.at(graph.edges.at(rightEnd)[0].Reverse()).size() != 1)
		{
			pos = leftEnd;
		}
		else
		{
			continue;
		}
		nodeToUnitig[pos].id = unitigSequences.size();
		nodeToUnitig[pos].end = true;
		std::string sequence = node.second;
		if (!pos.end) sequence = CommonUtils::ReverseComplement(sequence);
		while (graph.edges.count(pos) == 1 && graph.edges.at(pos).size() == 1)
		{
			pos = graph.edges.at(pos)[0];
			NodePos reversePos = pos.Reverse();
			assert(graph.edges.count(reversePos) == 1);
			if (graph.edges.at(reversePos).size() != 1) break;
			nodeToUnitig[pos].id = unitigSequences.size();
			nodeToUnitig[pos].end = true;
			std::string sequenceHere = graph.nodes.at(pos.id);
			if (!pos.end) sequenceHere = CommonUtils::ReverseComplement(sequenceHere);
			sequence += sequenceHere.substr(graph.edgeOverlap);
		}
		unitigSequences.push_back(sequence);
	}
	for (auto node : graph.nodes)
	{
		auto leftEnd = NodePos { node.first, true };
		auto rightEnd = NodePos { node.first, false };
		assert(nodeToUnitig.count(leftEnd) == 1 || nodeToUnitig.count(rightEnd) == 1);
		assert(!(nodeToUnitig.count(leftEnd) == 1 && nodeToUnitig.count(rightEnd) == 1));
	}
	std::ofstream file { filename };
	for (size_t i = 0; i < unitigSequences.size(); i++)
	{
		file << "S\t" << i << "\t" << unitigSequences[i] << std::endl;
	}
	for (auto edge : graph.edges)
	{
		NodePos unitigFrom;
		if (nodeToUnitig.count(edge.first) == 1)
		{
			unitigFrom = nodeToUnitig[edge.first];
		}
		else
		{
			assert(nodeToUnitig.count(edge.first.Reverse()) == 1);
			unitigFrom = nodeToUnitig[edge.first.Reverse()].Reverse();
		}
		for (auto target : edge.second)
		{
			NodePos unitigTo;
			if (nodeToUnitig.count(target) == 1)
			{
				unitigTo = nodeToUnitig[target];
			}
			else
			{
				assert(nodeToUnitig.count(target.Reverse()) == 1);
				unitigTo = nodeToUnitig[target.Reverse()].Reverse();
			}
			if (unitigFrom != unitigTo)
			{
				file << "L\t" << unitigFrom.id << "\t" << (unitigFrom.end ? "-" : "+") << "\t" << unitigTo.id << "\t" << (unitigTo.end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
			}
		}
	}
}

int main(int argc, char** argv)
{
	std::string inGraphFile { argv[1] };
	std::string outGraphFile { argv[2] };

	auto graph = GfaGraph::LoadFromFile(inGraphFile);
	graph.confirmDoublesidedEdges();
	for (auto edge : graph.edges)
	{
		NodePos revSource = edge.first.Reverse();
		for (auto target : edge.second)
		{
			NodePos revTarget = target.Reverse();
			assert(graph.edges.count(revTarget) == 1);
			bool found = false;
			for (auto check : graph.edges.at(revTarget))
			{
				if (check == revSource) found = true;
			}
			assert(found);
		}
	}
	// auto order = topoSort(graph);
	// auto revOrder = getOrderIndex(order);
	writeUnitigs(graph, outGraphFile);
}
