#include <fstream>
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

std::unordered_set<int> getLongNodes(const GfaGraph& graph, int minLen)
{
	std::unordered_set<int> result;
	for (auto node : graph.nodes)
	{
		if (node.second.size() >= minLen) result.insert(node.first);
	}
	return result;
}

std::vector<std::vector<NodePos>> splitAlnsToParts(const std::vector<vg::Alignment>& alns, const std::unordered_set<int>& longNodes, int minLen)
{
	std::vector<std::vector<NodePos>> result;
	for (auto aln : alns)
	{
		std::vector<bool> isLongNode;
		isLongNode.resize(aln.path().mapping_size(), false);
		if (longNodes.count(aln.path().mapping(0).position().node_id()) == 1 && aln.path().mapping(0).edit(0).from_length() >= minLen) isLongNode[0] = true;
		for (int i = 1; i < aln.path().mapping_size() - 1; i++)
		{
			if (longNodes.count(aln.path().mapping(i).position().node_id()) == 1) isLongNode[i] = true;
		}
		if (longNodes.count(aln.path().mapping(aln.path().mapping_size()-1).position().node_id()) == 1 && aln.path().mapping(aln.path().mapping_size()-1).edit(0).from_length() >= minLen) isLongNode[aln.path().mapping_size()-1] = true;
		size_t lastLongNode = isLongNode.size();
		for (size_t i = 0; i < isLongNode.size(); i++)
		{
			if (isLongNode[i])
			{
				if (lastLongNode == isLongNode.size())
				{
					lastLongNode = i;
					continue;
				}
				result.emplace_back();
				for (size_t j = lastLongNode; j <= i; j++)
				{
					result.back().emplace_back();
					result.back().back().id = aln.path().mapping(j).position().node_id();
					result.back().back().end = aln.path().mapping(j).position().is_reverse();
				}
				lastLongNode = i;
			}
		}
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, bool, std::string>> getConnectingNodes(const std::vector<std::vector<NodePos>>& parts, const std::unordered_set<int> longNodes, const GfaGraph& graph)
{
	std::vector<std::tuple<NodePos, NodePos, bool, std::string>> result;
	for (auto path : parts)
	{
		assert(path.size() >= 2);
		assert(longNodes.count(path[0].id) == 1);
		assert(longNodes.count(path.back().id) == 1);
		if (path.size() == 2)
		{
			result.emplace_back(path[0], path.back(), true, "");
			continue;
		}
		std::string sequence;
		std::string part = graph.nodes.at(path[1].id);
		if (path[1].end) part = CommonUtils::ReverseComplement(part);
		sequence += part.substr(0, graph.edgeOverlap);
		for (size_t i = 1; i < path.size()-1; i++)
		{
			assert(longNodes.count(path[i].id) == 0);
			part = graph.nodes.at(path[i].id);
			if (path[i].end) part = CommonUtils::ReverseComplement(part);
			sequence += part.substr(graph.edgeOverlap);
		}
		result.emplace_back(path[0], path.back(), false, sequence);
	}
	return result;
}

void makeGraphAndWrite(const std::vector<std::tuple<NodePos, NodePos, bool, std::string>>& connectors, const std::unordered_set<int> longNodes, const GfaGraph& graph, std::string filename)
{
	std::ofstream file { filename };
	for (auto node : longNodes)
	{
		file << "S\t" << node << "\t" << graph.nodes.at(node) << std::endl;
	}
	int nextId = graph.nodes.size()+1;
	for (auto path : connectors)
	{
		if (std::get<2>(path))
		{
			file << "L\t" << std::get<0>(path).id << "\t" << (std::get<0>(path).end ? "-" : "+") << "\t" << std::get<1>(path).id << "\t" << (std::get<1>(path).end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
			continue;
		}
		file << "S\t" << nextId << "\t" << std::get<3>(path) << std::endl;
		file << "L\t" << std::get<0>(path).id << "\t" << (std::get<0>(path).end ? "-" : "+") << "\t" << nextId << "\t+\t" << graph.edgeOverlap << "M" << std::endl;
		file << "L\t" << nextId << "\t+\t" << std::get<1>(path).id << "\t" << (std::get<1>(path).end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
		nextId++;
	}
}

std::vector<std::tuple<NodePos, NodePos, bool, std::string>> canonizePaths(const std::vector<std::tuple<NodePos, NodePos, bool, std::string>>& paths)
{
	std::vector<std::tuple<NodePos, NodePos, bool, std::string>> result;
	for (auto path : paths)
	{
		auto start = std::get<0>(path);
		auto end = std::get<1>(path);
		if (start > end)
		{
			end.end = !end.end;
			start.end = !start.end;
			assert(!(end > start));
			result.emplace_back(end, start, std::get<2>(path), CommonUtils::ReverseComplement(std::get<3>(path)));
		}
		else
		{
			result.push_back(path);
		}
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, bool, std::string>> pickUniquePaths(const std::vector<std::tuple<NodePos, NodePos, bool, std::string>>& paths)
{
	std::set<std::tuple<NodePos, NodePos, bool, std::string>> uniques { paths.begin(), paths.end() };
	std::vector<std::tuple<NodePos, NodePos, bool, std::string>> result { uniques.begin(), uniques.end() };
	return result;
}

int main(int argc, char** argv)
{
	std::string alnfile { argv[1] };
	std::string graphfile { argv[2] };
	int minLongnodeLength = std::stoi(argv[3]);
	int minLongnodeAlnlen = std::stoi(argv[4]);
	std::string outputGraphFile { argv[5] };

	auto alns = CommonUtils::LoadVGAlignments(alnfile);
	auto graph = GfaGraph::LoadFromFile(graphfile);
	auto longNodes = getLongNodes(graph, minLongnodeLength);
	auto parts = splitAlnsToParts(alns, longNodes, minLongnodeAlnlen);
	auto connectors = getConnectingNodes(parts, longNodes, graph);
	auto canon = canonizePaths(connectors);
	auto uniques = pickUniquePaths(canon);
	makeGraphAndWrite(uniques, longNodes, graph, outputGraphFile);
}