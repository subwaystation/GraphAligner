#include <fstream>
#include "GfaGraph.h"
#include "CommonUtils.h"

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

void makeGraphAndWrite(const std::vector<std::vector<NodePos>>& parts, const std::unordered_set<int> longNodes, const GfaGraph& graph, std::string filename)
{
	std::ofstream file { filename };
	for (auto node : longNodes)
	{
		file << "S\t" << node << "\t" << graph.nodes.at(node) << std::endl;
	}
	int nextId = graph.nodes.size()+1;
	for (auto path : parts)
	{
		assert(path.size() >= 2);
		assert(longNodes.count(path[0].id) == 1);
		assert(longNodes.count(path.back().id) == 1);
		if (path.size() == 2)
		{
			file << "L\t" << path[0].id << "\t" << (path[0].end ? "-" : "+") << "\t" << path[1].id << "\t" << (path[1].end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
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
		file << "S\t" << nextId << "\t" << sequence << std::endl;
		file << "L\t" << path[0].id << "\t" << (path[0].end ? "-" : "+") << "\t" << nextId << "\t+\t" << graph.edgeOverlap << "M" << std::endl;
		file << "L\t" << nextId << "\t+\t" << path.back().id << "\t" << (path.back().end ? "-" : "+") << "\t" << graph.edgeOverlap << "M" << std::endl;
		nextId++;
	}
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
	makeGraphAndWrite(parts, longNodes, graph, outputGraphFile);
}