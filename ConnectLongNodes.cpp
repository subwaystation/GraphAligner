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

std::vector<std::tuple<NodePos, NodePos, bool, std::string>> getSecondaryConnectors(const std::vector<vg::Alignment>& alns, const std::vector<FastQ>& reads, const std::unordered_set<int>& longNodes, size_t minLen, const GfaGraph& graph)
{
	std::vector<std::tuple<NodePos, NodePos, bool, std::string>> result;
	std::unordered_map<std::string, std::string> readSeqs;
	for (auto read : reads)
	{
		readSeqs[read.seq_id] = read.sequence;
	}
	std::unordered_map<std::string, std::vector<vg::Alignment>> alnsPerRead;
	for (auto aln : alns)
	{
		alnsPerRead[aln.name()].push_back(aln);
	}
	for (auto pair : alnsPerRead)
	{
		auto readseq = readSeqs[pair.first];
		std::sort(pair.second.begin(), pair.second.end(), [](const vg::Alignment& left, const vg::Alignment& right) { return left.query_position() < right.query_position(); });
		std::vector<std::pair<NodePos, size_t>> alnFirstValid;
		std::vector<std::pair<NodePos, size_t>> alnLastValid;
		alnFirstValid.resize(pair.second.size(), std::make_pair(NodePos{-1, false}, -1));
		alnLastValid.resize(pair.second.size(), std::make_pair(NodePos{-1, false}, -1));
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (pair.second[i].sequence().size() < minLen) continue;
			for (int j = 0; j < pair.second[i].path().mapping_size(); j++)
			{
				if (longNodes.count(pair.second[i].path().mapping(j).position().node_id()) == 1)
				{
					if (pair.second[i].path().mapping(j).position().offset() <= 64)
					{
						alnFirstValid[i].first.id = pair.second[i].path().mapping(j).position().node_id();
						alnFirstValid[i].first.end = pair.second[i].path().mapping(j).position().is_reverse();
						alnFirstValid[i].second = pair.second[i].query_position();
						break;
					}
				}
			}
			for (int j = pair.second[i].path().mapping_size()-1; j >= 0; j--)
			{
				if (longNodes.count(pair.second[i].path().mapping(j).position().node_id()) == 1)
				{
					if (pair.second[i].path().mapping(j).position().offset() + pair.second[i].path().mapping(j).edit(0).from_length() >= graph.nodes.at(pair.second[i].path().mapping(j).position().node_id()).size()-64)
					{
						alnLastValid[i].first.id = pair.second[i].path().mapping(j).position().node_id();
						alnLastValid[i].first.end = pair.second[i].path().mapping(j).position().is_reverse();
						alnLastValid[i].second = pair.second[i].query_position() + pair.second[i].sequence().size();
						break;
					}
				}
			}
		}
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (alnFirstValid[i].first.id == -1) continue;
			for (size_t j = i-1; j < pair.second.size(); j--)
			{
				if (alnLastValid[j].first.id == -1) continue;
				if (alnLastValid[j].first.id == alnFirstValid[i].first.id) continue;
				if (alnLastValid[j].second >= alnFirstValid[i].second) continue;
				result.emplace_back(alnLastValid[j].first, alnFirstValid[i].first, false, readseq.substr(alnLastValid[j].second, alnFirstValid[i].second - alnLastValid[j].second));
				break;
			}
		}
	}
	return result;
}

std::vector<std::tuple<NodePos, NodePos, bool, std::string>> pickPrimaryAndSecondaryConnectors(const std::vector<std::tuple<NodePos, NodePos, bool, std::string>>& primaries, const std::vector<std::tuple<NodePos, NodePos, bool, std::string>>& secondaries)
{
	std::unordered_set<int> hasPrimaryRightEdge;
	std::unordered_set<int> hasPrimaryLeftEdge;
	for (auto connector : primaries)
	{
		if (std::get<0>(connector).end)
		{
			hasPrimaryRightEdge.insert(std::get<0>(connector).id);
		}
		else
		{
			hasPrimaryLeftEdge.insert(std::get<0>(connector).id);
		}
		if (std::get<1>(connector).end)
		{
			hasPrimaryLeftEdge.insert(std::get<1>(connector).id);
		}
		else
		{
			hasPrimaryRightEdge.insert(std::get<1>(connector).id);
		}
	}
	std::vector<std::tuple<NodePos, NodePos, bool, std::string>> result = primaries;
	for (auto connector : secondaries)
	{
		if (std::get<0>(connector).end)
		{
			if (hasPrimaryRightEdge.count(std::get<0>(connector).id) == 1) continue;
		}
		else
		{
			if (hasPrimaryLeftEdge.count(std::get<0>(connector).id) == 1) continue;
		}
		if (std::get<1>(connector).end)
		{
			if (hasPrimaryLeftEdge.count(std::get<1>(connector).id) == 1) continue;
		}
		else
		{
			if (hasPrimaryRightEdge.count(std::get<1>(connector).id) == 1) continue;
		}
		result.push_back(connector);
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string alnfile { argv[1] };
	std::string graphfile { argv[2] };
	std::string readsfile { argv[3] };
	int minLongnodeLength = std::stoi(argv[4]);
	int minLongnodeAlnlen = std::stoi(argv[5]);
	std::string outputGraphFile { argv[6] };

	auto alns = CommonUtils::LoadVGAlignments(alnfile);
	auto reads = loadFastqFromFile(readsfile);
	auto graph = GfaGraph::LoadFromFile(graphfile);
	auto longNodes = getLongNodes(graph, minLongnodeLength);
	auto parts = splitAlnsToParts(alns, longNodes, minLongnodeAlnlen);
	auto connectors = getConnectingNodes(parts, longNodes, graph);
	auto canon = canonizePaths(connectors);
	auto uniques = pickUniquePaths(canon);
	auto secondaries = getSecondaryConnectors(alns, reads, longNodes, minLongnodeAlnlen, graph);
	auto canonSecondaries = canonizePaths(secondaries);
	auto merged = pickPrimaryAndSecondaryConnectors(uniques, canonSecondaries);
	makeGraphAndWrite(merged, longNodes, graph, outputGraphFile);
}