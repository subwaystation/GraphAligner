#include <string>
#include <cassert>
#include <unordered_map>
#include "GfaGraph.h"


GfaGraph splitGraph(const GfaGraph& input, size_t maxNodeSize)
{
	std::unordered_map<int, int> newStart;
	std::unordered_map<int, int> newEnd;
	int nextNode = 1;
	assert(maxNodeSize > input.edgeOverlap);
	GfaGraph result;
	result.edgeOverlap = input.edgeOverlap;

	for (auto pair : input.nodes)
	{
		int id = pair.first;
		std::string seq = pair.second;
		if (seq.size() < maxNodeSize)
		{
			newStart[id] = nextNode;
			newEnd[id] = nextNode;
			result.nodes[nextNode] = seq;
			nextNode += 1;
			continue;
		}
		size_t kmerCount = seq.size() - input.edgeOverlap;
		int splitNodeCount = (kmerCount + maxNodeSize - 1) / maxNodeSize;
		assert(splitNodeCount >= 1);
		size_t splitNodeSize = (kmerCount + splitNodeCount - 1) / splitNodeCount;
		assert(splitNodeSize * splitNodeCount + input.edgeOverlap >= seq.size());
		for (size_t i = 0; i < splitNodeCount; i++)
		{
			std::string splitSeq = seq.substr(i * splitNodeSize, splitNodeSize + input.edgeOverlap);
			assert(splitSeq.size() > 0);
			if (i == 0) newStart[id] = nextNode;
			if (i == splitNodeCount - 1) newEnd[id] = nextNode;
			if (i > 0) result.edges[NodePos { nextNode-1, true }].push_back(NodePos { nextNode, true });
			result.nodes[nextNode] = splitSeq;
			nextNode += 1;
		}
	}
	for (auto pair : input.edges)
	{
		auto source = pair.first;
		assert(newStart.count(source.id) == 1);
		assert(newEnd.count(source.id) == 1);
		for (auto target : pair.second)
		{
			assert(newStart.count(target.id) == 1);
			assert(newEnd.count(target.id) == 1);
			NodePos from, to;
			if (source.end)
			{
				from = NodePos { newEnd.at(source.id), true };
			}
			else
			{
				from = NodePos { newStart.at(source.id), false };
			}
			if (target.end)
			{
				to = NodePos { newStart.at(target.id), true };
			}
			else
			{
				to = NodePos { newEnd.at(target.id), false };
			}
			assert(result.nodes.count(from.id) == 1);
			assert(result.nodes.count(to.id) == 1);
			result.edges[from].push_back(to);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	size_t maxNodeSize = std::stol(argv[2]);
	std::string outputFile { argv[3] };

	auto graph = GfaGraph::LoadFromFile(graphFile);
	assert(maxNodeSize > graph.edgeOverlap);
	graph = splitGraph(graph, maxNodeSize);
	graph.SaveToFile(outputFile);
}