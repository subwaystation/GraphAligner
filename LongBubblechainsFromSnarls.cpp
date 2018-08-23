#include "GfaGraph.h"
#include "stream.hpp"
#include "vg.pb.h"

std::vector<std::vector<int>> LoadSnarlTraversalsFromFile(std::string filename)
{
	std::vector<std::vector<int>> result;
	std::ifstream file { filename, std::ios::in | std::ios::binary };
	std::function<void(vg::SnarlTraversal&)> lambda = [&result](vg::SnarlTraversal& s) {
		result.emplace_back();
		for (int i = 0; i < s.visits_size(); i++)
		{
			if (s.visits(i).has_snarl())
			{
				result.back().push_back(s.visits(i).snarl().start().node_id());
				result.back().push_back(s.visits(i).snarl().end().node_id());
			}
			else
			{
				result.back().push_back(s.visits(i).node_id());
			}
		}
	};
	stream::for_each(file, lambda);
	return result;
}

int find(std::unordered_map<int, int>& parent, int node)
{
	assert(parent.count(node) == 1);
	auto found = parent[node];
	if (found != node)
	{
		found = find(parent, found);
		parent[node] = found;
	}
	return found;
}

void merge(std::unordered_map<int, int>& parent, std::unordered_map<int, int>& rank, int left, int right)
{
	auto leftroot = find(parent, left);
	auto rightroot = find(parent, right);
	if (leftroot == rightroot) return;
	if (rank[leftroot] < rank[rightroot])
	{
		std::swap(leftroot, rightroot);
	}
	parent[rightroot] = leftroot;
	if (rank[leftroot] == rank[rightroot]) rank[leftroot]++;
}

std::vector<int> findLongChains(const GfaGraph& graph, const std::vector<std::vector<int>>& traversals, int minChainLength)
{
	std::unordered_map<int, int> parent;
	std::unordered_map<int, int> rank;
	for (auto node : graph.nodes)
	{
		parent[node.first] = node.first;
		rank[node.first] = 0;
	}
	while (true)
	{
		bool repeat = false;
		for (auto path : traversals)
		{
			for (size_t i = 1; i < path.size(); i++)
			{
				assert(graph.nodes.count(path[i-1]) == 1);
				assert(graph.nodes.count(path[i]) == 1);
				if (find(parent, path[i-1]) != find(parent, path[i]))
				{
					merge(parent, rank, path[i-1], path[i]);
					repeat = true;
				}
			}
		}
		if (!repeat) break;
	}
	std::unordered_map<int, int> size;
	for (auto pair : parent)
	{
		auto node = pair.first;
		auto parentnode = find(parent, node);
		size[parentnode] += graph.nodes.at(node).size() - graph.edgeOverlap;
	}
	std::cout << "node,component" << std::endl;
	std::vector<int> result;
	for (auto pair : parent)
	{
		auto node = pair.first;
		auto parentnode = find(parent, node);
		std::cout << node << "," << parentnode << std::endl;
		if (size[parentnode] >= minChainLength)
		{
			result.push_back(node);
		}
	}
	return result;
}

void writeResults(const std::vector<int>& nodes, std::string outfile)
{
	std::ofstream file {outfile};
	for (auto node : nodes)
	{
		file << node << std::endl;
	}
}

int main(int argc, char** argv)
{
	std::string graphfile { argv[1] };
	std::string snarltraversalfile { argv[2] };
	int minChainLength = std::stoi(argv[3]);
	std::string outfile { argv[4] };

	auto graph = GfaGraph::LoadFromFile(graphfile);
	auto traversals = LoadSnarlTraversalsFromFile(snarltraversalfile);
	auto result = findLongChains(graph, traversals, minChainLength);
	writeResults(result, outfile);
}