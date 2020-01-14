#include <limits>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include "GfaGraph.h"

class FlowGraph
{
public:
	struct Edge
	{
		Edge(NodePos originalNode, size_t from, size_t to, double cost, size_t capacity, size_t used) :
			originalNode(originalNode),
			from(from),
			to(to),
			cost(cost),
			capacity(capacity),
			used(used)
		{}
		NodePos originalNode;
		size_t from;
		size_t to;
		double cost;
		size_t capacity;
		size_t used;
	};
	std::vector<Edge> edges;
	size_t maxNode;
};

double contiguityBreakProbability()
{
	return 2.0 / 4600000.0; //approx from e. coli
}

double contiguityBreakScore()
{
	return -log(contiguityBreakProbability());
}

double copyCountLogOdd(double relativeCoverage)
{
	//assume erlang distribution, return log, simplify
	constexpr double k = 2;
	constexpr double lambda = 1;
	if (relativeCoverage <= 0.01) return -5;
	double logpmf = (log(relativeCoverage) * (k-1) - lambda * relativeCoverage);
	// if (logpmf < -5) logpmf = -5;
	return logpmf;
}

double copyCountScore(double relativeCoverage, double length)
{
	//arbitrarily take sqrt of length to flatten distribution
	return -copyCountLogOdd(relativeCoverage) * pow(length, 0.5);
}

std::pair<size_t, size_t> getNodeKmersLength(const GfaGraph& graph, int nodeId)
{
	assert(graph.tags.count(nodeId) == 1);
	std::stringstream sstr { graph.tags.at(nodeId) };
	size_t kmers = 0;
	size_t size = 0;
	while (sstr.good())
	{
		std::string tag;
		sstr >> tag;
		if (tag.substr(0, 5) == "RC:i:")
		{
			assert(kmers == 0);
			kmers = std::stol(tag.substr(5));
		}
		if (tag.substr(0, 5) == "LN:i:")
		{
			assert(size == 0);
			size = std::stol(tag.substr(5));
		}
	}
	assert(kmers > 0);
	assert(size > 0);
	return std::make_pair(kmers, size);
}

FlowGraph buildFlowGraph(const GfaGraph& graph, int numChromosomes, double onecopyCoverage)
{
	FlowGraph result;
	result.maxNode = 0;
	size_t nextNode = 0;
	std::unordered_map<int, size_t> nodeStart;
	for (auto node : graph.nodes)
	{
		nodeStart[node.first] = nextNode;
		size_t kmers;
		size_t length;
		std::tie(kmers, length) = getNodeKmersLength(graph, node.first);
		double coverage = (double)kmers / (double)length;
		assert(coverage > 0);
		result.edges.emplace_back(NodePos { -1, true }, -1, nextNode * 4, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, -1, nextNode * 4 + 2, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, nextNode * 4 + 1, -2, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, nextNode * 4 + 3, -2, 0, 100, 0);
		for (size_t i = onecopyCoverage; i < coverage * 2 + 3 * onecopyCoverage; i += onecopyCoverage)
		{
			double score = copyCountScore((double)i / coverage, length) - copyCountScore((double)(i-onecopyCoverage) / coverage, length);
			result.edges.emplace_back(NodePos { node.first, true }, nextNode * 4, nextNode * 4 + 1, score, 1, 0);
			result.edges.emplace_back(NodePos { node.first, false }, nextNode * 4 + 2, nextNode * 4 + 3, score, 1, 0);
		}
		result.maxNode = std::max(result.maxNode, (size_t)nextNode * 4 + 4);
		nextNode += 1;
	}
	for (auto edge : graph.edges)
	{
		for (auto target : edge.second)
		{
			result.edges.emplace_back(NodePos { -2, true }, nodeStart[edge.first.id] * 4 + 1 + (edge.first.end ? 0 : 2), nodeStart[target.id] * 4 + (target.end ? 0 : 2), 0, 1000000, 0);
			result.edges.emplace_back(NodePos { -2, true }, nodeStart[target.id] * 4 + 1 + (target.end ? 2 : 0), nodeStart[edge.first.id] * 4 + (edge.first.end ? 2 : 0), 0, 1000000, 0);
		}
	}
	result.edges.emplace_back(NodePos { -3, true }, -2, -3, 0, numChromosomes, 0);
	result.edges.emplace_back(NodePos { -4, true }, -2, -3, contiguityBreakScore(), 1000000, 0);
	result.edges.emplace_back(NodePos { -4, true }, -3, -1, 0, 1000000, 0);
	return result;
}

template <typename T>
class OffsetedVector
{
public:
	T& operator[](size_t index)
	{
		return data[index + offset];
	}
	void assign(size_t amount, size_t newoffset, T item)
	{
		offset = newoffset;
		data.assign(amount + newoffset, item);
	}
	auto begin()
	{
		return data.begin();
	}
	auto end()
	{
		return data.end();
	}
private:
	std::vector<T> data;
	size_t offset;
};

// https://en.wikipedia.org/wiki/Bellman%E2%80%93Ford_algorithm
std::vector<std::pair<size_t, bool>> getNegativeCycle(const FlowGraph& graph, size_t start, size_t end)
{
	static OffsetedVector<double> distances;
	static OffsetedVector<std::pair<size_t, bool>> from;
	static std::vector<size_t> posInPath;
	distances.assign(graph.maxNode, 3, std::numeric_limits<double>::max());
	from.assign(graph.maxNode, 3, std::pair<size_t, bool> { -5, true });
	posInPath.assign(graph.edges.size(), std::numeric_limits<size_t>::max());
	distances[start] = 0;
	for (size_t iter = 0; iter < graph.maxNode+5; iter++)
	{
		if (iter % 100 == 0) std::cout << "Bellman-Ford iteration " << iter << "/" << graph.maxNode << std::endl;
		bool changed = false;
		for (size_t i = 0; i < graph.edges.size(); i++)
		{
			if (graph.edges[i].used < graph.edges[i].capacity)
			{
				if (distances[graph.edges[i].to] > distances[graph.edges[i].from] + graph.edges[i].cost)
				{
					from[graph.edges[i].to] = std::make_pair(i, true);
					distances[graph.edges[i].to] = distances[graph.edges[i].from] + graph.edges[i].cost;
					changed = true;
				}
			}
			if (graph.edges[i].used > 0)
			{
				if (distances[graph.edges[i].from] > distances[graph.edges[i].to] - graph.edges[i].cost)
				{
					from[graph.edges[i].from] = std::make_pair(i, false);
					distances[graph.edges[i].from] = distances[graph.edges[i].to] - graph.edges[i].cost;
					changed = true;
				}
			}
		}
		if (!changed) break;
	}
	
	std::cout << "find path" << std::endl;
	std::vector<std::pair<size_t, bool>> result;
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		if (graph.edges[i].used < graph.edges[i].capacity)
		{
			if (distances[graph.edges[i].to] > distances[graph.edges[i].from] + graph.edges[i].cost)
			{
				result.emplace_back(i, true);
				break;
			}
		}
		if (graph.edges[i].used > 0)
		{
			if (distances[graph.edges[i].from] > distances[graph.edges[i].to] - graph.edges[i].cost)
			{
				result.emplace_back(i, false);
				break;
			}
		}
	}
	if (result.size() == 0)
	{
		assert(distances[start] == 0);
		return result;
	}
	assert(distances[start] < 0);
	assert(result.size() == 1);
	size_t lastNode = graph.edges[result[0].first].from;
	if (result[0].second) lastNode = graph.edges[result[0].first].to;
	posInPath[result[0].first] = 0;
	while (true)
	{
		if (result.back().second && graph.edges[result.back().first].from == lastNode) break;
		if (!result.back().second && graph.edges[result.back().first].to == lastNode) break;
		if (result.back().second)
		{
			std::pair<size_t, bool> edge = from[graph.edges[result.back().first].from];
			if (posInPath[edge.first] != std::numeric_limits<size_t>::max())
			{
				result.erase(result.begin(), result.begin() + posInPath[edge.first]);
				break;
			}
			result.push_back(edge);
			posInPath[edge.first] = result.size()-1;
		}
		else
		{
			std::pair<size_t, bool> edge = from[graph.edges[result.back().first].to];
			if (posInPath[edge.first] != std::numeric_limits<size_t>::max())
			{
				result.erase(result.begin(), result.begin() + posInPath[edge.first]);
				break;
			}
			result.push_back(edge);
			posInPath[edge.first] = result.size()-1;
		}
		assert(result.size() < graph.maxNode);
	}
	lastNode = graph.edges[result[0].first].from;
	if (result[0].second) lastNode = graph.edges[result[0].first].to;
	size_t firstNode = graph.edges[result.back().first].to;
	if (result.back().second) firstNode = graph.edges[result.back().first].from;
	assert(firstNode == lastNode);
	return result;
}

void solveFlow(FlowGraph& graph)
{
	size_t iteration = 0;
	while (true)
	{
		std::cout << "flow iteration " << iteration << std::endl;
		iteration++;
		auto path = getNegativeCycle(graph, -1, -3);
		if (path.size() == 0)
		{
			std::cout << "empty path, quit" << std::endl;
			break;
		}
		assert(path.size() >= 2);
		size_t lastNode = graph.edges[path[0].first].from;
		if (path[0].second) lastNode = graph.edges[path[0].first].to;
		size_t firstNode = graph.edges[path.back().first].to;
		if (path.back().second) firstNode = graph.edges[path.back().first].from;
		if (firstNode == -1 && lastNode == -3)
		{
			assert(graph.edges.back().from == -3);
			assert(graph.edges.back().to == -1);
			path.emplace_back(graph.edges.size()-1, true);
		}
		size_t available = 1000000;
		double cost = 0;
		for (auto pair : path)
		{
			if (pair.second)
			{
				assert(graph.edges[pair.first].used < graph.edges[pair.first].capacity);
				available = std::min(available, graph.edges[pair.first].capacity - graph.edges[pair.first].used);
				cost += graph.edges[pair.first].cost;
			}
			else
			{
				assert(graph.edges[pair.first].used > 0);
				available = std::min(available, graph.edges[pair.first].used);
				cost += graph.edges[pair.first].cost;
			}
		}
		assert(available == 1);
		std::cout << "cost " << cost << std::endl;
		if (cost >= 0)
		{
			break;
		}
		for (auto pair : path)
		{
			if (pair.second)
			{
				assert(graph.edges[pair.first].used + available <= graph.edges[pair.first].capacity);
				graph.edges[pair.first].used += available;
			}
			else
			{
				assert(graph.edges[pair.first].used >= available);
				graph.edges[pair.first].used -= available;
			}
		}
	}
	std::cout << "flow solved " << iteration << std::endl;
}

void outputFlow(const FlowGraph& flowGraph, const GfaGraph& graph, const std::string& filename)
{
	std::ofstream file { filename };
	std::unordered_map<NodePos, int> flowCounts;
	for (auto edge : flowGraph.edges)
	{
		if (edge.originalNode.id >= 0)
		{
			flowCounts[edge.originalNode] += edge.used;
		}
	}
	file << "node,count,kmers,length,coverage" << std::endl;
	for (auto node : graph.nodes)
	{
		size_t kmers, length;
		std::tie(kmers, length) = getNodeKmersLength(graph, node.first);
		size_t forward = flowCounts[NodePos { node.first, true }];
		size_t backward = flowCounts[NodePos { node.first, false }];
		file << node.first << "," << (double)(forward + backward) / 2 << "," << kmers << "," << length << "," << (double)kmers/(double)length << std::endl;
	}
}

double recalcCoverage(const FlowGraph& flowGraph, const GfaGraph& graph)
{
	std::unordered_map<NodePos, int> flowCounts;
	for (auto edge : flowGraph.edges)
	{
		if (edge.originalNode.id >= 0)
		{
			flowCounts[edge.originalNode] += edge.used;
		}
	}
	double totalKmers = 0;
	double totalLength = 0;
	for (auto node : graph.nodes)
	{
		double kmers, length;
		std::tie(kmers, length) = getNodeKmersLength(graph, node.first);
		double copyCount = ((double)flowCounts[NodePos { node.first, true }] + (double)flowCounts[NodePos { node.first, false }])/2.0;
		totalKmers += kmers * copyCount;
		totalLength += length * copyCount;
	}
	return totalKmers / totalLength;
}

int main(int argc, char** argv)
{
	std::string inGraph { argv[1] };
	int numChromosomes = std::stoi(argv[2]);
	double onecopyCoverage = std::stod(argv[3]);
	std::string outputCsv { argv[4] };

	auto graph = GfaGraph::LoadFromFile(inGraph, true);
	auto flowGraph = buildFlowGraph(graph, numChromosomes, onecopyCoverage);
	solveFlow(flowGraph);
	std::cout << "recalced coverage " << recalcCoverage(flowGraph, graph) << std::endl;
	outputFlow(flowGraph, graph, outputCsv);
}
