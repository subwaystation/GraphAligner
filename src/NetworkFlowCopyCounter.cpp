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
	std::unordered_map<int, size_t> nodeStart;
	size_t maxNode;
};

class CopycountConfidence
{
public:
	std::unordered_map<int, size_t> copyCount;
	std::unordered_map<int, double> equalConfidence;
	std::unordered_map<int, double> biggerConfidence;
	std::unordered_map<int, double> smallerConfidence;
};

double contiguityBreakLogOdds()
{
	//arbitrarily say 2 contiguity break per 3Gbp
	return (log(2.0) - log(3000000000.0));
}

double contiguityBreakScore()
{
	return -contiguityBreakLogOdds();
}

double copyCountLogOdd(double relativeCoverage)
{
	//assume erlang distribution, return log, simplify
	constexpr double k = 2;
	constexpr double lambda = 1;
	if (relativeCoverage <= 0.01) return -5;
	double logpmf = (log(relativeCoverage) * (k-1) - lambda * relativeCoverage);
	return logpmf;
}

double lengthLogistic(double length)
{
	//arbitrarily weight nodes based on length, short unreliable, approx 5000 bp reliable-ish, less than linearly for longer
	double result = 1.0 / (1.0 + pow(2.7, -0.001*length + 5)) * log(length);
	return result;
}

double copyCountScore(double relativeCoverage, double length)
{
	return -copyCountLogOdd(relativeCoverage) * lengthLogistic(length);
}

double zeroCountLogOdd(double relativeCoverage)
{
	//assume exponential distribution, return log
	constexpr double lambda = 10;
	double result = -lambda * relativeCoverage + log(lambda);
	return result;
}

double zeroCountScore(double relativeCoverage, double length)
{
	return -zeroCountLogOdd(relativeCoverage) * lengthLogistic(length);
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

double copyCountLogOdd(size_t copyCount, double onecopyCoverage, double coverage, size_t length)
{
	if (copyCount > 0)
	{
		return copyCountScore(coverage / (copyCount * onecopyCoverage), length);
	}
	return zeroCountScore(coverage / onecopyCoverage, length);
}

FlowGraph buildFlowGraph(const GfaGraph& graph, int numChromosomes, double onecopyCoverage)
{
	FlowGraph result;
	result.maxNode = 0;
	size_t nextNode = 0;
	for (auto node : graph.nodes)
	{
		result.nodeStart[node.first] = nextNode;
		size_t kmers;
		size_t length;
		std::tie(kmers, length) = getNodeKmersLength(graph, node.first);
		double coverage = (double)kmers / (double)length;
		assert(coverage > 0);
		result.edges.emplace_back(NodePos { -1, true }, -1, nextNode * 4, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, -1, nextNode * 4 + 2, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, nextNode * 4 + 1, -2, 0, 100, 0);
		result.edges.emplace_back(NodePos { -1, true }, nextNode * 4 + 3, -2, 0, 100, 0);
		for (size_t copyCount = 1; copyCount < coverage / onecopyCoverage * 2 + 3; copyCount += 1)
		{
			double score = copyCountLogOdd(copyCount, onecopyCoverage, coverage, length) - copyCountLogOdd(copyCount-1, onecopyCoverage, coverage, length);
			assert(!std::isnan(score));
			assert(std::isfinite(score));
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
			result.edges.emplace_back(NodePos { -2, true }, result.nodeStart[edge.first.id] * 4 + 1 + (edge.first.end ? 0 : 2), result.nodeStart[target.id] * 4 + (target.end ? 0 : 2), 0, 1000000, 0);
			result.edges.emplace_back(NodePos { -2, true }, result.nodeStart[target.id] * 4 + 1 + (target.end ? 2 : 0), result.nodeStart[edge.first.id] * 4 + (edge.first.end ? 2 : 0), 0, 1000000, 0);
		}
	}
	result.edges.emplace_back(NodePos { -3, true }, -2, -3, -contiguityBreakScore(), numChromosomes * 2, 0); //*2 because of forward and reverse flows
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
std::vector<std::pair<size_t, bool>> getNegativeCycle(const FlowGraph& graph, const size_t start, const size_t end)
{
	constexpr double delta = 0.0001;
	static OffsetedVector<double> distances;
	static OffsetedVector<std::pair<size_t, bool>> from;
	static std::vector<size_t> posInPath;
	distances.assign(graph.maxNode, 3, std::numeric_limits<double>::infinity());
	from.assign(graph.maxNode, 3, std::pair<size_t, bool> { -5, true });
	posInPath.assign(graph.edges.size(), std::numeric_limits<size_t>::max());
	distances[start] = 0;
	size_t iter = 0;
	for (iter = 0; iter < graph.maxNode+5; iter++)
	{
		bool changed = false;
		for (size_t i = 0; i < graph.edges.size(); i++)
		{
			if ((graph.edges[i].from == start && graph.edges[i].to == end) || (graph.edges[i].from == end && graph.edges[i].to == start)) continue;
			if (graph.edges[i].used < graph.edges[i].capacity && distances[graph.edges[i].from ]!= std::numeric_limits<double>::infinity())
			{
				if (distances[graph.edges[i].to] > distances[graph.edges[i].from] + graph.edges[i].cost + delta)
				{
					from[graph.edges[i].to] = std::make_pair(i, true);
					distances[graph.edges[i].to] = distances[graph.edges[i].from] + graph.edges[i].cost;
					changed = true;
				}
			}
			if (graph.edges[i].used > 0 && distances[graph.edges[i].to] != std::numeric_limits<double>::infinity())
			{
				if (distances[graph.edges[i].from] > distances[graph.edges[i].to] - graph.edges[i].cost + delta)
				{
					from[graph.edges[i].from] = std::make_pair(i, false);
					distances[graph.edges[i].from] = distances[graph.edges[i].to] - graph.edges[i].cost;
					changed = true;
				}
			}
		}
		if (!changed)
		{
			std::cout << "distances converged ";
			break;
		}
		if (distances[end] < 0)
		{
			std::cout << "negative path to end found " << distances[end] << " ";
			break;
		}
		if (distances[start] < 0)
		{
			std::cout << "negative loop to start found " << distances[start] << " ";
			break;
		}
	}
	std::cout << "Bellman-Ford took " << iter << " iterations" << std::endl;
	
	std::cout << "find path" << std::endl;
	std::vector<std::pair<size_t, bool>> result;
	if (distances[end] < 0)
	{
		std::cout << "path from end" << std::endl;
		result.emplace_back(from[end]);
	}
	else if (distances[start] < 0)
	{
		std::cout << "loop from start" << std::endl;
		result.emplace_back(from[start]);
	}
	else
	{
		for (size_t i = 0; i < graph.edges.size(); i++)
		{
			if ((graph.edges[i].from == start && graph.edges[i].to == end) || (graph.edges[i].from == end && graph.edges[i].to == start)) continue;
			if (graph.edges[i].used < graph.edges[i].capacity)
			{
				if (distances[graph.edges[i].to] > distances[graph.edges[i].from] + graph.edges[i].cost + delta)
				{
					from[graph.edges[i].to] = std::make_pair(i, true);
					distances[graph.edges[i].to] = distances[graph.edges[i].from] + graph.edges[i].cost;
					std::cout << "cycle from edge " << i << "+" << std::endl;
					result.emplace_back(i, true);
					break;
				}
			}
			if (graph.edges[i].used > 0)
			{
				if (distances[graph.edges[i].from] > distances[graph.edges[i].to] - graph.edges[i].cost + delta)
				{
					from[graph.edges[i].from] = std::make_pair(i, false);
					distances[graph.edges[i].from] = distances[graph.edges[i].to] - graph.edges[i].cost;
					std::cout << "cycle from edge " << i << "-" << std::endl;
					result.emplace_back(i, false);
					break;
				}
			}
		}
	}
	if (result.size() == 0)
	{
		std::cout << "positive weight path from end" << std::endl;
		result.emplace_back(from[end]);
	}
	assert(result.size() == 1);
	size_t lastNode = graph.edges[result[0].first].from;
	if (result[0].second) lastNode = graph.edges[result[0].first].to;
	posInPath[result[0].first] = 0;
	while (true)
	{
		if (result.back().second && graph.edges[result.back().first].from == start) break;
		if (!result.back().second && graph.edges[result.back().first].to == start) break;
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
	assert(firstNode == lastNode || (firstNode == start && lastNode == end));
	assert(result.size() >= 2);
	if (firstNode == start && lastNode == end)
	{
		double pathScore = 0;
		for (auto edge : result)
		{
			if (edge.second) pathScore += graph.edges[edge.first].cost;
			if (!edge.second) pathScore -= graph.edges[edge.first].cost;
		}
		assert(pathScore < 0 == distances[end] < 0);
	}
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
				cost -= graph.edges[pair.first].cost;
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

void outputFlow(const FlowGraph& flowGraph, const GfaGraph& graph, const CopycountConfidence& confidences, const std::string& filename)
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
	file << "node,count,prob,smallerprob,biggerprob,kmers,length,coverage" << std::endl;
	for (auto node : graph.nodes)
	{
		size_t kmers, length;
		std::tie(kmers, length) = getNodeKmersLength(graph, node.first);
		size_t forward = flowCounts[NodePos { node.first, true }];
		size_t backward = flowCounts[NodePos { node.first, false }];
		file << node.first << "," << confidences.copyCount.at(node.first) / 2 << "," << confidences.equalConfidence.at(node.first) << ",";
		if (confidences.smallerConfidence.count(node.first) == 1)
		{
			file << confidences.smallerConfidence.at(node.first) << ",";
		}
		else
		{
			file << "-" << ",";
		}
		if (confidences.biggerConfidence.count(node.first) == 1)
		{
			file << confidences.biggerConfidence.at(node.first) << ",";
		}
		else
		{
			file << "-" << ",";
		}
		file << kmers << "," << length << "," << (double)kmers/(double)length << std::endl;
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
		totalKmers += kmers;
		totalLength += length * copyCount;
	}
	return totalKmers / totalLength;
}

double shortestPathCost(const FlowGraph& graph, const size_t start, const size_t end)
{
	auto path = getNegativeCycle(graph, start, end);
	assert(path.size() >= 2);
	size_t lastNode = graph.edges[path[0].first].from;
	if (path[0].second) lastNode = graph.edges[path[0].first].to;
	size_t firstNode = graph.edges[path.back().first].to;
	if (path.back().second) firstNode = graph.edges[path.back().first].from;
	assert(lastNode == end);
	assert(firstNode == start);
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
			cost -= graph.edges[pair.first].cost;
		}
	}
	assert(available >= 1);
	return cost;
}

CopycountConfidence getCopycountConfidences(const FlowGraph& graph)
{
	std::unordered_map<NodePos, size_t> nodeStart;
	std::unordered_map<NodePos, size_t> nodeEnd;
	std::unordered_map<NodePos, double> minBigger;
	std::unordered_map<NodePos, double> minSmaller;
	CopycountConfidence result;
	for (auto edge : graph.edges)
	{
		if (edge.originalNode.id < 0) continue;
		assert(nodeStart.count(edge.originalNode) == 0 || nodeStart.at(edge.originalNode) == edge.from);
		assert(nodeEnd.count(edge.originalNode) == 0 || nodeEnd.at(edge.originalNode) == edge.to);
		nodeStart[edge.originalNode] = edge.from;
		nodeEnd[edge.originalNode] = edge.to;
		result.copyCount[edge.originalNode.id] += edge.used;
		if (edge.used < edge.capacity)
		{
			if (minBigger.count(edge.originalNode) == 0 || edge.cost < minBigger.at(edge.originalNode))
			{
				minBigger[edge.originalNode] = edge.cost;
			}
		}
		if (edge.used > 0)
		{
			if (minSmaller.count(edge.originalNode) == 0 || -edge.cost < minSmaller.at(edge.originalNode))
			{
				minSmaller[edge.originalNode] = -edge.cost;
			}
		}
	}
	for (auto node : nodeStart)
	{
		assert(nodeEnd.count(node.first) == 1);
		if (minBigger.count(node.first) == 1)
		{
			auto fwCost = shortestPathCost(graph, nodeEnd.at(node.first), node.second);
			// assert(fwCost >= -0.01);
			minBigger[node.first] += fwCost;
			assert(minBigger[node.first] >= -0.01);
		}
		if (minSmaller.count(node.first) == 1)
		{
			auto bwCost = shortestPathCost(graph, node.second, nodeEnd.at(node.first));
			// assert(bwCost >= -0.01);
			minSmaller[node.first] += bwCost;
			assert(minSmaller[node.first] >= -0.01);
		}
	}
	for (auto node : result.copyCount)
	{
		result.equalConfidence[node.first] = 0;
		bool hasSmaller = false;
		bool hasBigger = false;
		if (minSmaller.count(NodePos { node.first, true }) == 1)
		{
			result.smallerConfidence[node.first] = -minSmaller.at(NodePos { node.first, true });
			hasSmaller = true;
		}
		if (minSmaller.count(NodePos { node.first, false }) == 1)
		{
			double newSmaller = -minSmaller.at(NodePos { node.first, false });
			if (hasSmaller)
			{
				result.smallerConfidence[node.first] = std::min(result.smallerConfidence[node.first], newSmaller);
			}
			else
			{
				result.smallerConfidence[node.first] = newSmaller;
			}
		}
		if (minBigger.count(NodePos { node.first, true }) == 1)
		{
			result.biggerConfidence[node.first] = -minBigger.at(NodePos { node.first, true });
			hasBigger = true;
		}
		if (minBigger.count(NodePos { node.first, false }) == 1)
		{
			double newSmaller = -minBigger.at(NodePos { node.first, false });
			if (hasBigger)
			{
				result.biggerConfidence[node.first] = std::min(result.biggerConfidence[node.first], newSmaller);
			}
			else
			{
				result.biggerConfidence[node.first] = newSmaller;
			}
		}
		assert(hasSmaller || hasBigger);
		double totalConfidence = pow(2.7, result.equalConfidence[node.first]);
		if (hasSmaller) totalConfidence += pow(2.7, result.smallerConfidence[node.first]);
		if (hasBigger) totalConfidence += pow(2.7, result.biggerConfidence[node.first]);
		result.equalConfidence[node.first] = pow(2.7, result.equalConfidence[node.first]) / totalConfidence;
		if (hasSmaller) result.smallerConfidence[node.first] = pow(2.7, result.smallerConfidence[node.first]) / totalConfidence;
		if (hasBigger) result.biggerConfidence[node.first] = pow(2.7, result.biggerConfidence[node.first]) / totalConfidence;
	}
	return result;
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
	auto copycountConfidences = getCopycountConfidences(flowGraph);
	outputFlow(flowGraph, graph, copycountConfidences, outputCsv);
}
