#include <queue>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <cmath>
#include "MinimizerGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

constexpr int MAX_GRAPH_TOPO_DISTANCE = 512;
constexpr int MAX_SEQ_MINMER_DISTANCE = 512;
constexpr int SEQ_DIST_MODE_CUTOFF = 32;
constexpr int TREE_MIN_TRUNK_SIZE = 2;
constexpr int TREE_MAX_NODE_COUNT = 200000;

size_t charNum(char c)
{
	switch(c)
	{
		case 'A':
		case 'a':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'G':
		case 'g':
			return 2;
		case 'T':
		case 't':
			return 3;
	}
	assert(false);
	return 0;
}

std::string getHPC(const std::string& str)
{
	std::string result;
	result += str[0];
	for (size_t i = 1; i < str.size(); i++)
	{
		if (str[i] != str[i-1]) result += str[i];
	}
	return result;
}

MinimizerGraph::MinimizerGraph(size_t k, size_t w, const AlignmentGraph& graph) :
k(k),
w(w),
graph(graph)
{
	assert(k < w);
	assert(k < sizeof(size_t) * 8 / 2);
	if (!tryLoadOrdering("index.order.tmi"))
	{
		initRandomOrdering();
		storeOrdering("index.order.tmi");
	}
	if (!tryLoadTopology("index.topo.tmi"))
	{
		initTopology();
		storeTopology("index.topo.tmi");
	}
	if (!tryLoadTree("index.tree.tmi"))
	{
		initTree();
		storeTree("index.tree.tmi");
	}
}

void save(std::ofstream& file, size_t value);
template <typename T, typename U> void save(std::ofstream& file, std::pair<T, U> value);
template <typename T> void save(std::ofstream& file, const std::vector<T>& value);
template <typename T, typename U> void save(std::ofstream& file, const std::unordered_map<T, U>& value);

void load(std::ifstream& file, size_t& value);
template <typename T, typename U> void load(std::ifstream& file, std::pair<T, U>& value);
template <typename T> void load(std::ifstream& file, std::vector<T>& value);
template <typename T, typename U> void load(std::ifstream& file, std::unordered_map<T, U>& value);

void save(std::ofstream& file, size_t value)
{
	file << value << " ";
}

template <typename T, typename U>
void save(std::ofstream& file, std::pair<T, U> value)
{
	save(file, value.first);
	save(file, value.second);
}

void load(std::ifstream& file, size_t& result)
{
	file >> result;
}

template <typename T, typename U>
void load(std::ifstream& file, std::pair<T, U>& result)
{
	load(file, result.first);
	load(file, result.second);
}

template <typename T>
void load(std::ifstream& file, std::vector<T>& result)
{
	size_t count = 0;
	file >> count;
	for (size_t i = 0; i < count; i++)
	{
		T item;
		load(file, item);
		result.push_back(item);
	}
}

template <typename T, typename U>
void load(std::ifstream& file, std::unordered_map<T, U>& result)
{
	size_t count = 0;
	file >> count;
	for (size_t i = 0; i < count; i++)
	{
		T key;
		U value;
		load(file, key);
		load(file, value);
		result[key] = value;
	}
}

template <typename T>
void save(std::ofstream& file, const std::vector<T>& vec)
{
	file << vec.size() << " ";
	for (size_t i = 0; i < vec.size(); i++)
	{
		save(file, vec[i]);
	}
}

template <typename T, typename U>
void save(std::ofstream& file, const std::unordered_map<T, U>& map)
{
	file << map.size() << " ";
	for (auto pair : map)
	{
		save(file, pair.first);
		save(file, pair.second);
	}
}

void MinimizerGraph::storeOrdering(std::string filename)
{
	std::ofstream file {filename};
	save(file, minmerOrdering);
}

bool MinimizerGraph::tryLoadOrdering(std::string filename)
{
	std::ifstream file {filename};
	load(file, minmerOrdering);
	return file.good();
}

void MinimizerGraph::storeTopology(std::string filename)
{
	std::ofstream file {filename};
	save(file, minmers);
	save(file, minmerIndex);
	save(file, minmerPosition);
	save(file, topology);
	save(file, reverseTopology);
}

bool MinimizerGraph::tryLoadTopology(std::string filename)
{
	std::ifstream file {filename};
	load(file, minmers);
	load(file, minmerIndex);
	load(file, minmerPosition);
	load(file, topology);
	load(file, reverseTopology);
	return file.good();
}

void MinimizerGraph::storeTree(std::string filename)
{
	std::ofstream file {filename};
	save(file, tree.nodeLabels);
	save(file, tree.children);
	save(file, tree.nodeDepths);
	save(file, tree.leafMinmers);
	save(file, tree.parents);
}

bool MinimizerGraph::tryLoadTree(std::string filename)
{
	std::ifstream file {filename};
	load(file, tree.nodeLabels);
	load(file, tree.children);
	load(file, tree.nodeDepths);
	load(file, tree.leafMinmers);
	load(file, tree.parents);
	return file.good();
}

class NodeWithPriority
{
public:
	NodeWithPriority(size_t node, std::set<size_t> actives) : node(node), actives(actives) {}
	size_t node;
	std::set<size_t> actives;
	bool operator<(const NodeWithPriority& other) const
	{
		return actives.size() < other.actives.size();
	}
	bool operator>(const NodeWithPriority& other) const
	{
		return actives.size() > other.actives.size();
	}
};

void MinimizerGraph::initTree()
{
	std::set<size_t> uniqueMinmersSet;
	for (size_t i = 0; i < minmers.size(); i++)
	{
		uniqueMinmersSet.insert(minmers[i]);
	}
	std::vector<size_t> uniqueMinmers { uniqueMinmersSet.begin(), uniqueMinmersSet.end() };
	tree.nodeLabels.push_back(0);
	tree.nodeDepths.push_back(0);
	tree.children.emplace_back();
	tree.leafMinmers.emplace_back();
	tree.parents.push_back(0);
	std::priority_queue<NodeWithPriority> nodeQueue;
	std::cerr << "build roots" << std::endl;
	for (size_t i = 0; i < uniqueMinmers.size(); i++)
	{
		size_t minmer = uniqueMinmers[i];
		size_t newIndex = tree.nodeLabels.size();
		tree.children.emplace_back();
		tree.children[0][minmer] = newIndex;
		tree.leafMinmers.emplace_back();
		tree.nodeLabels.push_back(minmer);
		tree.nodeDepths.push_back(1);
		tree.parents.push_back(0);
		std::set<size_t> actives { minmerIndex[minmer].begin(), minmerIndex[minmer].end() };
		nodeQueue.emplace(newIndex, actives);
	}
	std::cerr << "roots " << tree.nodeLabels.size() << " trunksize " << nodeQueue.top().actives.size() << std::endl;
	std::cerr << "expand tree" << std::endl;
	while (tree.nodeLabels.size() < TREE_MAX_NODE_COUNT && nodeQueue.size() > 0)
	{
		auto top = nodeQueue.top();
		if (nodeQueue.top().actives.size() < TREE_MIN_TRUNK_SIZE) break;
		nodeQueue.pop();
		for (size_t i = 0; i < uniqueMinmers.size(); i++)
		{
			size_t minmer = uniqueMinmers[i];
			std::set<size_t> activesHere;
			for (auto active : top.actives)
			{
				if (reverseTopology[active].count(minmer) == 1)
				{
					for (auto neighbor : reverseTopology[active].at(minmer))
					{
						activesHere.insert(neighbor.first);
					}
				}
			}
			if (activesHere.size() == 0) continue;
			size_t newIndex = tree.nodeLabels.size();
			tree.children.emplace_back();
			tree.children[top.node][minmer] = newIndex;
			tree.leafMinmers.emplace_back();
			tree.nodeLabels.push_back(minmer);
			tree.nodeDepths.push_back(1);
			tree.parents.push_back(top.node);
			nodeQueue.emplace(newIndex, activesHere);
			if (tree.nodeLabels.size() % 100 == 0) std::cerr << "tree " << tree.nodeLabels.size() << " trunksize " << top.actives.size() << std::endl;
		}
	}
	assert(nodeQueue.size() > 0);
	std::cerr << "treenodes: " << tree.nodeLabels.size() << std::endl;
	std::cerr << "max leaf:" << nodeQueue.top().actives.size() << std::endl;
	while (nodeQueue.size() > 0)
	{
		auto top = nodeQueue.top();
		for (auto active : top.actives)
		{
			tree.leafMinmers[top.node].emplace_back(active, tree.nodeDepths[top.node]);
		}
		nodeQueue.pop();
	}
}

void MinimizerGraph::getTotalPathsRec(size_t node, std::vector<size_t>& totalPathsPerNode) const
{
	totalPathsPerNode[node] += tree.leafMinmers[node].size();
	for (auto child : tree.children[node])
	{
		getTotalPathsRec(child.second, totalPathsPerNode);
		totalPathsPerNode[node] += totalPathsPerNode[child.second];
	}
}

void MinimizerGraph::initRandomOrdering()
{
	minmerOrdering.resize(pow(2, k * 2));
	for (size_t i = 0; i < minmerOrdering.size(); i++)
	{
		minmerOrdering[i] = i;
	}
	for (size_t i = 0; i < minmerOrdering.size(); i++)
	{
		std::swap(minmerOrdering[i], minmerOrdering[i+(rand() % (minmerOrdering.size()-i))]);
	}
}

bool MinimizerGraph::minmerCompare(size_t left, size_t right) const
{
	return minmerOrdering[left] < minmerOrdering[right];
}

bool addUnique(std::vector<size_t>& vec, size_t minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return false;
	}
	vec.emplace_back(minmer);
	return true;
}

bool addUnique(std::vector<GraphPos>& vec, GraphPos minmer)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i] == minmer) return false;
	}
	vec.emplace_back(minmer);
	return true;
}

bool addUnique(std::vector<std::pair<size_t, size_t>>& vec, size_t minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return false;
	}
	vec.emplace_back(minmer, length);
	return true;
}

bool addUnique(std::vector<std::pair<GraphPos, size_t>>& vec, GraphPos minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return false;
	}
	vec.emplace_back(minmer, length);
	return true;
}

void MinimizerGraph::addMinmerIfUnique(const std::vector<size_t>& posStack, size_t minmer, size_t pos, std::unordered_map<size_t, std::vector<size_t>>& endPositions)
{
	assert(pos >= k-1);
	for (auto index : endPositions[posStack[pos]])
	{
		if (minmers[index] == minmer) return;
	}
	endPositions[posStack[pos]].push_back(minmers.size());
	minmerIndex[minmer].push_back(minmers.size());
	minmers.push_back(minmer);
	minmerPosition.push_back(posStack[pos]);
}

void MinimizerGraph::addMinmersRec(size_t pos, std::string& currentMer, std::vector<size_t>& posStack, std::unordered_map<size_t, std::vector<size_t>>& endPositions)
{
	assert(posStack.size() == currentMer.size());
	posStack.push_back(pos);
	currentMer += graph.NodeSequences(pos);
	if (posStack.size() == w)
	{
		assert(posStack.size() == currentMer.size());
		assert(posStack.size() == w);
		auto found = findOneMinimizer(currentMer);
		std::string foundMinmer;
		for (size_t i = 0; i < k; i++)
		{
			foundMinmer += "ACGT"[(found.first >> (i*2)) % 4];
		}
		for (size_t i = 0; i < currentMer.size() - foundMinmer.size() + 1; i++)
		{
			assert((i+k-1) != found.second || currentMer.substr(i, k) == foundMinmer);
			if (currentMer.substr(i, k) == foundMinmer)
			{
				addMinmerIfUnique(posStack, found.first, i + k - 1, endPositions);
			}
		}
	}
	else
	{
		if (pos < graph.NodeEnd(graph.IndexToNode(pos))-1)
		{
			addMinmersRec(pos+1, currentMer, posStack, endPositions);
		}
		else
		{
			for (auto neighbor : graph.outNeighbors[graph.IndexToNode(pos)])
			{
				addMinmersRec(graph.NodeStart(neighbor), currentMer, posStack, endPositions);
			}
		}
	}
	currentMer.pop_back();
	posStack.pop_back();
}

void MinimizerGraph::initTopoRec(size_t minmerPos, size_t pos, size_t len, const std::unordered_map<size_t, std::vector<size_t>>& endPositions)
{
	assert(len < (k + w) * 2);
	if (endPositions.count(pos) == 1)
	{
		assert(endPositions.at(pos).size() > 0);
		for (auto previousIndex : endPositions.at(pos))
		{
			addUnique(topology[minmerPos], previousIndex, len);
		}
	}
	else
	{
		if (pos == graph.NodeStart(graph.IndexToNode(pos)))
		{
			for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
			{
				initTopoRec(minmerPos, graph.NodeEnd(neighbor)-1, len+1, endPositions);
			}
		}
		else
		{
			initTopoRec(minmerPos, pos-1, len+1, endPositions);
		}
	}
}

void MinimizerGraph::expandTopoRec(size_t start, size_t current, size_t len)
{
	for (auto inNeighbor : topology[current])
	{
		assert(inNeighbor.second > 0);
		if (len + inNeighbor.second <= MAX_GRAPH_TOPO_DISTANCE)
		{
			if (addUnique(topology[start], inNeighbor.first, len + inNeighbor.second))
			{
				expandTopoRec(start, inNeighbor.first, len + inNeighbor.second);
			}
		}
	}
}

bool MinimizerGraph::checkMinmerExistsRec(const std::string& minmer, size_t pos, size_t minmerPos) const
{
	if (graph.NodeSequences(pos) != minmer[minmerPos]) return false;
	if (minmerPos == 0) return true;
	if (pos == graph.NodeStart(graph.IndexToNode(pos)))
	{
		for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
		{
			if (checkMinmerExistsRec(minmer, graph.NodeEnd(neighbor)-1, minmerPos-1)) return true;
		}
		return false;
	}
	else
	{
		return checkMinmerExistsRec(minmer, pos-1, minmerPos-1);
	}
}

void MinimizerGraph::verifyMinmerPositions(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const
{
	size_t checked = 0;
	for (auto position : endPositions)
	{
		for (size_t index : position.second)
		{
			checked++;
			size_t minmerI = minmers[index];
			std::string minmerStr;
			for (size_t i = 0; i < k; i++)
			{
				minmerStr += "ACGT"[minmerI % 4];
				minmerI >>= 2;
			}
			assert(checkMinmerExistsRec(minmerStr, position.first, minmerStr.size()-1));
		}
	}
	assert(checked == minmers.size());
}

void MinimizerGraph::verifyDirectTopology(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const
{
	size_t checked = 0;
	for (auto position : endPositions)
	{
		for (auto index : position.second)
		{
			checked++;
			std::set<size_t> neighborsHere;
			if (position.first == graph.NodeStart(graph.IndexToNode(position.first)))
			{
				for (auto neighbor : graph.inNeighbors[graph.IndexToNode(position.first)])
				{
					checkDirectTopologyRec(endPositions, graph.NodeEnd(neighbor)-1, MAX_GRAPH_TOPO_DISTANCE-1, neighborsHere);
				}
			}
			else
			{
				checkDirectTopologyRec(endPositions, position.first-1, MAX_GRAPH_TOPO_DISTANCE-1, neighborsHere);
			}
			assert(neighborsHere.size() == topology[index].size());
			for (auto neighbor : topology[index])
			{
				assert(neighborsHere.count(neighbor.first) == 1);
			}
		}
	}
	assert(checked == minmers.size());
}

void MinimizerGraph::checkDirectTopologyRec(const std::unordered_map<size_t, std::vector<size_t>>& endPositions, size_t pos, size_t lengthLeft, std::set<size_t>& foundNeighbors) const
{
	if (endPositions.count(pos) == 1)
	{
		for (auto here : endPositions.at(pos))
		{
			foundNeighbors.insert(here);
		}
		return;
	}
	if (lengthLeft == 0) return;
	if (pos == graph.NodeStart(graph.IndexToNode(pos)))
	{
		for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
		{
			checkDirectTopologyRec(endPositions, graph.NodeEnd(neighbor)-1, lengthLeft-1, foundNeighbors);
		}
	}
	else
	{
		checkDirectTopologyRec(endPositions, pos-1, lengthLeft-1, foundNeighbors);
	}
}

void MinimizerGraph::verifyExpandedTopology(const std::unordered_map<size_t, std::vector<size_t>>& endPositions) const
{
	size_t checked = 0;
	for (auto position : endPositions)
	{
		for (auto index : position.second)
		{
			checked++;
			std::set<size_t> neighborsHere;
			if (position.first == graph.NodeStart(graph.IndexToNode(position.first)))
			{
				for (auto neighbor : graph.inNeighbors[graph.IndexToNode(position.first)])
				{
					checkExpandedTopologyRec(endPositions, graph.NodeEnd(neighbor)-1, MAX_GRAPH_TOPO_DISTANCE-1, neighborsHere);
				}
			}
			else
			{
				checkExpandedTopologyRec(endPositions, position.first-1, MAX_GRAPH_TOPO_DISTANCE-1, neighborsHere);
			}
			// assert(neighborsHere.size() == topology[index].size());
			for (auto neighbor : topology[index])
			{
				assert(neighborsHere.count(neighbor.first) == 1);
			}
		}
	}
	assert(checked == minmers.size());
}

void MinimizerGraph::checkExpandedTopologyRec(const std::unordered_map<size_t, std::vector<size_t>>& endPositions, size_t pos, size_t lengthLeft, std::set<size_t>& foundNeighbors) const
{
	if (endPositions.count(pos) == 1)
	{
		for (auto here : endPositions.at(pos))
		{
			foundNeighbors.insert(here);
		}
	}
	if (lengthLeft == 0) return;
	if (pos == graph.NodeStart(graph.IndexToNode(pos)))
	{
		for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
		{
			checkExpandedTopologyRec(endPositions, graph.NodeEnd(neighbor)-1, lengthLeft-1, foundNeighbors);
		}
	}
	else
	{
		checkExpandedTopologyRec(endPositions, pos-1, lengthLeft-1, foundNeighbors);
	}
}

void MinimizerGraph::initTopology()
{
	minmerIndex.resize(pow(2, k * 2));
	std::string kmerStack;
	std::vector<size_t> posStack;
	std::cerr << "get minimizers" << std::endl;
	std::unordered_map<size_t, std::vector<size_t>> endPositions;
	size_t processed = 0;
	for (size_t nodeIndex = 0; nodeIndex < graph.nodeStart.size(); nodeIndex++)
	{
		for (size_t i = graph.NodeStart(nodeIndex); i < graph.NodeEnd(nodeIndex); i++)
		{
			addMinmersRec(i, kmerStack, posStack, endPositions);
			assert(posStack.size() == 0);
			assert(kmerStack.size() == 0);
		}
		if (processed % 10000 == 0) std::cerr << "node " << processed << "/" << graph.nodeStart.size() << " minmers " << minmers.size() << " est " << (graph.nodeStart.size() * minmers.size()) / (processed+1) << std::endl;
		processed++;
	}
	verifyMinmerPositions(endPositions);
	topology.resize(minmers.size());
	std::cerr << "get topology" << std::endl;
	processed = 0;
	for (auto pos : endPositions)
	{
		for (auto minmerIndex : pos.second)
		{
			if (pos.first == graph.NodeStart(graph.IndexToNode(pos.first)))
			{
				for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos.first)])
				{
					initTopoRec(minmerIndex, graph.NodeEnd(neighbor)-1, 1, endPositions);
				}
			}
			else
			{
				initTopoRec(minmerIndex, pos.first-1, 1, endPositions);
			}
		}
		if (processed % 100000 == 0) std::cerr << "init topo " << processed << "/" << endPositions.size() << std::endl;
		processed++;
	}
	verifyDirectTopology(endPositions);
	processed = 0;
	for (size_t i = 0; i < minmers.size(); i++)
	{
		for (auto neighbor : topology[i])
		{
			expandTopoRec(i, neighbor.first, neighbor.second);
		}
		if (processed % 100000 == 0) std::cerr << "expand topo " << processed << "/" << minmers.size() << std::endl;
		processed++;
	}
	verifyExpandedTopology(endPositions);
	reverseTopology.resize(topology.size());
	for (size_t i = 0; i < topology.size(); i++)
	{
		for (auto neighbor : topology[i])
		{
			reverseTopology[neighbor.first][minmers[i]].emplace_back(i, neighbor.second);
		}
	}
	std::cerr << "topo finished" << std::endl;
	std::cerr << "k: " << k << " w: " << w << std::endl;
	std::cerr << "number of minimizers: " << minmers.size() << std::endl;
	size_t uniqueMinmers = 0;
	for (size_t i = 0; i < minmerIndex.size(); i++)
	{
		if (minmerIndex[i].size() > 0) uniqueMinmers++;
	}
	std::cerr << "unique minmers: " << uniqueMinmers << std::endl;
	size_t numEdges = 0;
	for (size_t i = 0; i < topology.size(); i++)
	{
		numEdges += topology[i].size();
	}
	std::cerr << "number of edges: " << numEdges << std::endl;
}

size_t MinimizerGraph::minmerize(std::string kmer) const
{
	assert(kmer.size() == k);
	size_t result = 0;
	for (size_t i = 0; i < kmer.size(); i++)
	{
		result |= charNum(kmer[i]) << (i * 2);
	}
	assert(result < pow(4, k));
	return result;
}

size_t MinimizerGraph::nextminmer(size_t minmer, char newChar) const
{
	minmer >>= 2;
	minmer |= charNum(newChar) << ((k-1) * 2);
	assert(minmer < pow(4, k));
	return minmer;
}

std::vector<std::pair<size_t, size_t>> MinimizerGraph::inNeighbors(size_t index) const
{
	return topology[index];
}

bool MinimizerGraph::findPath(const std::string& seq, std::vector<size_t>& posStack, const AlignmentGraph& graph) const
{
	if (posStack.size() == seq.size()) return true;
	if (graph.NodeSequences(posStack.back()) != seq[posStack.size()-1]) return false;
	size_t node = graph.IndexToNode(posStack.back());
	if (posStack.back() == graph.NodeEnd(node)-1)
	{
		for (auto neighbor : graph.outNeighbors[node])
		{
			posStack.push_back(graph.NodeStart(neighbor));
			if (findPath(seq, posStack, graph)) return true;
			posStack.pop_back();
		}
		return false;
	}
	else
	{
		posStack.push_back(posStack.back()+1);
		bool result = findPath(seq, posStack, graph);
		if (!result) posStack.pop_back();
		return result;
	}
}

// void MinimizerGraph::align(const std::string& originalSeq) const
// {
// 	std::vector<size_t> posInGraph;
// 	posInGraph.reserve(originalSeq.size());
// 	posInGraph.push_back(graph.NodeStart(graph.nodeLookup.at(246)[0]));
// 	assert(findPath(originalSeq, posInGraph, graph));
// 	assert(posInGraph.size() == originalSeq.size());
// 	for (size_t i = 0; i < originalSeq.size(); i++)
// 	{
// 		assert(graph.NodeSequences(posInGraph[i]) == originalSeq[i]);
// 	}
// 	doForMinimizers(originalSeq, [this, &posInGraph](size_t minmer, size_t pos)
// 	{
// 		bool found = false;
// 		for (auto index : minmerIndex[minmer])
// 		{
// 			if (minmerPosition[index] == posInGraph[pos])
// 			{
// 				found = true;
// 				break;
// 			}
// 		}
// 		assert(found);
// 	});
// }

void MinimizerGraph::align(const std::string& readname, const std::string& originalSeq) const
{
	std::cerr << "align" << std::endl;
	std::string seq = originalSeq;
	assert(seq.size() > k + w);
// #ifdef PRINTLENS
	auto timeStart = std::chrono::system_clock::now();
	std::vector<size_t> lens;
	lens.resize(seq.size(), 0);
// #endif
	std::vector<size_t> activeTreeNodes;
	std::vector<size_t> activeNextTreeNodes;
	std::vector<size_t> treeNodeLastActivePos;
	std::vector<size_t> treeNodeLastActiveLen;
	std::unordered_map<size_t, std::pair<size_t, size_t>> activePathMinmers;
	std::unordered_map<size_t, std::pair<size_t, size_t>> nextActivePathMinmers;
	activeTreeNodes.push_back(0);
	treeNodeLastActivePos.resize(tree.nodeLabels.size(), 0);
	treeNodeLastActiveLen.resize(tree.nodeLabels.size(), 0);
	doForMinimizers(seq, [this, &activePathMinmers, &nextActivePathMinmers, &activeTreeNodes, &lens, &activeNextTreeNodes, &treeNodeLastActivePos, &treeNodeLastActiveLen](size_t minmer, size_t seqPos)
	{
		assert(activeTreeNodes.size() > 0);
		assert(activeTreeNodes[0] == 0);
		treeNodeLastActivePos[0] = seqPos;
		activeNextTreeNodes.push_back(0);
		assert(nextActivePathMinmers.size() == 0);
		// std::cerr << "pos " << seqPos << " treeactives " << activeTreeNodes.size() << " minimizeractives " << activePathMinmers.size() << " minimizers " << minmerIndex[minmer].size();
		for (auto pathMinmer : activePathMinmers)
		{
			size_t index = pathMinmer.first;
			size_t length = pathMinmer.second.first;
			size_t position = pathMinmer.second.second;
			if (seqPos - position > MAX_SEQ_MINMER_DISTANCE) continue;
			if (nextActivePathMinmers[index].first < length) nextActivePathMinmers[index] = std::make_pair(length, position);
			auto found = reverseTopology[index].find(minmer);
			if (found == reverseTopology[index].end()) continue;
			for (auto neighbor : found->second)
			{
				if (nextActivePathMinmers[neighbor.first].first < length + neighbor.second) nextActivePathMinmers[neighbor.first] = std::make_pair(length + neighbor.second, seqPos);
// #ifdef PRINTLENS
				for (size_t i = position; i <= seqPos; i++)
				{
					lens[i] = std::max(lens[i], length + neighbor.second);
				}
// #endif
			}
		}
		// std::cerr << " minimizers added " << nextActivePathMinmers.size();
		std::swap(activePathMinmers, nextActivePathMinmers);
		nextActivePathMinmers.clear();
		size_t added = 0;
		for (auto node : activeTreeNodes)
		{
			assert(seqPos >= treeNodeLastActivePos[node]);
			if (seqPos - treeNodeLastActivePos[node] > MAX_SEQ_MINMER_DISTANCE) continue;
			for (auto leaf : tree.leafMinmers[node])
			{
				added++;
				if (activePathMinmers[leaf.first].first < leaf.second) activePathMinmers[leaf.first] = std::make_pair(leaf.second, seqPos);
			}
			auto found = tree.children[node].find(minmer);
			if (found == tree.children[node].end()) continue;
			activeNextTreeNodes.push_back(found->second);
			treeNodeLastActivePos[found->second] = seqPos;
			treeNodeLastActiveLen[found->second] = treeNodeLastActiveLen[node] + 1;
// // #ifdef PRINTLENS
			for (size_t i = treeNodeLastActivePos[node]; i <= seqPos; i++)
			{
				lens[i] = std::max(lens[i], treeNodeLastActiveLen[found->second]);
			}
// // #endif
		}
		// std::cerr << " tree added " << added << std::endl;
		std::swap(activeTreeNodes, activeNextTreeNodes);
		activeNextTreeNodes.clear();
	});
// #ifdef PRINTLENS
	auto timeEnd = std::chrono::system_clock::now();
	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::cerr << "read: " << readname << std::endl;
	std::cerr << "time: " << time << std::endl;
	size_t maxlen = 0;
	for (size_t i = 0; i < lens.size(); i++)
	{
		maxlen = std::max(maxlen, lens[i]);
	}
	std::cerr << "maxlen: " << maxlen << std::endl;
	// std::cerr << "lens: " << std::endl;
	// for (size_t i = 0; i < lens.size(); i++)
	// {
	// 	std::cerr << lens[i] << std::endl;
	// }
// #endif
}

// void MinimizerGraph::align(const std::string& originalSeq) const
// {
// 	std::cerr << "align" << std::endl;
// 	// std::string seq = getHPC(originalSeq);
// 	std::string seq = originalSeq;
// 	assert(seq.size() > k + w);
// // #ifdef PRINTLENS
// 	auto timeStart = std::chrono::system_clock::now();
// 	std::vector<size_t> lens;
// 	lens.resize(seq.size(), 0);
// // #endif
// 	std::vector<size_t> lastMinmerSeqPos;
// 	std::vector<size_t> lastMinmerMatchLen;
// 	lastMinmerSeqPos.resize(minmers.size(), 0);
// 	lastMinmerMatchLen.resize(minmers.size(), 0);
// 	size_t lastprint = 0;
// 	doForMinimizers(seq, [this, &seq, &lens, &lastMinmerSeqPos, &lastMinmerMatchLen](size_t minimizer, size_t seqPos)
// 	{
// 		for (auto graphPos : minmerIndex[minimizer])
// 		{
// 			auto neighbors = inNeighbors(graphPos);
// 			size_t currentLen = 0;
// // #ifdef PRINTLENS
// 			size_t lastPos = 0;
// // #endif
// 			for (size_t i = 0; i < neighbors.size(); i++)
// 			{
// 				size_t lastLen = lastMinmerMatchLen[neighbors[i].first];
// 				int seqDist = seqPos - lastMinmerSeqPos[neighbors[i].first];
// 				int topoDistance = neighbors[i].second;
// 				// assert(topoDistance <= k+w);
// 				bool valid = true;
// 				if (seqDist > SEQ_DIST_MODE_CUTOFF || topoDistance > SEQ_DIST_MODE_CUTOFF)
// 				{
// 					if (seqDist < topoDistance * 0.6 || seqDist > topoDistance / 0.6) valid = false;
// 				}
// 				// if (seqDist > MAX_SEQ_MINMER_DISTANCE) valid = false;
// 				if (valid && topoDistance + lastLen > currentLen)
// 				{
// 					currentLen = topoDistance + lastLen;
// // #ifdef PRINTLENS
// 					lastPos = lastMinmerSeqPos[neighbors[i].first];
// // #endif
// 				}
// 			}
// 			lastMinmerSeqPos[graphPos] = seqPos;
// 			lastMinmerMatchLen[graphPos] = currentLen;
// // #ifdef PRINTLENS
// 			for (size_t i = lastPos; i < seqPos; i++)
// 			{
// 				lens[i] = std::max(lens[i], currentLen);
// 			}
// // #endif
// 		}
// 		std::cerr << "align " << seqPos << "/" << seq.size() << " matches " << minmerIndex[minimizer].size() << std::endl;
// 	});
// // #ifdef PRINTLENS
// 	auto timeEnd = std::chrono::system_clock::now();
// 	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
// 	std::cerr << "time: " << time << std::endl;
// 	size_t maxlen = 0;
// 	for (size_t i = 0; i < lens.size(); i++)
// 	{
// 		maxlen = std::max(maxlen, lens[i]);
// 	}
// 	std::cerr << "maxlen: " << maxlen << std::endl;
// 	std::cerr << "lens: " << std::endl;
// 	for (size_t i = 0; i < lens.size(); i++)
// 	{
// 		std::cerr << lens[i] << std::endl;
// 	}
// // #endif
// }

std::pair<size_t, size_t> MinimizerGraph::findOneMinimizer(const std::string& seq) const
{
	size_t current = minmerize(seq.substr(0, k));
	size_t smallest = current;
	size_t pos = k-1;
	for (size_t i = k; i < seq.size(); i++)
	{
		current = nextminmer(current, seq[i]);
		if (minmerCompare(current, smallest))
		{
			smallest = current;
			pos = i;
		}
	}
	assert(minmerOrdering[smallest] != std::numeric_limits<size_t>::max());
	return std::make_pair(smallest, pos);
}

GraphPos::GraphPos(int node, size_t pos) :
node(node),
pos(pos)
{
}

bool GraphPos::operator==(const GraphPos& other) const
{
	return node == other.node && pos == other.pos;
}

bool GraphPos::operator!=(const GraphPos& other) const
{
	return !(*this == other);
}