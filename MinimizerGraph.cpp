#include <fstream>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <cmath>
#include "MinimizerGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

constexpr int MAX_GRAPH_TOPO_DISTANCE = 15;
constexpr int MAX_SEQ_MINMER_DISTANCE = 15;
// constexpr int SEQ_DIST_MODE_CUTOFF = 32;
constexpr int MAX_TREE_DEPTH = 20;
constexpr int MAX_TREE_PATH_LENGTH = 256;

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
w(w)
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
		initTopology(graph);
		storeTopology("index.topo.tmi");
	}
	if (!tryLoadTree("index.tree.tmi"))
	{
		initTree();
		storeTree("index.tree.tmi");
	}
}

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
void save(std::ofstream& file, const std::unordered_map<T, U>& map)
{
	file << map.size() << " ";
	for (auto pair : map)
	{
		save(file, pair.first);
		save(file, pair.second);
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
	save(file, topology);
}

bool MinimizerGraph::tryLoadTopology(std::string filename)
{
	std::ifstream file {filename};
	load(file, minmers);
	load(file, minmerIndex);
	load(file, topology);
	return file.good();
}

void MinimizerGraph::storeTree(std::string filename)
{
	std::ofstream file {filename};
	save(file, tree.nodeLabels);
	save(file, tree.children);
	save(file, tree.nodeDepths);
}

bool MinimizerGraph::tryLoadTree(std::string filename)
{
	std::ifstream file {filename};
	load(file, tree.nodeLabels);
	load(file, tree.children);
	load(file, tree.nodeDepths);
	return file.good();
}

void MinimizerGraph::buildTreeRec(size_t parent, const std::vector<std::pair<size_t, size_t>>& activeIndices, const std::vector<size_t>& uniqueMinmers)
{
	if (tree.nodeDepths[parent] == MAX_TREE_DEPTH) return;
	for (size_t i = 0; i < uniqueMinmers.size(); i++)
	{
		size_t minmer = uniqueMinmers[i];
		std::vector<std::pair<size_t, size_t>> nextActives;
		std::vector<size_t> leafMinmers;
		for (auto index : activeIndices)
		{
			if (reverseTopology[index.first].count(minmer) == 0) continue;
			for (auto neighbor : reverseTopology[index.first].at(minmer))
			{
				size_t length = neighbor.second + index.second;
				if (length <= MAX_TREE_PATH_LENGTH && tree.nodeDepths[parent] < MAX_TREE_DEPTH-1)
				{
					nextActives.emplace_back(neighbor.first, length);
				}
				else
				{
					leafMinmers.emplace_back(index.first);
				}
			}
		}
		if (nextActives.size() == 0 && leafMinmers.size() == 0) continue;

		size_t newIndex = tree.children.size();
		tree.children.emplace_back();
		tree.children[parent][minmer] = newIndex;
		tree.nodeLabels.push_back(minmer);
		tree.nodeDepths.push_back(tree.nodeDepths[parent] + 1);
		tree.leafMinmers.push_back(leafMinmers);

		if (nextActives.size() == 1)
		{
			tree.leafMinmers.back().push_back(nextActives[0].first);
			continue;
		}
		if (nextActives.size() == 0) continue;
		buildTreeRec(newIndex, nextActives, uniqueMinmers);
	}
}

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
	for (size_t i = 0; i < uniqueMinmers.size(); i++)
	{
		std::cerr << "branch " << i << "/" << uniqueMinmers.size() << " nodes " << tree.nodeLabels.size() << std::endl;
		size_t minmer = uniqueMinmers[i];
		size_t newIndex = tree.nodeLabels.size();
		tree.children.emplace_back();
		tree.children[0][minmer] = newIndex;
		tree.leafMinmers.emplace_back();
		tree.nodeLabels.push_back(minmer);
		tree.nodeDepths.push_back(1);
		std::vector<std::pair<size_t, size_t>> activeIndices;
		for (auto index : minmerIndex[minmer])
		{
			activeIndices.emplace_back(index, 0);
		}
		buildTreeRec(newIndex, activeIndices, uniqueMinmers);
	}
	std::cerr << "tree size " << tree.nodeLabels.size() << std::endl;
	std::vector<size_t> totalPathsPerNode;
	totalPathsPerNode.resize(tree.nodeLabels.size(), 0);
	getTotalPathsRec(0, totalPathsPerNode);
	std::vector<size_t> numNodesPerDepth;
	std::vector<size_t> numPathsPerDepth;
	std::vector<size_t> maxPathsPerDepth;
	numNodesPerDepth.resize(MAX_TREE_DEPTH+1, 0);
	numPathsPerDepth.resize(MAX_TREE_DEPTH+1, 0);
	maxPathsPerDepth.resize(MAX_TREE_DEPTH+1, 0);
	for (size_t i = 0; i < tree.nodeLabels.size(); i++)
	{
		numNodesPerDepth[tree.nodeDepths[i]]++;
		numPathsPerDepth[tree.nodeDepths[i]] += totalPathsPerNode[i];
		maxPathsPerDepth[tree.nodeDepths[i]] = std::max(maxPathsPerDepth[tree.nodeDepths[i]], totalPathsPerNode[i]);
	}
	for (size_t i = 0; i < MAX_TREE_DEPTH+1; i++)
	{
		std::cerr << "depth " << i << " nodes " << numNodesPerDepth[i] << " paths " << numPathsPerDepth[i] << " avg " << (double)numPathsPerDepth[i] / (double)numNodesPerDepth[i] << " max " << maxPathsPerDepth[i] << std::endl;
	}
}

void MinimizerGraph::getTotalPathsRec(size_t node, std::vector<size_t>& totalPathsPerNode)
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
}

void MinimizerGraph::addMinmersRec(size_t pos, std::string& currentMer, std::vector<size_t>& posStack, std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph)
{
	assert(posStack.size() == currentMer.size());
	posStack.push_back(pos);
	currentMer += graph.NodeSequences(pos);
	if (posStack.size() == w)
	{
		assert(posStack.size() == currentMer.size());
		assert(posStack.size() == w);
		auto found = findOneMinimizer(currentMer);
		addMinmerIfUnique(posStack, found.first, found.second, endPositions);
	}
	else
	{
		if (pos < graph.NodeEnd(graph.IndexToNode(pos))-1)
		{
			addMinmersRec(pos+1, currentMer, posStack, endPositions, graph);
		}
		else
		{
			for (auto neighbor : graph.outNeighbors[graph.IndexToNode(pos)])
			{
				addMinmersRec(graph.NodeStart(neighbor), currentMer, posStack, endPositions, graph);
			}
		}
	}
	currentMer.pop_back();
	posStack.pop_back();
}

void MinimizerGraph::initTopoRec(size_t minmerIndex, size_t pos, size_t len, const std::unordered_map<size_t, std::vector<size_t>>& endPositions, const AlignmentGraph& graph)
{
	assert(len < (k + w) * 2);
	if (endPositions.count(pos) == 1)
	{
		assert(endPositions.at(pos).size() > 0);
		for (auto previousIndex : endPositions.at(pos))
		{
			addUnique(topology[minmerIndex], previousIndex, len);
		}
	}
	else
	{
		if (pos == graph.NodeStart(graph.IndexToNode(pos)))
		{
			for (auto neighbor : graph.inNeighbors[graph.IndexToNode(pos)])
			{
				initTopoRec(minmerIndex, graph.NodeEnd(neighbor)-1, len+1, endPositions, graph);
			}
		}
		else
		{
			initTopoRec(minmerIndex, pos-1, len+1, endPositions, graph);
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

void MinimizerGraph::initTopology(const AlignmentGraph& graph)
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
			addMinmersRec(i, kmerStack, posStack, endPositions, graph);
			assert(posStack.size() == 0);
			assert(kmerStack.size() == 0);
		}
		if (processed % 10000 == 0) std::cerr << "node " << processed << "/" << graph.nodeStart.size() << " minmers " << minmers.size() << " est " << (graph.nodeStart.size() * minmers.size()) / (processed+1) << std::endl;
		processed++;
	}
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
					initTopoRec(minmerIndex, graph.NodeEnd(neighbor)-1, 1, endPositions, graph);
				}
			}
			else
			{
				initTopoRec(minmerIndex, pos.first-1, 1, endPositions, graph);
			}
		}
		if (processed % 100000 == 0) std::cerr << "init topo " << processed << "/" << endPositions.size() << std::endl;
		processed++;
	}
	processed = 0;
	for (size_t i = 0; i < minmers.size(); i++)
	{
		expandTopoRec(i, i, 0);
		if (processed % 100000 == 0) std::cerr << "expand topo " << processed << "/" << minmers.size() << std::endl;
		processed++;
	}
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

void MinimizerGraph::align(const std::string& originalSeq) const
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
	activeTreeNodes.push_back(0);
	treeNodeLastActivePos.resize(tree.nodeLabels.size(), 0);
	treeNodeLastActiveLen.resize(tree.nodeLabels.size(), 0);
	doForMinimizers(seq, [this, &activeTreeNodes, &lens, &activeNextTreeNodes, &treeNodeLastActivePos, &treeNodeLastActiveLen](size_t minmer, size_t seqPos)
	{
		assert(activeTreeNodes.size() > 0);
		assert(activeTreeNodes[0] == 0);
		treeNodeLastActivePos[0] = seqPos;
		activeNextTreeNodes.push_back(0);
		// std::cerr << "pos " << seqPos << " actives " << activeTreeNodes.size() << " minimizers " << minmerIndex[minmer].size();
		// if (tree.children[0].count(minmer) == 0) std::cerr << " NOT IN THE ROOT!";
		// std::cerr << std::endl;
		for (auto node : activeTreeNodes)
		{
			assert(seqPos >= treeNodeLastActivePos[node]);
			if (seqPos - treeNodeLastActivePos[node] > MAX_SEQ_MINMER_DISTANCE) continue;
			auto found = tree.children[node].find(minmer);
			if (found == tree.children[node].end()) continue;
			activeNextTreeNodes.push_back(found->second);
			treeNodeLastActivePos[found->second] = seqPos;
			treeNodeLastActiveLen[found->second] = treeNodeLastActiveLen[node] + 1;
// #ifdef PRINTLENS
			for (size_t i = treeNodeLastActivePos[node]; i <= seqPos; i++)
			{
				lens[i] = std::max(lens[i], treeNodeLastActiveLen[found->second]);
			}
// #endif
		}
		std::swap(activeTreeNodes, activeNextTreeNodes);
		activeNextTreeNodes.clear();
	});
// #ifdef PRINTLENS
	auto timeEnd = std::chrono::system_clock::now();
	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::cerr << "time: " << time << std::endl;
	size_t maxlen = 0;
	for (size_t i = 0; i < lens.size(); i++)
	{
		maxlen = std::max(maxlen, lens[i]);
	}
	std::cerr << "maxlen: " << maxlen << std::endl;
	std::cerr << "lens: " << std::endl;
	for (size_t i = 0; i < lens.size(); i++)
	{
		std::cerr << lens[i] << std::endl;
	}
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