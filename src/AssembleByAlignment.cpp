#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "GfaGraph.h"

struct Path
{
	std::vector<NodePos> position;
	std::vector<size_t> nodeSize;
};

struct Alignment
{
	size_t leftPath;
	size_t rightPath;
	std::vector<std::pair<size_t, size_t>> alignedPairs;
	size_t alignmentLength;
	double alignmentIdentity;
};

struct TransitiveClosureMapping
{
	std::map<std::pair<size_t, size_t>, size_t> mapping;
};

std::vector<Path> loadAlignmentsAsPaths(std::string fileName)
{
	std::unordered_map<std::string, std::vector<vg::Alignment>> alnsPerRead;
	{
		auto vgAlns = CommonUtils::LoadVGAlignments(fileName);
		for (auto aln : vgAlns)
		{
			alnsPerRead[aln.name()].push_back(aln);
		}
	}
	std::vector<Path> result;
	for (auto pair : alnsPerRead)
	{
		assert(pair.second.size() > 0);
		std::sort(pair.second.begin(), pair.second.end(), [](const vg::Alignment& left, const vg::Alignment& right) { return left.query_position() < right.query_position(); });
		result.emplace_back();
		for (auto aln : pair.second)
		{
			for (int i = 0; i < aln.path().mapping_size(); i++)
			{
				result.back().position.emplace_back(aln.path().mapping(i).position().node_id(), !aln.path().mapping(i).position().is_reverse());
			}
		}
		if (result.back().position.size() == 0) result.pop_back();
	}
	std::cerr << result.size() << " paths" << std::endl;
	return result;
}

std::vector<Alignment> align(const std::vector<NodePos>& leftPath, const std::vector<NodePos>& rightPath, const std::vector<size_t>& leftNodeSize, const std::vector<size_t>& rightNodeSize, size_t left, size_t right, double mismatchPenalty)
{
	enum BacktraceType
	{
		Insertion,
		Deletion,
		Match,
		Mismatch,
		Start
	};
	static std::vector<std::vector<double>> DPscores;
	static std::vector<std::vector<BacktraceType>> DPtrace;
	if (DPscores.size() < leftPath.size()+1) DPscores.resize(leftPath.size()+1);
	if (DPtrace.size() < leftPath.size()+1) DPtrace.resize(leftPath.size()+1);
	if (DPscores.back().size() < rightPath.size()+1)
	{
		for (size_t i = 0; i < DPscores.size(); i++)
		{
			DPscores[i].resize(rightPath.size()+1, 0);
			DPtrace[i].resize(rightPath.size()+1, Start);
		}
	}
	std::set<std::pair<size_t, size_t>> maxima;
	for (size_t i = 0; i < leftPath.size(); i++)
	{
		for (size_t j = 0; j < rightPath.size(); j++)
		{
			DPscores[i+1][j+1] = 0;
			DPtrace[i+1][j+1] = Start;
			bool match = (leftPath[i] == rightPath[j]);
			size_t leftSize = leftNodeSize[i];
			size_t rightSize = rightNodeSize[j];
			double insertionCost = leftSize * mismatchPenalty;
			double deletionCost = rightSize * mismatchPenalty;
			double mismatchCost = std::max(insertionCost, deletionCost);
			double matchScore = leftSize;
			assert(!match || leftSize == rightSize);
			if (DPscores[i][j+1] - insertionCost > DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i][j+1] - insertionCost;
				DPtrace[i+1][j+1] = Insertion;
			}
			if (DPscores[i+1][j] - deletionCost > DPscores[i+1][j+1])
			{
				DPscores[i+1][j+1] = DPscores[i+1][j] - deletionCost;
				DPtrace[i+1][j+1] = Deletion;
			}
			if (match && DPscores[i][j] + matchScore > DPscores[i][j])
			{
				DPscores[i+1][j+1] = DPscores[i][j] + matchScore;
				DPtrace[i+1][j+1] = Match;
				maxima.erase(std::make_pair(i, j));
				maxima.emplace(i+1, j+1);
			}
			if (!match && DPscores[i][j] + mismatchCost > DPscores[i][j])
			{
				DPscores[i+1][j+1] = DPscores[i][j] + mismatchCost;
				DPtrace[i+1][j+1] = Mismatch;
			}
		}
	}
	size_t matchLen = 0;
	size_t mismatchLen = 0;
	std::vector<Alignment> result;
	for (auto maximum : maxima)
	{
		size_t maxI = maximum.first;
		size_t maxJ = maximum.second;
		result.emplace_back();
		result.back().leftPath = left;
		result.back().rightPath = right;
		result.back().alignmentLength = 0;
		while (DPtrace[maxI][maxJ] != Start)
		{
			assert(maxI > 0);
			assert(maxJ > 0);
			size_t leftSize = leftNodeSize[maxI-1];
			size_t rightSize = rightNodeSize[maxJ-1];
			switch(DPtrace[maxI][maxJ])
			{
				case Insertion:
					mismatchLen += leftSize;
					maxI -= 1;
					continue;
				case Deletion:
					mismatchLen += rightSize;
					maxJ -= 1;
					continue;
				case Match:
					assert(leftSize == rightSize);
					result.back().alignedPairs.emplace_back(maxI-1, maxJ-1);
					matchLen += leftSize;
					maxI -= 1;
					maxJ -= 1;
					continue;
				case Mismatch:
					mismatchLen += std::max(leftSize, rightSize);
					maxI -= 1;
					maxJ -= 1;
					continue;
				case Start:
				default:
					assert(false);
			}
		}
		result.back().alignmentLength = matchLen + mismatchLen;
		if (result.back().alignmentLength == 0)
		{
			result.back().alignmentIdentity = 0;
		}
		else
		{
			result.back().alignmentIdentity = (double)matchLen / ((double)matchLen + (double)mismatchLen);
		}
	}
	return result;
}

std::vector<Alignment> induceOverlaps(const std::vector<Path>& paths, double mismatchPenalty, size_t minAlnLength, double minAlnIdentity)
{
	std::vector<Alignment> result;
	std::unordered_map<NodePos, std::vector<size_t>> crossesNode;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (auto node : paths[i].position)
		{
			crossesNode[node].push_back(i);
		}
	}
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << i << "/" << paths.size() << std::endl;
		std::unordered_map<size_t, size_t> possibleMatches;
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			auto node = paths[i].position[j];
			size_t nodeSize = paths[i].nodeSize[j];
			for (auto other : crossesNode[node])
			{
				if (other <= i) continue;
				possibleMatches[other] += nodeSize;
			}
		}
		for (auto pair : possibleMatches)
		{
			size_t j = pair.first;
			if (pair.second < minAlnLength) continue;
			if (i == j) continue;
			auto alns = align(paths[i].position, paths[j].position, paths[i].nodeSize, paths[j].nodeSize, i, j, mismatchPenalty);
			for (auto aln : alns)
			{
				if (aln.alignmentLength < minAlnLength) continue;
				if (aln.alignmentIdentity < minAlnIdentity) continue;
				result.push_back(aln);
			}
		}
	}
	std::cerr << result.size() << " induced alignments" << std::endl;
	return result;
}

template <typename T>
T find(std::map<T, T>& parent, T key)
{
	if (parent.count(key) == 0)
	{
		parent[key] = key;
		return key;
	}
	if (parent.at(key) == key)
	{
		return key;
	}
	auto result = find(parent, parent.at(key));
	parent[key] = result;
	return result;
}

template <typename T>
void set(std::map<T, T>& parent, T key, T target)
{
	auto found = find(parent, key);
	parent[found] = find(parent, target);
}

TransitiveClosureMapping getTransitiveClosures(const std::vector<Alignment>& alns)
{
	TransitiveClosureMapping result;
	std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> parent;
	for (auto aln : alns)
	{
		for (auto pair : aln.alignedPairs)
		{
			if (pair.first == -1 || pair.second == -1) continue;
			std::pair<size_t, size_t> leftKey { aln.leftPath, pair.first };
			std::pair<size_t, size_t> rightKey { aln.rightPath, pair.second };
			set(parent, leftKey, rightKey);
		}
	}
	std::map<std::pair<size_t, size_t>, size_t> closureNumber;
	size_t nextClosure = 0;
	for (auto key : parent)
	{
		auto found = find(parent, key.first);
		if (closureNumber.count(found) == 0)
		{
			closureNumber[found] = nextClosure;
			nextClosure += 1;
		}
		result.mapping[key.first] = closureNumber.at(found);
	}
	std::cerr << nextClosure << " transitive closure sets" << std::endl;
	std::cerr << result.mapping.size() << " transitive closure items" << std::endl;
	return result;
}

GfaGraph getGraph(const TransitiveClosureMapping& transitiveClosures, const std::vector<Path>& paths, const GfaGraph& graph)
{
	std::unordered_set<size_t> outputtedClosures;
	GfaGraph result;
	result.edgeOverlap = graph.edgeOverlap;
	for (auto pair : transitiveClosures.mapping)
	{
		if (outputtedClosures.count(pair.second) == 1) continue;
		NodePos pos = paths[pair.first.first].position[pair.first.second];
		auto seq = graph.nodes.at(pos.id);
		if (!pos.end) seq = CommonUtils::ReverseComplement(seq);
		result.nodes[pair.second] = seq;
		outputtedClosures.insert(pair.second);
	}
	std::cerr << outputtedClosures.size() << " outputted closures" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 1; j < paths[i].position.size(); j++)
		{
			std::pair<size_t, size_t> previousKey { i, j-1 };
			std::pair<size_t, size_t> currentKey { i, j };
			if (transitiveClosures.mapping.count(previousKey) == 0 || transitiveClosures.mapping.count(currentKey) == 0) continue;
			assert(transitiveClosures.mapping.count(previousKey) == 1);
			assert(transitiveClosures.mapping.count(currentKey) == 1);
			assert(result.nodes.count(transitiveClosures.mapping.at(previousKey)) == 1);
			assert(result.nodes.count(transitiveClosures.mapping.at(currentKey)) == 1);
			NodePos from { transitiveClosures.mapping.at(previousKey), true };
			NodePos to { transitiveClosures.mapping.at(currentKey), true };
			result.edges[from].push_back(to);
		}
	}
	return result;
}

TransitiveClosureMapping insertMiddles(const TransitiveClosureMapping& original, const std::vector<Path>& paths)
{
	TransitiveClosureMapping result;
	result = original;
	size_t nextNum = 0;
	for (auto pair : original.mapping)
	{
		nextNum = std::max(nextNum, pair.second);
	}
	nextNum =+ 1;
	for (size_t i = 0; i < paths.size(); i++)
	{
		size_t firstExisting = paths[i].position.size();
		size_t lastExisting = paths[i].position.size();
		for (size_t j = 0; j < paths[i].position.size(); j++)
		{
			std::pair<size_t, size_t> key { i, j };
			if (original.mapping.count(key) == 1)
			{
				if (firstExisting == paths[i].position.size()) firstExisting = j;
				lastExisting = j;
			}
		}
		for (size_t j = firstExisting; j < lastExisting; j++)
		{
			std::pair<size_t, size_t> key { i, j };
			if (original.mapping.count(key) == 1) continue;
			result.mapping[key] = nextNum;
			nextNum += 1;
		}
	}
	std::cerr << nextNum << " transitive closure sets after inserting middles" << std::endl;
	std::cerr << result.mapping.size() << " transitive closure items after inserting middles" << std::endl;
	return result;
}

std::vector<Path> addNodeLengths(const std::vector<Path>& original, const GfaGraph& graph)
{
	std::vector<Path> result = original;
	for (size_t i = 0; i < original.size(); i++)
	{
		result[i].nodeSize.reserve(result[i].position.size());
		for (size_t j = 0; j < original[i].position.size(); j++)
		{
			result[i].nodeSize.push_back(graph.nodes.at(original[i].position[j].id).size() - graph.edgeOverlap);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string inputAlns { argv[2] };
	size_t alignmentEndCut = std::stol(argv[3]);
	double mismatchPenalty = std::stod(argv[4]);
	size_t minAlnLength = std::stol(argv[5]);
	double minAlnIdentity = std::stod(argv[6]);
	std::string outputGraph { argv[7] };

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	std::cerr << "load alignments" << std::endl;
	auto paths = loadAlignmentsAsPaths(inputAlns);
	std::cerr << "add node lengths" << std::endl;
	paths = addNodeLengths(paths, graph);
	std::cerr << "induce overlaps" << std::endl;
	auto alns = induceOverlaps(paths, mismatchPenalty, minAlnLength, minAlnIdentity);
	std::cerr << "get transitive closure" << std::endl;
	auto transitiveClosures = getTransitiveClosures(alns);
	std::cerr << "insert middles" << std::endl;
	transitiveClosures = insertMiddles(transitiveClosures, paths);
	std::cerr << "graphify" << std::endl;
	auto result = getGraph(transitiveClosures, paths, graph);
	std::cerr << "output" << std::endl;
	result.SaveToFile(outputGraph);
}