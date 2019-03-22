#include <mutex>
#include <thread>
#include <iostream>
#include <regex>
#include <fstream>
#include <string>
#include "GraphAlignerBitvectorBanded.h"
#include "BigraphToDigraph.h"
#include "vg.pb.h"
#include "AlignmentGraph.h"
#include "Aligner.h"
#include "GfaGraph.h"
#include "CommonUtils.h"
#include "fastqloader.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "MummerSeeder.h"

struct FusionAlignment
{
	FusionAlignment() {}
	std::string readName;
	double scoreFraction;
	std::string leftExon;
	bool leftReverse;
	std::string rightExon;
	bool rightReverse;
	size_t leftLen;
	size_t rightLen;
	std::string leftGene;
	std::string rightGene;
	int scoreDifference;
	std::string corrected;
};

struct BestPrefixAlignment
{
	std::vector<int> scores;
	std::vector<size_t> geneIndex;
};

struct PrefixSuffixAlignment
{
	std::string leftGene;
	std::string rightGene;
	size_t breakpoint;
};

std::string geneFromTranscript(std::string transcript)
{
	std::regex generegex("[_ ]gene:(ENSG\\d{11}\\.\\d{1,2})[_ ]");
	std::smatch match;
	std::regex_search(transcript, match, generegex);
	assert(match.ready());
	assert(!match.empty());
	assert(match.size() >= 2);
	assert(match[1].matched);
	return std::string { match[1].first, match[1].second };
}

std::unordered_map<std::string, std::unordered_set<int>> getGeneBelongers(const std::vector<vg::Alignment>& alns, const GfaGraph& graph)
{
	std::unordered_map<std::string, std::unordered_set<int>> result;
	for (auto aln : alns)
	{
		std::string gene = geneFromTranscript(aln.name());
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			assert(graph.nodes.count(aln.path().mapping(i).position().node_id()) == 1);
			result[gene].insert(aln.path().mapping(i).position().node_id());
		}
	}
	return result;
}

GfaGraph getNonfusionGraph(std::string gene, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers)
{
	assert(geneBelongers.count(gene) == 1);
	return graph.GetSubgraph(geneBelongers.at(gene));
}

void addBestPrefixAlns(BestPrefixAlignment& bestAlns, size_t geneIndex, const GraphAlignerBitvectorBanded<size_t, int32_t, uint64_t>& tableGetter, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, const FastQ& read)
{
	auto table = tableGetter.getTableFromFullStart(read.sequence[0], read.sequence.substr(1), true, reusableState);
	if (bestAlns.scores.size() == 0)
	{
		bestAlns.scores.resize(read.sequence.size()-1, std::numeric_limits<int>::max() - read.sequence.size());
		bestAlns.geneIndex.resize(read.sequence.size()-1, 0);
	}
	assert(bestAlns.scores.size() == read.sequence.size()-1);
	assert(bestAlns.geneIndex.size() == read.sequence.size()-1);
	for (size_t i = 1; i < table.slices.size(); i++)
	{
		for (size_t o = 0; o < 64; o++)
		{
			size_t j = table.slices[i].j + o;
			if (j >= bestAlns.scores.size()) break;
			auto score = table.slices[i].minPerRow.getValue(o);
			if (score < bestAlns.scores[j])
			{
				bestAlns.scores[j] = score;
				bestAlns.geneIndex[j] = geneIndex;
			}
		}
	}
}

void mergeBestPrefixAlns(BestPrefixAlignment& to, const BestPrefixAlignment& from)
{
	assert(to.scores.size() == from.scores.size());
	assert(to.geneIndex.size() == from.geneIndex.size());
	assert(to.scores.size() == to.geneIndex.size());
	for (size_t i = 0; i < to.scores.size(); i++)
	{
		if (from.scores[i] < to.scores[i])
		{
			to.scores[i] = from.scores[i];
			to.geneIndex[i] = from.geneIndex[i];
		}
	}
}

std::unordered_map<std::string, BestPrefixAlignment> getBestPrefixAlignments(const std::unordered_map<std::string, std::unordered_set<size_t>>& hasSeeds, const std::vector<std::string>& checkables, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers, const std::vector<FastQ>& reads, size_t numThreads)
{
	std::cerr << "get best prefix alignments" << std::endl;
	std::vector<std::thread> threads;
	std::vector<std::unordered_map<size_t, BestPrefixAlignment>> resultPerThread;
	resultPerThread.resize(numThreads);
	size_t nextPair = 0;
	std::mutex nextPairMutex;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([thread, &hasSeeds, &nextPair, &nextPairMutex, &checkables, &resultPerThread, &graph, &geneBelongers, &reads]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextPair == checkables.size()) break;
					i = nextPair;
					std::cerr << "prefix gene " << nextPair << "/" << checkables.size() << std::endl;
					nextPair += 1;
				}
				auto nonfusiongraph = getNonfusionGraph(checkables[i], graph, geneBelongers);
				auto alignmentGraph = DirectedGraph::BuildFromGFA(nonfusiongraph, true);
				GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, 1000, false };
				GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params { 1000, 0, alignmentGraph, std::numeric_limits<size_t>::max(), true, false, false, true, true };
				GraphAlignerBitvectorBanded<size_t, int32_t, uint64_t> tableGetter { params };
				for (auto index : hasSeeds.at(checkables[i]))
				{
					addBestPrefixAlns(resultPerThread[thread][index], i, tableGetter, reusableState, reads[index]);
				}
			}
		});
	}
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads[thread].join();
	}
	threads.clear();
	std::unordered_map<std::string, BestPrefixAlignment> result;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		for (auto pair : resultPerThread[thread])
		{
			if (result.count(reads[pair.first].seq_id) == 0)
			{
				result[reads[pair.first].seq_id] = pair.second;
			}
			else
			{
				mergeBestPrefixAlns(result[reads[pair.first].seq_id], pair.second);
			}
		}
	}
	return result;
}

void writeFusions(const std::vector<FusionAlignment>& result, std::string filename)
{
	std::ofstream file { filename };
	for (auto aln : result)
	{
		file << aln.readName << "\t" << aln.scoreFraction << "\t" << aln.scoreDifference << "\t" << aln.leftGene << "\t" << aln.rightGene << "\t" << aln.leftLen << "\t" << aln.leftExon << "\t" << (aln.leftReverse ? "-" : "+") << "\t" << aln.rightExon << "\t" << (aln.rightReverse ? "-" : "+") << "\t" << aln.rightLen << std::endl;
	}
}

void writeCorrected(const std::vector<FusionAlignment>& result, const GfaGraph& graph, std::string filename)
{
	std::ofstream file { filename };
	for (auto aln : result)
	{
		file << ">" << aln.readName << std::endl;
		file << aln.corrected << std::endl;
	}
}

std::unordered_map<std::string, std::unordered_set<size_t>> getExtraGeneMatches(const std::vector<vg::Alignment>& transcripts, const std::vector<FastQ>& reads, size_t numThreads, int seedSize)
{
	GfaGraph fakeGraph;
	std::vector<std::string> nameMapping;
	for (auto transcript : transcripts)
	{
		fakeGraph.nodes[nameMapping.size()] = transcript.sequence();
		nameMapping.push_back(geneFromTranscript(transcript.name()));
	}
	auto seeder = MummerSeeder(fakeGraph, "");
	std::vector<std::unordered_map<std::string, std::unordered_set<size_t>>> resultPerThread;
	resultPerThread.resize(numThreads);
	size_t nextIndex = 0;
	std::mutex nextPairMutex;
	std::vector<std::thread> threads;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&nameMapping, &reads, &nextIndex, &nextPairMutex, &resultPerThread, thread, &seeder, seedSize]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextIndex == reads.size()) break;
					i = nextIndex;
					std::cerr << "seed " << nextIndex << "/" << reads.size() << std::endl;
					nextIndex += 1;
				}
				auto seeds = seeder.getMemSeeds(reads[i].sequence, -1, seedSize);
				for (auto seed : seeds)
				{
					resultPerThread[thread][nameMapping[seed.nodeID]].insert(i);
				}
			}
		});
	}
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads[thread].join();
	}
	threads.clear();
	std::unordered_map<std::string, std::unordered_set<size_t>> result;
	for (size_t i = 0; i < resultPerThread.size(); i++)
	{
		for (auto pair : resultPerThread[i])
		{
			result[pair.first].insert(pair.second.begin(), pair.second.end());
		}
	}
	return result;
}

std::vector<FastQ> reverseComplements(const std::vector<FastQ>& reads)
{
	std::vector<FastQ> result;
	result.reserve(reads.size());
	for (size_t i = 0; i < reads.size(); i++)
	{
		result.push_back(reads[i].reverseComplement());
	}
	return result;
}

std::vector<std::pair<std::string, PrefixSuffixAlignment>> mergePrefixSuffix(const std::unordered_map<std::string, BestPrefixAlignment>& prefixAlns, const std::unordered_map<std::string, BestPrefixAlignment>& suffixAlns, const std::vector<std::string>& geneOrder)
{
	std::vector<std::pair<std::string, PrefixSuffixAlignment>> result;
	assert(prefixAlns.size() == suffixAlns.size());
	for (auto pair : prefixAlns)
	{
		auto prefix = pair.second;
		auto suffix = suffixAlns.at(pair.first);
		size_t bestIndex = 0;
		int bestScore = prefix.scores[prefix.scores.size()-2] + suffix.scores[suffix.scores.size()-2];
		//the scores can be different because of banding
		// assert(prefix.scores.back() - suffix.scores.back() >= -1 && prefix.scores.back() - suffix.scores.back() <= 1); // todo fix off by one error somewhere
		assert(prefix.scores.size() == suffix.scores.size());
		for (size_t i = 2; i < prefix.scores.size()-2; i++)
		{
			int scoreHere = prefix.scores[i] + suffix.scores[suffix.scores.size()-i];
			assert(scoreHere <= prefix.scores.back() + suffix.scores.back());
			if (scoreHere < bestScore)
			{
				bestScore = scoreHere;
				bestIndex = i;
			}
		}
		assert(bestIndex != 0);
		result.emplace_back();
		result.back().first = pair.first;
		result.back().second.leftGene = geneOrder[prefix.geneIndex[bestIndex]];
		result.back().second.rightGene = geneOrder[suffix.geneIndex[suffix.scores.size()-bestIndex]];
		result.back().second.breakpoint = bestIndex;
	}
	return result;
}

std::vector<std::string> getGeneOrder(const std::unordered_map<std::string, std::unordered_set<size_t>>& hasSeeds)
{
	std::vector<std::string> result;
	for (auto pair : hasSeeds)
	{
		result.push_back(pair.first);
	}
	return result;
}

std::string getCorrected(const GfaGraph& graph, const vg::Alignment& leftAln, const vg::Alignment& rightAln)
{
	std::string result;
	for (int i = 0; i < leftAln.path().mapping_size(); i++)
	{
		for (int j = 0; j < leftAln.path().mapping(i).edit_size(); j++)
		{
			std::string n = graph.nodes.at(leftAln.path().mapping(i).position().node_id());
			if (leftAln.path().mapping(i).position().is_reverse()) n = CommonUtils::ReverseComplement(n);
			result += n.substr(leftAln.path().mapping(i).position().offset(), leftAln.path().mapping(i).edit(j).from_length());
		}
	}
	result += 'N';
	for (int i = 0; i < rightAln.path().mapping_size(); i++)
	{
		for (int j = 0; j < rightAln.path().mapping(i).edit_size(); j++)
		{
			std::string n = graph.nodes.at(rightAln.path().mapping(i).position().node_id());
			if (rightAln.path().mapping(i).position().is_reverse()) n = CommonUtils::ReverseComplement(n);
			result += n.substr(rightAln.path().mapping(i).position().offset(), rightAln.path().mapping(i).edit(j).from_length());
		}
	}
	return result;
}

std::pair<std::string, bool> getLastExon(const GfaGraph& graph, const vg::Alignment& aln)
{
	if (aln.path().mapping_size() == 1)
	{
		return std::make_pair(graph.originalNodeName.at(aln.path().mapping(0).position().node_id()), aln.path().mapping(0).position().is_reverse());
	}
	return std::make_pair(graph.originalNodeName.at(aln.path().mapping(aln.path().mapping_size()-2).position().node_id()), aln.path().mapping(aln.path().mapping_size()-2).position().is_reverse());
}

std::pair<std::string, bool> getFirstExon(const GfaGraph& graph, const vg::Alignment& aln)
{
	if (aln.path().mapping_size() == 1)
	{
		return std::make_pair(graph.originalNodeName.at(aln.path().mapping(0).position().node_id()), aln.path().mapping(0).position().is_reverse());
	}
	return std::make_pair(graph.originalNodeName.at(aln.path().mapping(1).position().node_id()), aln.path().mapping(1).position().is_reverse());
}

std::vector<FusionAlignment> getAlnFromMerge(const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers, const std::vector<std::pair<std::string, PrefixSuffixAlignment>>& pairs, const std::unordered_map<std::string, BestPrefixAlignment>& prefixAlns, const std::vector<FastQ>& reads, size_t numThreads)
{
	std::unordered_map<std::string, FastQ> readmap;
	for (auto read : reads)
	{
		readmap[read.seq_id] = read;
	}
	std::vector<std::vector<FusionAlignment>> resultPerThread;
	resultPerThread.resize(numThreads);
	std::vector<std::thread> threads;
	size_t nextIndex = 0;
	std::mutex nextPairMutex;
	for (int thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&resultPerThread, thread, &nextIndex, &nextPairMutex, &pairs, &readmap, &graph, &geneBelongers, &reads, &prefixAlns](){
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextIndex == pairs.size()) break;
					i = nextIndex;
					std::cerr << "fusion " << nextIndex << "/" << pairs.size() << std::endl;
					nextIndex += 1;
				}
				auto pair = pairs[i];

				auto read = readmap.at(pair.first);
				auto leftGraph = getNonfusionGraph(pair.second.leftGene, graph, geneBelongers);
				auto leftAlnGraph = DirectedGraph::BuildFromGFA(leftGraph, true);
				GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState leftReusableState { leftAlnGraph, 1000, false };
				auto leftAlns = AlignOneWay(leftAlnGraph, read.seq_id, read.sequence.substr(0, pair.second.breakpoint), 1000, 1000, true, leftReusableState, false, true, false);
				auto leftAln = *leftAlns.alignments[0].alignment;
				replaceDigraphNodeIdsWithOriginalNodeIds(leftAln, leftAlnGraph);

				auto rightGraph = getNonfusionGraph(pair.second.rightGene, graph, geneBelongers);
				auto rightAlnGraph = DirectedGraph::BuildFromGFA(rightGraph, true);
				GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState rightReusableState { rightAlnGraph, 1000, false };
				auto rightAlns = AlignOneWay(rightAlnGraph, read.seq_id, read.sequence.substr(pair.second.breakpoint), 1000, 1000, true, rightReusableState, false, true, false);
				auto rightAln = *rightAlns.alignments[0].alignment;
				replaceDigraphNodeIdsWithOriginalNodeIds(rightAln, rightAlnGraph);

				resultPerThread[thread].emplace_back();
				resultPerThread[thread].back().readName = pair.first;
				resultPerThread[thread].back().scoreFraction = (double)(leftAln.score() + rightAln.score()) / (double)read.sequence.size();
				std::tie(resultPerThread[thread].back().leftExon, resultPerThread[thread].back().leftReverse) = getLastExon(graph, leftAln);
				std::tie(resultPerThread[thread].back().rightExon, resultPerThread[thread].back().rightReverse) = getFirstExon(graph, rightAln);
				resultPerThread[thread].back().leftLen = pair.second.breakpoint;
				resultPerThread[thread].back().rightLen = read.sequence.size() - pair.second.breakpoint;
				resultPerThread[thread].back().leftGene = pair.second.leftGene;
				resultPerThread[thread].back().rightGene = pair.second.rightGene;
				resultPerThread[thread].back().scoreDifference = leftAln.score() + rightAln.score() - prefixAlns.at(read.seq_id).scores.back();
				resultPerThread[thread].back().corrected = getCorrected(graph, leftAln, rightAln);
			}
		});
	}
	std::vector<FusionAlignment> result;
	for (int i = 0; i < numThreads; i++)
	{
		threads[i].join();
		result.insert(result.end(), resultPerThread[i].begin(), resultPerThread[i].end());
	}
	return result;
}

int main(int argc, char** argv)
{
	std::cerr << "Fusion finder " << VERSION << std::endl;

	std::string graphFile { argv[1] };
	std::string transcriptAlignmentFile { argv[2] };
	std::string readFile { argv[3] };
	int numThreads = std::stoi(argv[4]);
	int seedSize = std::stoi(argv[5]);
	std::string resultFusionFile { argv[6] };
	std::string correctedReadsFile { argv[7] };

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(graphFile);
	std::cerr << "load reads" << std::endl;
	auto reads = loadFastqFromFile(readFile);
	std::cerr << "load transcript alignments" << std::endl;
	auto transcripts = CommonUtils::LoadVGAlignments(transcriptAlignmentFile);
	std::cerr << "get gene belongers" << std::endl;
	auto geneBelongers = getGeneBelongers(transcripts, graph);
	std::cerr << "get extra gene-matches" << std::endl;
	auto extraGeneMatches = getExtraGeneMatches(transcripts, reads, numThreads, seedSize);
	std::cerr << "get gene order" << std::endl;
	auto geneOrder = getGeneOrder(extraGeneMatches);
	std::cerr << "get prefix alignments" << std::endl;
	auto prefixAlns = getBestPrefixAlignments(extraGeneMatches, geneOrder, graph, geneBelongers, reads, numThreads);
	std::cerr << "get suffix alignments" << std::endl;
	auto suffixAlns = getBestPrefixAlignments(extraGeneMatches, geneOrder, graph, geneBelongers, reverseComplements(reads), numThreads);
	std::cerr << "merge prefix-suffix alignments" << std::endl;
	auto pairs = mergePrefixSuffix(prefixAlns, suffixAlns, geneOrder);
	std::cerr << "get fusions" << std::endl;
	auto bestAlns = getAlnFromMerge(graph, geneBelongers, pairs, prefixAlns, reads, numThreads);
	std::cerr << "write fusions" << std::endl;
	writeFusions(bestAlns, resultFusionFile);
	std::cerr << "write corrected reads" << std::endl;
	writeCorrected(bestAlns, graph, correctedReadsFile);
}