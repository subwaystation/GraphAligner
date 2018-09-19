#include <iostream>
#include <cmath>
#include "GfaGraph.h"
#include "CommonUtils.h"

constexpr int MAX_REPEAT_COUNT = 30;

std::pair<std::unordered_map<int, size_t>, size_t> getCoverages(const GfaGraph& graph, const std::vector<vg::Alignment>& alns)
{
	std::unordered_map<int, size_t> result;
	size_t totalCoverage = 0;
	for (auto aln : alns)
	{
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			result[aln.path().mapping(i).position().node_id()] += aln.path().mapping(i).edit(0).from_length();
			totalCoverage += aln.path().mapping(i).edit(0).from_length();
		}
	}
	return std::make_pair(result, totalCoverage);
}

std::vector<std::pair<int, std::vector<double>>> estimateRepeatCounts(const std::unordered_map<int, size_t>& nodeCoverages, const GfaGraph& graph, double avgCoverage)
{
	std::vector<std::pair<int, std::vector<double>>> result;
	for (auto node : graph.nodes)
	{
		if (nodeCoverages.count(node.first) == 0)
		{
			std::vector<double> posteriors;
			posteriors.push_back(1);
			for (int i = 1; i < MAX_REPEAT_COUNT; i++)
			{
				posteriors.push_back(0);
			}
			continue;
		}
		size_t nodeSize = node.second.size() - graph.edgeOverlap;
		double k = nodeCoverages.at(node.first);
		double logConditional[MAX_REPEAT_COUNT];

		double kmin1sum = 0;
		for (int i = 1; i < k; i++)
		{
			kmin1sum += log(i);
		}

		// assume repeat count 0 have a poisson distribution with mean 0.1 * avgcov
		double lambda0 = nodeSize * avgCoverage * 0.1;
		logConditional[0] = k * log(lambda0) - lambda0 - kmin1sum;

		for (int i = 1; i < MAX_REPEAT_COUNT; i++)
		{
			//assume repeat count > 0 nodes have a poisson distribution with mean repeatcount*avgcov
			double lambda = nodeSize * avgCoverage * i;
			logConditional[i] = k * log(lambda) - lambda - kmin1sum;
		}
		//arbitrarily assume ~exponentially decreasing repeat count, in reality almost certainly decreases faster
		//arbitrarily ~10% false positive kmer rate in the graph
		double priors[MAX_REPEAT_COUNT];
		priors[0] = 0.10;
		priors[1] = 0.50;
		for (int i = 2; i < MAX_REPEAT_COUNT; i++)
		{
			priors[i] = priors[i-1] * 0.5;
		}
		double logOdds[MAX_REPEAT_COUNT];
		for (int i = 0; i < MAX_REPEAT_COUNT; i++)
		{
			logOdds[i] = log(priors[i]) + logConditional[i];
		}
		double maxLogOdd = logOdds[0];
		for (int i = 0; i < MAX_REPEAT_COUNT; i++)
		{
			maxLogOdd = std::max(maxLogOdd, logOdds[i]);
		}
		double oddsSum = 0;
		for (int i = 0; i < MAX_REPEAT_COUNT; i++)
		{
			oddsSum += exp(logOdds[i] - maxLogOdd);
		}

		std::vector<double> posteriors;
		for (int i = 0; i < MAX_REPEAT_COUNT; i++)
		{
			posteriors.push_back(exp(logOdds[i] - maxLogOdd) / oddsSum);
		}

		result.emplace_back(node.first, posteriors);
	}
	return result;
}

std::vector<std::pair<int, int>> pickSafeCounts(const std::vector<std::pair<int, std::vector<double>>>& estimatedProbabilities, double safeCutoff)
{
	std::vector<std::pair<int, int>> result;
	for (auto pair : estimatedProbabilities)
	{
		double max = pair.second[0];
		int maxCount = 0;
		double sum = 0;
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (pair.second[i] > max)
			{
				max = pair.second[i];
				maxCount = i;
			}
			sum += pair.second[i];
		}
		//numeric issues..
		assert(sum >= 0.98 && sum <= 1.02);
		if (max >= sum * safeCutoff) result.emplace_back(pair.first, maxCount);
	}
	return result;
}

std::vector<std::pair<int, int>> pickTotallySafeCounts(const std::vector<std::pair<int, std::vector<double>>>& estimatedProbabilities)
{
	std::vector<std::pair<int, int>> result;
	for (auto pair : estimatedProbabilities)
	{
		int exactlyOne = 0;
		int countNotZero = 0;
		int countExactlyOne = 0;
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (pair.second[i] == 1)
			{
				exactlyOne = i;
				countExactlyOne += 1;
			}
			if (pair.second[i] != 0)
			{
				countNotZero += 1;
			}
		}
		if (countExactlyOne == 1 && countNotZero == 1) result.emplace_back(pair.first, exactlyOne);
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string alnsFile { argv[2] };
	size_t estimatedGenomeSize = std::stoi(argv[3]);

	auto graph = GfaGraph::LoadFromFile(graphFile);
	auto alns = CommonUtils::LoadVGAlignments(alnsFile);
	auto coverages = getCoverages(graph, alns);
	double avgCoverage = (double)coverages.second / (double)estimatedGenomeSize;
	std::cerr << "avg coverage " << avgCoverage << std::endl;
	auto nodeRepeatCounts = estimateRepeatCounts(coverages.first, graph, avgCoverage);
	// auto safeRepeats = pickSafeCounts(nodeRepeatCounts, 0.999);
	auto safeRepeats = pickTotallySafeCounts(nodeRepeatCounts);
	std::cerr << safeRepeats.size() << " estimated known repeat counts" << std::endl;

	std::cout << "node\tcount" << std::endl;
	for (auto pair : safeRepeats)
	{
		std::cout << pair.first << "\t" << pair.second << std::endl;
	}
}