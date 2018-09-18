#include <iostream>
#include <cmath>
#include "GfaGraph.h"
#include "CommonUtils.h"

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
			for (int i = 1; i < 10; i++)
			{
				posteriors.push_back(0);
			}
			continue;
		}
		size_t nodeSize = node.second.size() - graph.edgeOverlap;
		double k = nodeCoverages.at(node.first);
		double logConditional[10];
		// assume repeat count 0 nodes have normally distributed coverages with mean 0 and deviation 0.1*avgcov
		// arbitrary hack, increase the deviation with the sqrt of the node size instead of linearly, discourage long nodes from turning into repeat-0
		double deviation = 0.1 * avgCoverage * sqrt(nodeSize);
		double variance = deviation * deviation;
		logConditional[0] = -log(sqrt(2 * 3.14 * variance)) - (k*k)/(2*variance*variance);
		//extra penalty, assume only 10% of kmers are errors
		logConditional[0] += nodeSize * log(0.1);

		double kmin1sum = 0;
		for (int i = 1; i < k; i++)
		{
			kmin1sum += log(i);
		}
		//arbitrarily assume max repeat count 10
		for (int i = 1; i < 10; i++)
		{
			//assume repeat count > 0 nodes have a poisson distribution with mean repeatcount*avgcov
			double lambda = nodeSize * avgCoverage * i;
			logConditional[i] = k * log(lambda) - lambda - kmin1sum;
		}
		//arbitrarily assume ~exponentially decreasing repeat count, in reality almost certainly decreases faster
		//arbitrarily 10% false positive kmer rate in the graph
		double priors[10] { 0.10, 0.45, 0.20, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125 };
		double logOdds[10];
		for (int i = 0; i < 10; i++)
		{
			logOdds[i] = log(priors[i]) + logConditional[i];
		}
		double maxLogOdd = logOdds[0];
		for (int i = 0; i < 10; i++)
		{
			maxLogOdd = std::max(maxLogOdd, logOdds[i]);
		}
		double oddsSum = 0;
		for (int i = 0; i < 10; i++)
		{
			oddsSum += exp(logOdds[i] - maxLogOdd);
		}

		std::vector<double> posteriors;
		for (int i = 0; i < 10; i++)
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
	auto safeRepeats = pickSafeCounts(nodeRepeatCounts, 0.999);
	std::cerr << safeRepeats.size() << " estimated known repeat counts" << std::endl;
	for (auto pair : safeRepeats)
	{
		std::cout << pair.first << "," << pair.second << std::endl;
	}
}