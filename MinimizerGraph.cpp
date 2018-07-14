#include <unordered_map>
#include <iostream>
#include <chrono>
#include "MinimizerGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

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

MinimizerGraph::MinimizerGraph(size_t k, size_t w, const GfaGraph& graph) :
k(k),
w(w)
{
	assert(k < 32);
	assert(graph.edgeOverlap > k+w);
	initTopology(graph);
}

bool MinimizerGraph::minmerCompare(size_t left, size_t right) const
{
	return left < right;
}

void add(std::vector<std::pair<size_t, size_t>>& vec, size_t minmer, size_t length)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (vec[i].first == minmer) return;
	}
	vec.emplace_back(minmer, length);
}

void MinimizerGraph::initTopology(const GfaGraph& graph)
{
	std::vector<size_t> window;
	window.resize(w);
	for (auto source : graph.nodes)
	{
		for (int i = 0; i < 2; i++)
		{
			std::string combo;
			// if (i == 0) combo = getHPC(source.second); else combo = getHPC(CommonUtils::ReverseComplement(source.second));
			if (i == 0) combo = source.second; else combo = CommonUtils::ReverseComplement(source.second);
			assert(combo.size() > k+w);
			window[0] = minmerize(combo.substr(0, k));
			size_t smallestPos = 0;
			for (size_t i = 1; i < w; i++)
			{
				window[i] = nextminmer(window[i-1], combo[k+i-1]);
				if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
			}
			std::vector<std::pair<size_t, size_t>> minmers;
			minmers.emplace_back(window[smallestPos], smallestPos);
			for (size_t combopos = k+w; combopos < combo.size(); combopos++)
			{
				size_t windowpos = (combopos - k) % w;
				size_t newMinmer = nextminmer(window[(windowpos + w - 1) % w], combo[combopos]);
				if (smallestPos == windowpos)
				{
					size_t oldSmallest = window[smallestPos];
					size_t oldSmallestPos = smallestPos;
					window[windowpos] = newMinmer;
					for (size_t i = 0; i < w; i++)
					{
						if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
					}
					size_t dist = (smallestPos + w - oldSmallestPos) % w;
					if (dist == 0) dist = w;
					minmers.emplace_back(window[smallestPos], combopos);
					// add(minmerTopology[window[smallestPos]], oldSmallest, dist);
				}
				else if (minmerCompare(newMinmer, window[smallestPos]))
				{
					minmers.emplace_back(window[smallestPos], combopos);
					// add(minmerTopology[newMinmer], window[smallestPos], (windowpos + w - smallestPos) % w);
					smallestPos = windowpos;
					window[windowpos] = newMinmer;
				}
			}
			for (size_t i = 1; i < minmers.size(); i++)
			{
				for (size_t j = i-1; j < i; j++)
				{
					add(minmerTopology[minmers[i].first], minmers[j].first, minmers[i].second - minmers[j].second);
				}
			}
		}
	}
	size_t numEdges = 0;
	for (auto pair : minmerTopology)
	{
		numEdges += pair.second.size();
	}
	std::cerr << "number of minimizers: " << minmerTopology.size() << std::endl;
	std::cerr << "number of edges: " << numEdges << std::endl;
	std::cerr << "graph density: " << ((double)numEdges / (double)minmerTopology.size() / (double)minmerTopology.size()) << std::endl;
}

size_t MinimizerGraph::minmerize(std::string kmer) const
{
	assert(kmer.size() == k);
	size_t result;
	for (size_t i = 0; i < kmer.size(); i++)
	{
		result |= charNum(kmer[i]) << (i * 2);
	}
	return result;
}

size_t MinimizerGraph::nextminmer(size_t minmer, char newChar) const
{
	minmer >>= 2;
	minmer |= charNum(newChar) << (k * 2);
	return minmer;
}

std::vector<std::pair<size_t, size_t>> MinimizerGraph::inNeighbors(size_t minmer) const
{
	auto found = minmerTopology.find(minmer);
	if (found == minmerTopology.end())
	{
		return std::vector<std::pair<size_t, size_t>> {};
	}
	return found->second;
}

void MinimizerGraph::align(const std::string& originalSeq) const
{
	// std::string seq = getHPC(originalSeq);
	std::string seq = originalSeq;
	assert(seq.size() > k + w);
// #ifdef PRINTLENS
	auto timeStart = std::chrono::system_clock::now();
	std::vector<size_t> lens;
	lens.resize(seq.size(), 0);
// #endif
	std::unordered_map<size_t, size_t> lastMinmerPos;
	std::vector<size_t> window;
	window.resize(w);
	window[0] = minmerize(seq.substr(0, k));
	size_t smallestPos = 0;
	std::unordered_map<size_t, size_t> lastMinmerLength;
	for (size_t i = 1; i < w; i++)
	{
		window[i] = nextminmer(window[i-1], seq[i]);
		if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
	}
	for (size_t seqPos = k+w; seqPos < seq.size(); seqPos++)
	{
		size_t windowpos = (seqPos - k) % w;
		size_t newMinmer = nextminmer(window[(windowpos + w - 1) % w], seq[seqPos]);
		bool check = false;
		window[windowpos] = newMinmer;
		if (smallestPos == windowpos)
		{
			check = true;
			for (size_t i = 0; i < w; i++)
			{
				if (minmerCompare(window[i], window[smallestPos])) smallestPos = i;
			}
		}
		else if (minmerCompare(newMinmer, window[smallestPos]))
		{
			check = true;
		}
		if (check)
		{
			auto neighbors = inNeighbors(newMinmer);
			size_t currentLen = 0;
// #ifdef PRINTLENS
			size_t lastPos = 0;
// #endif
			for (size_t i = 0; i < neighbors.size(); i++)
			{
				size_t lastLen = lastMinmerLength[neighbors[i].first];
				int seqDist = seqPos - lastMinmerPos[neighbors[i].first];
				int topoDistance = neighbors[i].second;
				assert(topoDistance <= k+w);
				bool valid = true;
				if (seqDist > (k+w)*2) valid = false;
				if (valid && topoDistance + lastLen > currentLen)
				{
					currentLen = topoDistance + lastLen;
// #ifdef PRINTLENS
					lastPos = lastMinmerPos[neighbors[i].first];
// #endif
				}
			}
			lastMinmerPos[newMinmer] = seqPos;
			lastMinmerLength[newMinmer] = currentLen;
// #ifdef PRINTLENS
			for (size_t i = lastPos; i < seqPos; i++)
			{
				lens[i] = std::max(lens[i], currentLen);
			}
// #endif
		}
	}
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
