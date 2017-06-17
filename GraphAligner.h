#ifndef GraphAligner_H
#define GraphAligner_H

//http://biorxiv.org/content/early/2017/04/06/124941
#include <queue>
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/container/flat_set.hpp>
#include <unordered_set>
#include "vg.pb.h"
#include "2dArray.h"
#include "SliceRow.h"
#include "SparseBoolMatrix.h"
#include "SparseMatrix.h"

using namespace boost;

void printtime(const char* msg)
{
	static auto time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	auto newtime = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << msg << " " << newtime << " (" << (newtime - time) << ")" << std::endl;
	time = newtime;
}

template <typename LengthType, typename ScoreType>
class GraphAligner
{
public:
	class AlignmentResult
	{
	public:
		AlignmentResult(vg::Alignment alignment, int maxDistanceFromBand, bool alignmentFailed) :
		alignment(alignment),
		maxDistanceFromBand(maxDistanceFromBand),
		alignmentFailed(alignmentFailed)
		{
		}
		vg::Alignment alignment;
		int maxDistanceFromBand;
		bool alignmentFailed;
	};
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> M;
		std::vector<ScoreType> Q;
		std::vector<ScoreType> R;
		std::vector<MatrixPosition> Rbacktrace;
		std::vector<MatrixPosition> Qbacktrace;
		std::vector<LengthType> maxScorePositionPerRow;
	};
	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};

	GraphAligner() :
	nodeStart(),
	indexToNode(),
	nodeLookup(),
	nodeIDs(),
	inNeighbors(),
	nodeSequences(),
	gapStartPenalty(1),
	gapContinuePenalty(1)
	{
		//add the start dummy node as the first node
		dummyNodeStart = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(false);
		nodeSequences.push_back('N');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
	}
	
	void AddNode(int nodeId, std::string sequence, bool reverseNode)
	{
		//subgraph extraction might produce different subgraphs with common nodes
		//don't add duplicate nodes
		if (nodeLookup.count(nodeId) != 0) return;

		assert(std::numeric_limits<LengthType>::max() - sequence.size() > nodeSequences.size());
		nodeLookup[nodeId] = nodeStart.size();
		nodeIDs.push_back(nodeId);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(reverseNode);
		nodeSequences.insert(nodeSequences.end(), sequence.begin(), sequence.end());
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeIDs.size() == nodeStart.size());
		assert(nodeStart.size() == inNeighbors.size());
		assert(inNeighbors.size() == nodeEnd.size());
		assert(nodeEnd.size() == notInOrder.size());
		assert(nodeSequences.size() == indexToNode.size());
		assert(inNeighbors.size() == outNeighbors.size());
	}
	
	void AddEdgeNodeId(int node_id_from, int node_id_to)
	{
		assert(nodeLookup.count(node_id_from) > 0);
		assert(nodeLookup.count(node_id_to) > 0);
		auto from = nodeLookup[node_id_from];
		auto to = nodeLookup[node_id_to];
		assert(to >= 0);
		assert(from >= 0);
		assert(to < inNeighbors.size());
		assert(from < nodeStart.size());

		//subgraph extraction might produce different subgraphs with common edges
		//don't add duplicate edges
		if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) != inNeighbors[to].end()) return;

		inNeighbors[to].push_back(from);
		outNeighbors[from].push_back(to);
		if (from >= to)
		{
			notInOrder[to] = true;
		}
	}

	void Finalize()
	{
		//add the end dummy node as the last node
		dummyNodeEnd = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		reverse.push_back(false);
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		nodeSequences.push_back('N');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		finalized = true;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, const std::vector<SeedHit>& seedHits) const
	{
		assert(finalized);
		auto seedsInMatrix = getSeedHitPositionsInMatrix(sequence, seedHits);
		auto trace = getBacktrace(sequence, seedsInMatrix);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace));
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

	std::vector<MatrixPosition> getSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const
	{
		std::vector<MatrixPosition> result;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			assert(nodeLookup.count(seedHits[i].nodeId) > 0);
			result.emplace_back(nodeStart[nodeLookup.at(seedHits[i].nodeId)] + seedHits[i].nodePos, seedHits[i].sequencePosition);
		}
		std::sort(result.begin(), result.end(), [](auto& left, auto& right) { return left.second < right.second; });
		return result;
	}

	AlignmentResult emptyAlignment() const
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::min());
		return AlignmentResult { result, 0, true };
	}

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, int maxDistanceFromBand) const
	{
		vg::Alignment result;
		result.set_name(seq_id);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		size_t pos = 0;
		size_t oldNode = indexToNode[trace[0].first];
		while (oldNode == dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment();
			assert(pos < trace.size());
			oldNode = indexToNode[trace[pos].first];
			assert(oldNode < nodeIDs.size());
		}
		if (oldNode == dummyNodeEnd) return emptyAlignment();
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(nodeIDs[oldNode]);
		position->set_is_reverse(reverse[oldNode]);
		for (; pos < trace.size(); pos++)
		{
			if (indexToNode[trace[pos].first] == dummyNodeEnd) break;
			if (indexToNode[trace[pos].first] == oldNode) continue;
			oldNode = indexToNode[trace[pos].first];
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(nodeIDs[oldNode]);
			position->set_is_reverse(reverse[oldNode]);
		}
		result.set_score(score);
		result.set_sequence(sequence);
		return AlignmentResult { result, maxDistanceFromBand, false };
	}

	template <typename MatrixType>
	std::tuple<ScoreType, int, std::vector<MatrixPosition>> backtraceExpandoThingy(LengthType position, const SparseMatrix<MatrixPosition>& backtraceMatrix, const MatrixType& band, const std::string& sequence) const
	{
		assert(backtraceMatrix.sizeRows() == sequence.size()+1);
		assert(backtraceMatrix.sizeColumns() == nodeSequences.size());
		std::vector<MatrixPosition> trace;
		MatrixPosition currentPosition = std::make_pair(position, sequence.size());
		assert(band(currentPosition.first, currentPosition.second));
		trace.push_back(currentPosition);
		ScoreType score = 0;
		while (currentPosition.second > 0)
		{
			// std::cerr << currentPosition.first << ", " << currentPosition.second << ", " << nodeIDs[indexToNode[currentPosition.first]] << std::endl;
			assert(band(currentPosition.first, currentPosition.second));
			assert(currentPosition.second >= 0);
			assert(currentPosition.second < sequence.size()+1);
			assert(currentPosition.first >= 0);
			assert(currentPosition.first < nodeSequences.size());
			//If we're at the dummy node, we have to stay there
			if (currentPosition.first == 0) break;
			assert(backtraceMatrix.exists(currentPosition.first, currentPosition.second));
			auto newPos = backtraceMatrix(currentPosition.first, currentPosition.second);
			if (newPos.first == currentPosition.first || newPos.second == currentPosition.second)
			{
				score--;
			}
			else if (sequence[currentPosition.second] == nodeSequences[currentPosition.first])
			{
				score++;
			}
			else
			{
				score--;
			}
			// assert(newPos.second < currentPosition.second || (newPos.second == currentPosition.second && newPos.first < currentPosition.first));
			currentPosition = newPos;
			trace.push_back(currentPosition);
		}
		std::reverse(trace.begin(), trace.end());
		return std::make_tuple(score, sequence.size() - score, trace);
	}

	class ExpandoCell
	{
	public:
		ExpandoCell(LengthType w, LengthType j, LengthType btw, LengthType btj, size_t seedHit) :
		position(w, j),
		backtrace(btw, btj),
		seedHit(seedHit)
		{}
		ExpandoCell(MatrixPosition pos, MatrixPosition bt, size_t seedHit) :
		position(pos),
		backtrace(bt),
		seedHit(seedHit)
		{}
		MatrixPosition position;
		MatrixPosition backtrace;
		size_t seedHit;
	};

	class SeedhitEdge
	{
	public:
		SeedhitEdge(const std::vector<MatrixPosition>& trace, size_t start, size_t end) :
		trace(trace),
		start(start),
		end(end)
		{}
		std::vector<MatrixPosition> trace;
		size_t start;
		size_t end;
	};

	MatrixPosition makeBacktrace(SparseMatrix<MatrixPosition>& resultBacktrace, const std::vector<SeedhitEdge>& edges, size_t startHit, size_t endHit, size_t seedhitCount) const
	{
		std::vector<const SeedhitEdge*> hitBacktrace;
		std::vector<LengthType> hitDistance;
		std::vector<std::vector<const SeedhitEdge*>> outEdges;
		assert(startHit < seedhitCount + 2);
		assert(endHit < seedhitCount + 2);
		outEdges.resize(seedhitCount+2);
		for (size_t i = 0; i < edges.size(); i++)
		{
			assert(edges[i].start < seedhitCount + 2);
			assert(edges[i].end < seedhitCount + 2);
			outEdges[edges[i].start].push_back(&(edges[i]));
		}
		hitBacktrace.resize(seedhitCount+2, nullptr);
		hitDistance.resize(seedhitCount+2, std::numeric_limits<LengthType>::max());
		hitBacktrace[startHit] = nullptr;
		hitDistance[startHit] = 0;
		bool foundOne = true;
		for (size_t j = 0; j < outEdges[startHit].size(); j++)
		{
			assert(outEdges[startHit][j]->trace.size() < hitDistance[outEdges[startHit][j]->end]);
			hitDistance[outEdges[startHit][j]->end] = hitDistance[startHit] + outEdges[startHit][j]->trace.size();
			hitBacktrace[outEdges[startHit][j]->end] = outEdges[startHit][j];
		}
		while (foundOne)
		{
			foundOne = false;
			for (size_t i = 0; i < seedhitCount; i++)
			{
				if (hitDistance[i] == std::numeric_limits<LengthType>::max()) continue;
				for (size_t j = 0; j < outEdges[i].size(); j++)
				{
					if (hitDistance[i] + outEdges[i][j]->trace.size() < hitDistance[outEdges[i][j]->end])
					{
						hitDistance[outEdges[i][j]->end] = hitDistance[i] + outEdges[i][j]->trace.size();
						hitBacktrace[outEdges[i][j]->end] = outEdges[i][j];
						foundOne = true;
					}
				}
			}
		}
		assert(hitBacktrace[endHit] != nullptr);
		assert(hitDistance[endHit] != std::numeric_limits<LengthType>::max());
		size_t currentHit = endHit;
		MatrixPosition finalPosition = std::make_pair(0, 0);
		while (currentHit != startHit)
		{
			assert(hitBacktrace[currentHit] != nullptr);
			assert(hitDistance[currentHit] != std::numeric_limits<LengthType>::max());
			if (finalPosition.first == 0 && finalPosition.second == 0 && hitBacktrace[currentHit]->trace.size() > 0)
			{
				finalPosition = hitBacktrace[currentHit]->trace.back();
			}
			for (size_t i = 0; i < hitBacktrace[currentHit]->trace.size()-1; i++)
			{
				if (hitBacktrace[currentHit]->trace[i] != hitBacktrace[currentHit]->trace[i+1])
				{
					resultBacktrace.set(hitBacktrace[currentHit]->trace[i+1].first, hitBacktrace[currentHit]->trace[i+1].second, hitBacktrace[currentHit]->trace[i]);
				}
			}
			if (hitBacktrace[currentHit]->start == startHit)
			{
				resultBacktrace.set(hitBacktrace[currentHit]->trace[0].first, hitBacktrace[currentHit]->trace[0].second, {0, 0});
			}
			currentHit = hitBacktrace[currentHit]->start;
		}
		assert(finalPosition.first != 0 && finalPosition.second != 0);
		return finalPosition;
	}

	std::vector<MatrixPosition> getBacktracePath(const SparseMatrix<MatrixPosition>& trace, LengthType w, LengthType j, bool reverse) const
	{
		MatrixPosition pos = std::make_pair(w, j);
		std::vector<MatrixPosition> result;
		result.push_back(pos);
		while (true)
		{
			assert(trace.exists(pos.first, pos.second));
			auto newPos = trace.get(pos.first, pos.second);
			if (newPos == pos) break;
			pos = newPos;
			result.push_back(pos);
		}
		if (reverse) std::reverse(result.begin(), result.end());
		return result;
	}

	void propagateReachability(std::vector<bool>& reachableFromStart, const std::vector<SeedhitEdge>& edges) const
	{
		bool foundOne = true;
		while (foundOne)
		{
			foundOne = false;
			for (size_t i = 0; i < edges.size(); i++)
			{
				if (reachableFromStart[edges[i].start] && !reachableFromStart[edges[i].end])
				{
					foundOne = true;
					reachableFromStart[edges[i].end] = true;
				}
			}
		}
	}

	template <typename MatrixType>
	std::pair<LengthType, ScoreType> getScoreAndPositionWithExpandoThingy(const std::string& sequence, SparseMatrix<MatrixPosition>& resultBacktrace, MatrixType& visited, const std::vector<MatrixPosition>& seedHits) const
	{
		std::vector<SparseBoolMatrix<SliceRow<LengthType>>> visitedBySeedHitForward;
		std::vector<SparseBoolMatrix<SliceRow<LengthType>>> visitedBySeedHitBackward;
		SparseBoolMatrix<SliceRow<LengthType>> visitedByAnyForward {nodeSequences.size() + 1, sequence.size()+1};
		SparseBoolMatrix<SliceRow<LengthType>> visitedByAnyBackward {nodeSequences.size() + 1, sequence.size()+1};
		visitedBySeedHitForward.resize(seedHits.size(), {nodeSequences.size() + 1, sequence.size()+1});
		visitedBySeedHitBackward.resize(seedHits.size(), {nodeSequences.size() + 1, sequence.size()+1});
		std::vector<SparseMatrix<MatrixPosition>> forwardBacktraces;
		forwardBacktraces.resize(seedHits.size(), {nodeSequences.size() + 1, sequence.size()+1});
		std::vector<SparseMatrix<MatrixPosition>> backwardBacktraces;
		backwardBacktraces.resize(seedHits.size(), {nodeSequences.size() + 1, sequence.size()+1});
		std::vector<SeedhitEdge> seedhitEdges;
		std::vector<bool> foundForward;
		std::vector<bool> foundBackward;
		foundForward.resize(seedHits.size(), false);
		foundBackward.resize(seedHits.size(), false);
		auto startHit = seedHits.size();
		auto endHit = seedHits.size()+1;
		std::vector<bool> reachableFromStart;
		reachableFromStart.resize(seedHits.size()+2, false);
		reachableFromStart[startHit] = true;
		std::vector<ExpandoCell> currentDistanceQueueForward;
		std::vector<ExpandoCell> plusOneDistanceQueueForward;
		std::vector<ExpandoCell> currentDistanceQueueBackward;
		std::vector<ExpandoCell> plusOneDistanceQueueBackward;
		for (size_t seedI = 0; seedI < seedHits.size(); seedI++)
		{
			auto w = seedHits[seedI].first;
			auto j = seedHits[seedI].second;
			visited.set(w, j);
			forwardBacktraces[seedI].set(w, j, {w, j});
			backwardBacktraces[seedI].set(w, j, {w, j});
			visitedBySeedHitForward[seedI].set(w, j);
			visitedBySeedHitBackward[seedI].set(w, j);
			auto nodeIndex = indexToNode[w];
			plusOneDistanceQueueForward.emplace_back(w, j+1, w, j, seedI);
			if (w == nodeEnd[nodeIndex]-1)
			{
				for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
				{
					auto u = nodeStart[outNeighbors[nodeIndex][i]];
					plusOneDistanceQueueForward.emplace_back(u, j, w, j, seedI);
					if (nodeSequences[u] == sequence[j])
					{
						currentDistanceQueueForward.emplace_back(u, j+1, w, j, seedI);
					}
					else
					{
						plusOneDistanceQueueForward.emplace_back(u, j+1, w, j, seedI);
					}
				}
			}
			else
			{
				auto u = w+1;
				plusOneDistanceQueueForward.emplace_back(u, j, w, j, seedI);
				if (nodeSequences[u] == sequence[j])
				{
					currentDistanceQueueForward.emplace_back(u, j+1, w, j, seedI);
				}
				else
				{
					plusOneDistanceQueueForward.emplace_back(u, j+1, w, j, seedI);
				}
			}
			plusOneDistanceQueueBackward.emplace_back(w, j-1, w, j, seedI);
			if (w == nodeStart[nodeIndex])
			{
				for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
				{
					auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
					plusOneDistanceQueueBackward.emplace_back(u, j, w, j, seedI);
					if (nodeSequences[u] == sequence[j])
					{
						currentDistanceQueueBackward.emplace_back(u, j-1, w, j, seedI);
					}
					else
					{
						plusOneDistanceQueueBackward.emplace_back(u, j-1, w, j, seedI);
					}
				}
			}
			else
			{
				auto u = w - 1;
				plusOneDistanceQueueBackward.emplace_back(u, j, w, j, seedI);
				if (nodeSequences[u] == sequence[j])
				{
					currentDistanceQueueBackward.emplace_back(u, j-1, w, j, seedI);
				}
				else
				{
					plusOneDistanceQueueBackward.emplace_back(u, j-1, w, j, seedI);
				}
			}
		}
		ScoreType currentDistance = 0;
		size_t processed = 0;
		while (true)
		{
			if (currentDistanceQueueForward.size() > 0)
			{
				auto picked = currentDistanceQueueForward.back();
				currentDistanceQueueForward.pop_back();
				processed++;
				auto w = picked.position.first;
				auto j = picked.position.second;
				auto seedHit = picked.seedHit;
				if (foundForward[seedHit]) continue;
				if (visitedBySeedHitForward[seedHit].get(w, j)) continue;
				visitedBySeedHitForward[seedHit].set(w, j);
				if (seedHits.size() > 1) visitedByAnyForward.set(w, j);
				forwardBacktraces[seedHit].set(w, j, picked.backtrace);
				visited.set(w, j);
				if (j == sequence.size())
				{
					foundForward[seedHit] = true;
					auto bt = getBacktracePath(forwardBacktraces[seedHit], w, j, true);
					seedhitEdges.emplace_back(bt, seedHit, endHit);
					propagateReachability(reachableFromStart, seedhitEdges);
					continue;
				}
				if (seedHits.size() > 1 && visitedByAnyBackward.get(w, j))
				{
					for (size_t i = 0; i < seedHits.size(); i++)
					{
						if (i == seedHit) continue;
						if (visitedBySeedHitBackward[i].get(w, j))
						{
							foundForward[seedHit] = true;
							auto bt1 = getBacktracePath(forwardBacktraces[seedHit], w, j, true);
							assert(backwardBacktraces[i].exists(w, j));
							auto bt2 = getBacktracePath(backwardBacktraces[i], w, j, false);
							bt1.insert(bt1.end(), bt2.begin(), bt2.end());
							seedhitEdges.emplace_back(bt1, seedHit, i);
							propagateReachability(reachableFromStart, seedhitEdges);
							continue;
						}
					}
				}
				plusOneDistanceQueueForward.emplace_back(w, j+1, w, j, seedHit);
				auto nodeIndex = indexToNode[w];
				if (w == nodeEnd[nodeIndex]-1)
				{
					for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeStart[outNeighbors[nodeIndex][i]];
						plusOneDistanceQueueForward.emplace_back(u, j, w, j, seedHit);
						if (sequence[j] == nodeSequences[u])
						{
							currentDistanceQueueForward.emplace_back(u, j+1, w, j, seedHit);
						}
						else
						{
							plusOneDistanceQueueForward.emplace_back(u, j+1, w, j, seedHit);
						}
					}
				}
				else
				{
					auto u = w+1;
					plusOneDistanceQueueForward.emplace_back(u, j, w, j, seedHit);
					if (sequence[j] == nodeSequences[u])
					{
						currentDistanceQueueForward.emplace_back(u, j+1, w, j, seedHit);
					}
					else
					{
						plusOneDistanceQueueForward.emplace_back(u, j+1, w, j, seedHit);
					}
				}
			}
			else if (currentDistanceQueueBackward.size() > 0)
			{
				auto picked = currentDistanceQueueBackward.back();
				currentDistanceQueueBackward.pop_back();
				processed++;
				auto w = picked.position.first;
				auto j = picked.position.second;
				auto seedHit = picked.seedHit;
				if (foundBackward[seedHit]) continue;
				if (visitedBySeedHitBackward[seedHit].get(w, j)) continue;
				visitedBySeedHitBackward[seedHit].set(w, j);
				if (seedHits.size() > 1) visitedByAnyBackward.set(w, j);
				backwardBacktraces[seedHit].set(w, j, picked.backtrace);
				visited.set(w, j);
				if (j <= 1)
				{
					foundBackward[seedHit] = true;
					auto bt = getBacktracePath(backwardBacktraces[seedHit], w, j, false);
					seedhitEdges.emplace_back(bt, startHit, seedHit);
					propagateReachability(reachableFromStart, seedhitEdges);
					continue;
				}
				if (seedHits.size() > 1 && visitedByAnyForward.get(w, j))
				{
					for (size_t i = 0; i < seedHits.size(); i++)
					{
						if (i == seedHit) continue;
						if (visitedBySeedHitForward[i].get(w, j))
						{
							foundBackward[seedHit] = true;
							auto bt1 = getBacktracePath(backwardBacktraces[seedHit], w, j, false);
							assert(forwardBacktraces[i].exists(w, j));
							auto bt2 = getBacktracePath(forwardBacktraces[i], w, j, true);
							bt2.insert(bt2.end(), bt1.begin(), bt1.end());
							seedhitEdges.emplace_back(bt2, i, seedHit);
							propagateReachability(reachableFromStart, seedhitEdges);
							continue;
						}
					}
				}
				plusOneDistanceQueueBackward.emplace_back(w, j-1, w, j, seedHit);
				auto nodeIndex = indexToNode[w];
				if (w == nodeStart[nodeIndex])
				{
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
						plusOneDistanceQueueBackward.emplace_back(u, j, w, j, seedHit);
						if (sequence[j-1] == nodeSequences[w])
						{
							currentDistanceQueueBackward.emplace_back(u, j-1, w, j, seedHit);
						}
						else
						{
							plusOneDistanceQueueBackward.emplace_back(u, j-1, w, j, seedHit);
						}
					}
				}
				else
				{
					auto u = w-1;
					plusOneDistanceQueueBackward.emplace_back(u, j, w, j, seedHit);
					if (sequence[j-1] == nodeSequences[w])
					{
						currentDistanceQueueBackward.emplace_back(u, j-1, w, j, seedHit);
					}
					else
					{
						plusOneDistanceQueueBackward.emplace_back(u, j-1, w, j, seedHit);
					}
				}
			}
			else
			{
				std::swap(currentDistanceQueueForward, plusOneDistanceQueueForward);
				plusOneDistanceQueueForward.clear();
				std::swap(currentDistanceQueueBackward, plusOneDistanceQueueBackward);
				plusOneDistanceQueueBackward.clear();
				currentDistance += 1;
				assert(currentDistance < sequence.size());
			}
			if (reachableFromStart[endHit]) break;
		}

		auto finalPosition = makeBacktrace(resultBacktrace, seedhitEdges, startHit, endHit, seedHits.size());

		assert(currentDistance < sequence.size());
		return std::make_pair(finalPosition.first, currentDistance);
	}

	std::tuple<ScoreType, int, std::vector<MatrixPosition>> getBacktrace(const std::string& sequence, const std::vector<MatrixPosition>& seedHits) const
	{
		SparseMatrix<MatrixPosition> backtraceMatrix {nodeSequences.size(), sequence.size() + 1};
		SparseBoolMatrix<SliceRow<LengthType>> band {nodeSequences.size() + 1, sequence.size() + 1};
		auto expandoResult = getScoreAndPositionWithExpandoThingy(sequence, backtraceMatrix, band, seedHits);
		auto result = backtraceExpandoThingy(std::get<0>(expandoResult), backtraceMatrix, band, sequence);
		auto expandoCells = band.totalOnes();
		std::cerr << "number of expando cells: " << expandoCells << std::endl;
		return result;
	}

	ScoreType gapPenalty(LengthType length) const
	{
		if (length == 0) return 0;
		return gapStartPenalty + gapContinuePenalty * (length - 1);
	}

	ScoreType matchScore(char graph, char sequence) const
	{
		return graph == sequence ? 1 : -1;
	}

	std::vector<bool> notInOrder;
	std::vector<LengthType> nodeStart;
	std::vector<LengthType> nodeEnd;
	std::vector<LengthType> indexToNode;
	std::map<int, LengthType> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<LengthType>> inNeighbors;
	std::vector<std::vector<LengthType>> outNeighbors;
	std::vector<bool> reverse;
	std::string nodeSequences;
	ScoreType gapStartPenalty;
	ScoreType gapContinuePenalty;
	LengthType dummyNodeStart = 0;
	LengthType dummyNodeEnd = 1;
	bool finalized;
};

#endif