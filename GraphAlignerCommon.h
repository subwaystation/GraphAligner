#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

#include <vector>
#include "AlignmentGraph.h"
#include "ArrayPriorityQueue.h"
#include "NodeSlice.h"
#include "WordSlice.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerCommon
{
public:
	class NodeWithPriority
	{
	public:
		NodeWithPriority(LengthType node, size_t offset, size_t endOffset, int priority) : node(node), offset(offset), endOffset(endOffset), priority(priority) {}
		bool operator>(const NodeWithPriority& other) const
		{
		return priority > other.priority;
		}
		bool operator<(const NodeWithPriority& other) const
		{
		return priority < other.priority;
		}
		LengthType node;
		size_t offset;
		size_t endOffset;
		int priority;
	};
	class EdgeWithPriority
	{
	public:
		EdgeWithPriority(LengthType target, size_t offset, size_t endOffset, int priority, WordSlice<LengthType, ScoreType, Word> incoming, WordSlice<LengthType, ScoreType, Word> incomingUp) : target(target), offset(offset), endOffset(endOffset), priority(priority), incoming(incoming), incomingUp(incomingUp) {}
		bool operator>(const EdgeWithPriority& other) const
		{
			return priority > other.priority;
		}
		bool operator<(const EdgeWithPriority& other) const
		{
			return priority < other.priority;
		}
		LengthType target;
		size_t offset;
		size_t endOffset;
		int priority;
		WordSlice<LengthType, ScoreType, Word> incoming;
		WordSlice<LengthType, ScoreType, Word> incomingUp;
	};
	class AlignerGraphsizedState
	{
	public:
		AlignerGraphsizedState(const AlignmentGraph& graph, int maxBandwidth) :
		calculableQueue(WordConfiguration<Word>::WordSize + maxBandwidth + 1),
		nodesliceMap(),
		currentBand(),
		previousBand()
		{
			nodesliceMap.resize(graph.NodeSize(), {0, 0, 0, 0, 0});
			currentBand.resize(graph.NodeSize(), false);
			previousBand.resize(graph.NodeSize(), false);
		}
		void clear()
		{
			calculableQueue.clear();
			nodesliceMap.assign(nodesliceMap.size(), {0, 0, 0, 0, 0});
			currentBand.assign(currentBand.size(), false);
			previousBand.assign(previousBand.size(), false);
		}
		ArrayPriorityQueue<EdgeWithPriority> calculableQueue;
		std::vector<typename NodeSlice<typename WordContainer<LengthType, ScoreType, Word>::Slice>::MapItem> nodesliceMap;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
	};
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class Params
	{
	public:
		//cutoff for doing the backtrace in the sqrt-slice pass
		//"bulges" in the band are responsible for almost all of the time spent aligning,
		//and this way they don't need to be recalculated, saving about half of the time.
		//semi-arbitrarily twenty five thousand, empirically a good enough cutoff
		static constexpr size_t BacktraceOverrideCutoff = 25000;
		Params(LengthType initialBandwidth, LengthType rampBandwidth, const AlignmentGraph& graph, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations) :
		initialBandwidth(initialBandwidth),
		rampBandwidth(rampBandwidth),
		graph(graph),
		maxCellsPerSlice(maxCellsPerSlice),
		quietMode(quietMode),
		sloppyOptimizations(sloppyOptimizations)
		{
		}
		const size_t maxCellsPerSlice;
		const LengthType initialBandwidth;
		const LengthType rampBandwidth;
		const AlignmentGraph& graph;
		const bool quietMode;
		const bool sloppyOptimizations;
	};
	class SeedHit
	{
	public:
		SeedHit(int nodeID, size_t seqPos, bool reverse) :
		nodeID(nodeID),
		seqPos(seqPos),
		reverse(reverse)
		{
		}
		int nodeID;
		size_t seqPos;
		bool reverse;
	};
	class OnewayTrace
	{
	public:
		OnewayTrace() :
		trace(),
		score(0)
		{
		}
		static OnewayTrace TraceFailed()
		{
			OnewayTrace result;
			result.score = std::numeric_limits<ScoreType>::max();
			return result;
		}
		bool failed() const
		{
			return score == std::numeric_limits<ScoreType>::max();
		}
		std::vector<MatrixPosition> trace;
		ScoreType score;
	};
	class Trace
	{
	public:
		OnewayTrace forward;
		OnewayTrace backward;
	};
	static bool characterMatch(char sequenceCharacter, char graphCharacter)
	{
		assert(graphCharacter == 'A' || graphCharacter == 'T' || graphCharacter == 'C' || graphCharacter == 'G');
		switch(sequenceCharacter)
		{
			case 'A':
			case 'a':
			return graphCharacter == 'A';
			break;
			case 'T':
			case 't':
			return graphCharacter == 'T';
			break;
			case 'C':
			case 'c':
			return graphCharacter == 'C';
			break;
			case 'G':
			case 'g':
			return graphCharacter == 'G';
			break;
			case 'N':
			case 'n':
			return true;
			break;
			case 'R':
			case 'r':
			return graphCharacter == 'A' || graphCharacter == 'G';
			break;
			case 'Y':
			case 'y':
			return graphCharacter == 'C' || graphCharacter == 'T';
			break;
			case 'K':
			case 'k':
			return graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'M':
			case 'm':
			return graphCharacter == 'C' || graphCharacter == 'A';
			break;
			case 'S':
			case 's':
			return graphCharacter == 'C' || graphCharacter == 'G';
			break;
			case 'W':
			case 'w':
			return graphCharacter == 'A' || graphCharacter == 'T';
			break;
			case 'B':
			case 'b':
			return graphCharacter == 'C' || graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'D':
			case 'd':
			return graphCharacter == 'A' || graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'H':
			case 'h':
			return graphCharacter == 'A' || graphCharacter == 'C' || graphCharacter == 'T';
			break;
			case 'V':
			case 'v':
			return graphCharacter == 'A' || graphCharacter == 'C' || graphCharacter == 'G';
			break;
			default:
			assert(false);
			std::abort();
			return false;
			break;
		}
	}
};

#endif
