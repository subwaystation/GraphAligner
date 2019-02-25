#include "GfaGraph.h"
#include "AlignmentGraph.h"
#include "BigraphToDigraph.h"

size_t graphNodeSize = 0;
size_t graphEdgeSize = 0;
bool originalGraphDeterministic;

struct Automaton
{
	std::vector<std::vector<std::pair<size_t, size_t>>> transitions;
};

struct Digraph
{
	std::string seq;
	std::vector<std::vector<size_t>> outNeighbors;
	size_t numEdges() const
	{
		size_t result = 0;
		for (size_t i = 0; i < outNeighbors.size(); i++)
		{
			result += outNeighbors[i].size();
		}
		return result;
	}
};

std::vector<size_t> getPathReachables(const std::vector<std::vector<size_t>>& outNeighbors, size_t maxPathLen)
{
	std::vector<std::tuple<size_t, int, size_t>> callStack;
	size_t i = 0;
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	std::vector<size_t> stack;
	index.resize(outNeighbors.size(), std::numeric_limits<size_t>::max());
	lowlink.resize(outNeighbors.size(), std::numeric_limits<size_t>::max());
	onStack.resize(outNeighbors.size(), false);
	size_t checknode = 0;
	std::vector<size_t> result;
	result.resize(outNeighbors.size(), std::numeric_limits<size_t>::max());
	while (true)
	{
		if (callStack.size() == 0)
		{
			while (checknode < outNeighbors.size() && index[checknode] != std::numeric_limits<size_t>::max())
			{
				checknode++;
			}
			if (checknode == outNeighbors.size()) break;
			callStack.emplace_back(checknode, 0, 0);
			checknode++;
		}
		auto top = callStack.back();
		const size_t v = std::get<0>(top);
		int state = std::get<1>(top);
		size_t w;
		size_t neighborI = std::get<2>(top);
		callStack.pop_back();
		switch(state)
		{
			case 0:
				assert(index[v] == std::numeric_limits<size_t>::max());
				assert(lowlink[v] == std::numeric_limits<size_t>::max());
				assert(!onStack[v]);
				index[v] = i;
				lowlink[v] = i;
				i += 1;
				stack.push_back(v);
				onStack[v] = true;
				[[fallthrough]];
			startloop:
			case 1:
				if (neighborI >= outNeighbors[v].size()) goto endloop;
				assert(neighborI < outNeighbors[v].size());
				w = outNeighbors[v][neighborI];
				if (index[w] == std::numeric_limits<size_t>::max())
				{
					assert(lowlink[w] == std::numeric_limits<size_t>::max());
					assert(!onStack[w]);
					callStack.emplace_back(v, 2, neighborI);
					callStack.emplace_back(w, 0, 0);
					continue;
				}
				else if (onStack[w])
				{
					lowlink[v] = std::min(lowlink[v], index[w]);
					neighborI += 1;
					goto startloop;
				}
				else
				{
					neighborI += 1;
					goto startloop;
				}
			case 2:
				assert(neighborI < outNeighbors[v].size());
				w = outNeighbors[v][neighborI];
				assert(index[w] != std::numeric_limits<size_t>::max());
				assert(lowlink[w] != std::numeric_limits<size_t>::max());
				lowlink[v] = std::min(lowlink[v], lowlink[w]);
				neighborI++;
				goto startloop;
			endloop:
			case 3:
				if (lowlink[v] == index[v])
				{
					std::vector<size_t> thisComponent;
					do
					{
						w = stack.back();
						stack.pop_back();
						onStack[w] = false;
						thisComponent.push_back(w);
					} while (w != v);
					if (thisComponent.size() >= 2)
					{
						for (auto node : thisComponent)
						{
							result[node] = maxPathLen;
						}
					}
					else
					{
						result[thisComponent[0]] = 1;
						for (auto neighbor : outNeighbors[thisComponent[0]])
						{
							if (neighbor == thisComponent[0])
							{
								result[thisComponent[0]] = maxPathLen;
								break;
							}
							assert(result[neighbor] != std::numeric_limits<size_t>::max());
							result[thisComponent[0]] = std::max(result[thisComponent[0]], result[neighbor]+1);
						}
						result[thisComponent[0]] = std::min(result[thisComponent[0]], maxPathLen);
					}
				}
		}
	}
	return result;
}

template <typename T>
void addIfUnique(std::vector<T>& c, T item)
{
	for (auto i : c)
	{
		if (i == item) return;
	}
	c.push_back(item);
}

Digraph simplerGraph(const AlignmentGraph& graph)
{
	Digraph result;
	std::unordered_map<int, size_t> nodeStart;
	std::unordered_map<int, size_t> nodeEnd;
	for (size_t i = 0; i < graph.NodeSize(); i++)
	{
		nodeStart[i] = result.seq.size();
		for (size_t j = 0; j < graph.NodeLength(i); j++)
		{
			assert(graph.NodeSequences(i, j) == 'A' || graph.NodeSequences(i, j) == 'C' || graph.NodeSequences(i, j) == 'G' || graph.NodeSequences(i, j) == 'T');
			result.seq.push_back(graph.NodeSequences(i, j));
			result.outNeighbors.emplace_back();
			if (j > 0) result.outNeighbors[result.outNeighbors.size()-2].push_back(result.outNeighbors.size()-1);
		}
		nodeEnd[i] = result.seq.size()-1;
	}
	assert(result.outNeighbors.size() == result.seq.size());
	for (size_t source = 0; source < graph.outNeighbors.size(); source++)
	{
		for (auto target : graph.outNeighbors[source])
		{
			result.outNeighbors[nodeEnd[source]].push_back(nodeStart[target]);
		}
	}
	return result;
}

Automaton toNxM(const Digraph& graph, size_t numColumns)
{
	graphNodeSize = graph.seq.size();
	graphEdgeSize = 0;
	for (auto list : graph.outNeighbors)
	{
		graphEdgeSize += list.size();
	}
	originalGraphDeterministic = true;
	for (size_t i = 0; i < graph.outNeighbors.size(); i++)
	{
		for (size_t j = 0; j < graph.outNeighbors[i].size(); j++)
		{
			for (size_t k = j+1; k < graph.outNeighbors[i].size(); k++)
			{
				if (graph.seq[graph.outNeighbors[i][j]] == graph.seq[graph.outNeighbors[i][k]]) originalGraphDeterministic = false;
			}
		}
	}

	auto pathLens = getPathReachables(graph.outNeighbors, numColumns);

	Automaton result;
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions[0].emplace_back(1, 's');
	result.transitions[3].emplace_back(2, 'e');
	std::vector<size_t> currentEqClass;
	std::vector<size_t> previousEqClass;
	currentEqClass.resize(graph.seq.size(), 3);
	previousEqClass.resize(graph.seq.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < graph.seq.size(); i++)
	{
		assert(pathLens[i] != std::numeric_limits<size_t>::max());
		assert(pathLens[i] <= numColumns);
		// if (pathLens[i] < numColumns) continue;
		previousEqClass[i] = 3;
		if (numColumns == 1)
		{
			addIfUnique(result.transitions[1], std::make_pair(currentEqClass[i], (size_t)graph.seq[i]));
		}
	}
	for (size_t j = numColumns-1; j > 0; j--)
	{
		std::map<std::tuple<size_t, size_t, size_t, size_t>, size_t> eqClasses;
		for (size_t i = 0; i < graph.seq.size(); i++)
		{
			// if (pathLens[i] < j)
			// {
			// 	//not reachable
			// 	currentEqClass[i] = std::numeric_limits<size_t>::max();
			// 	continue;
			// }
			size_t A = std::numeric_limits<size_t>::max();
			size_t C = std::numeric_limits<size_t>::max();
			size_t G = std::numeric_limits<size_t>::max();
			size_t T = std::numeric_limits<size_t>::max();
			for (auto neighbor : graph.outNeighbors[i])
			{
				if (graph.seq[neighbor] == 'A') A = previousEqClass[neighbor];
				if (graph.seq[neighbor] == 'C') C = previousEqClass[neighbor];
				if (graph.seq[neighbor] == 'G') G = previousEqClass[neighbor];
				if (graph.seq[neighbor] == 'T') T = previousEqClass[neighbor];
			}
			if (A == std::numeric_limits<size_t>::max() && C == A && G == A && T == A)
			{
				//doesn't reach a final state
				currentEqClass[i] = std::numeric_limits<size_t>::max();
				continue;
			}
			std::tuple<size_t, size_t, size_t, size_t> eqClass { A, C, G, T };
			if (eqClasses.count(eqClass) == 1)
			{
				currentEqClass[i] = eqClasses[eqClass];
				if (j == 1)
				{
					addIfUnique(result.transitions[1], std::make_pair(currentEqClass[i], (size_t)graph.seq[i]));
				}
				continue;
			}
			currentEqClass[i] = result.transitions.size();
			eqClasses[eqClass] = result.transitions.size();
			if (j == 1) result.transitions[1].emplace_back(result.transitions.size(), graph.seq[i]);
			result.transitions.emplace_back();
			if (A != std::numeric_limits<size_t>::max()) result.transitions.back().emplace_back(A, 'A');
			if (C != std::numeric_limits<size_t>::max()) result.transitions.back().emplace_back(C, 'C');
			if (G != std::numeric_limits<size_t>::max()) result.transitions.back().emplace_back(G, 'G');
			if (T != std::numeric_limits<size_t>::max()) result.transitions.back().emplace_back(T, 'T');
		}
		std::swap(currentEqClass, previousEqClass);
	}
	return result;
}

Automaton toNxMWithoutEquivalence(const Digraph& graph, size_t numColumns)
{
	auto pathLens = getPathReachables(graph.outNeighbors, numColumns);

	Automaton result;
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions.emplace_back();
	result.transitions[0].emplace_back(1, 's');
	result.transitions[2].emplace_back(3, 'e');
	std::vector<size_t> currentEqClass;
	std::vector<size_t> previousEqClass;
	currentEqClass.resize(graph.seq.size(), std::numeric_limits<size_t>::max());
	previousEqClass.resize(graph.seq.size(), std::numeric_limits<size_t>::max());
	for (size_t j = numColumns-1; j < numColumns; j--)
	{
		for (size_t i = 0; i < graph.seq.size(); i++)
		{
			currentEqClass[i] = result.transitions.size();
			result.transitions.emplace_back();
			if (j == 0) result.transitions[1].emplace_back(currentEqClass[i], graph.seq[i]);
			if (j == numColumns-1) result.transitions[currentEqClass[i]].emplace_back(2, i + 256);
			if (j < numColumns-1)
			{
				for (auto neighbor : graph.outNeighbors[i])
				{
					assert(previousEqClass[neighbor] != std::numeric_limits<size_t>::max());
					result.transitions[currentEqClass[i]].emplace_back(previousEqClass[neighbor], graph.seq[neighbor]);
				}
			}
		}
		std::swap(currentEqClass, previousEqClass);
	}
	return result;
}

Automaton toNxM(const AlignmentGraph& graph, size_t numColumns)
{
	return toNxM(simplerGraph(graph), numColumns);
}

Automaton toNxMWithoutEquivalence(const AlignmentGraph& graph, size_t numColumns)
{
	return toNxMWithoutEquivalence(simplerGraph(graph), numColumns);
}

Automaton reverse(const Automaton& original)
{
	Automaton result;
	result.transitions.resize(original.transitions.size());
	for (size_t source = 0; source < original.transitions.size(); source++)
	{
		for (auto target : original.transitions[source])
		{
			result.transitions[target.first].emplace_back(source, target.second);
		}
	}
	return result;
}

std::unordered_set<size_t> intersection(const std::unordered_set<size_t>& left, const std::unordered_set<size_t>& right)
{
	std::unordered_set<size_t> result;
	for (auto item : left)
	{
		if (right.count(item) == 1) result.insert(item);
	}
	return result;
}

std::unordered_set<size_t> difference(const std::unordered_set<size_t>& left, const std::unordered_set<size_t>& right)
{
	std::unordered_set<size_t> result;
	for (auto item : left)
	{
		if (right.count(item) == 0) result.insert(item);
	}
	return result;
}

bool equal(const std::unordered_set<size_t>& left, const std::unordered_set<size_t>& right)
{
	return left.size() == right.size() && intersection(left, right).size() == left.size();
}

// Automaton reduceGrid(const Automaton& original)
// {
// 	// assume the automaton is a topologically sorted DAG
// 	// and that 0 and 1 are the start/end states
// 	std::vector<size_t> equivalenceClass;
// 	equivalenceClass.resize(original.transitions.size(), std::numeric_limits<size_t>::max());
// 	equivalenceClass[0] = 0;
// 	equivalenceClass[1] = 1;
// 	std::mapt<std::tuple<size_t, size_t, size_t, size_t, size_t>, size_t> eqClasses;
// 	for (size_t i = 2; i < original.transitions.size(); i++)
// 	{
// 		size_t A = std::numeric_limits;
// 	}
// }

Automaton reduce(const Automaton& original)
{
	// https://en.wikipedia.org/wiki/DFA_minimization#Hopcroft's_algorithm
	// it's not really a DFA but it is "almost deterministic"
	// the equivalence class invariant is valid so it still works
	auto bw = reverse(original);
	std::vector<std::unordered_set<size_t>> P;
	std::vector<std::unordered_set<size_t>> W;
	{
		std::unordered_set<size_t> finalStates;
		std::unordered_set<size_t> nonFinalStates;
		for (size_t i = 0; i < original.transitions.size(); i++)
		{
			if (original.transitions[i].size() == 0)
			{
				finalStates.insert(i);
			}
			else
			{
				nonFinalStates.insert(i);
			}
		}
		P.push_back(nonFinalStates);
		assert(finalStates.size() > 0);
		assert(nonFinalStates.size() > 0);
		P.push_back(finalStates);
		W.push_back(finalStates);
	}
	std::vector<char> chars { 'A', 'C', 'T', 'G', 's', 'e' };
	assert(W.size() > 0);
	while (W.size() > 0)
	{
		std::unordered_set<size_t> A = W.back();
		assert(A.size() > 0);
		W.pop_back();
		for (auto c : chars)
		{
			std::unordered_set<size_t> X;
			for (auto node : A)
			{
				for (auto inNeighbor : bw.transitions[node])
				{
					if (inNeighbor.second == c) X.insert(inNeighbor.first);
				}
			}
			if (X.size() == 0) continue;
			std::vector<std::unordered_set<size_t>> newP;
			for (auto Y : P)
			{
				auto inter = intersection(Y, X);
				auto diff = difference(Y, X);
				if (inter.size() == 0 || diff.size() == 0)
				{
					newP.push_back(Y);
					continue;
				}
				newP.push_back(inter);
				newP.push_back(diff);
				bool existed = false;
				for (size_t i = 0; i < W.size(); i++)
				{
					if (equal(Y, W[i]))
					{
						std::swap(W[i], W.back());
						W.pop_back();
						W.push_back(inter);
						W.push_back(diff);
						existed = true;
						break;
					}
				}
				if (!existed)
				{
					if (inter.size() < diff.size())
					{
						W.push_back(inter);
					}
					else
					{
						W.push_back(diff);
					}
				}
			}
			P = newP;
		}
	}
	Automaton result;
	std::vector<size_t> belongsToEquivalenceClass;
	belongsToEquivalenceClass.resize(original.transitions.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < P.size(); i++)
	{
		for (auto node : P[i])
		{
			assert(belongsToEquivalenceClass[node] == std::numeric_limits<size_t>::max());
			belongsToEquivalenceClass[node] = i;
		}
	}
	result.transitions.resize(P.size());
	for (size_t i = 0; i < P.size(); i++)
	{
		assert(P[i].size() > 0);
		std::set<std::pair<size_t, char>> inNeighbors;
		for (auto node : P[i])
		{
			for (auto inNeighbor : bw.transitions[node])
			{
				assert(belongsToEquivalenceClass[inNeighbor.first] != i);
				inNeighbors.emplace(belongsToEquivalenceClass[inNeighbor.first], inNeighbor.second);
			}
		}
		for (auto neighbor : inNeighbors)
		{
			result.transitions[neighbor.first].emplace_back(i, neighbor.second);
		}
	}
	return result;
}

bool isUniqueStart(const Automaton& graph)
{
	std::unordered_set<char> starts;
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		if (graph.transitions[i].size() == 0) continue;
		if (graph.transitions[i][0].second != 's') continue;
		return graph.transitions[i].size() == 1;
	}
	assert(false);
	return false;
}

size_t nonDeterministicForks(const Automaton& graph)
{
	std::unordered_set<size_t> result;
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		for (size_t j = 0; j < graph.transitions[i].size(); j++)
		{
			for (size_t k = j+1; k < graph.transitions[i].size(); k++)
			{
				if (graph.transitions[i][j].second == graph.transitions[i][k].second)
				{
					result.insert(i);
				}
			}
		}
	}
	return result.size();
}

size_t maxWidth(const Automaton& graph)
{
	std::vector<size_t> depth;
	depth.resize(graph.transitions.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		for (auto transition : graph.transitions[i])
		{
			if (transition.second == 's') depth[i] = 0;
		}
	}
	size_t maxDepth = 0;
	while (true)
	{
		bool found = false;
		for (size_t i = 0; i < graph.transitions.size(); i++)
		{
			if (depth[i] == std::numeric_limits<size_t>::max()) continue;
			for (auto neighbor : graph.transitions[i])
			{
				if (depth[neighbor.first] == std::numeric_limits<size_t>::max())
				{
					depth[neighbor.first] = depth[i]+1;
					maxDepth = std::max(maxDepth, depth[neighbor.first]);
					found = true;
				}
				assert(depth[neighbor.first] == depth[i]+1);
			}
		}
		if (!found) break;
	}
	std::vector<size_t> width;
	width.resize(maxDepth+1, 0);
	size_t maxWidth = 0;
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		assert(depth[i] != std::numeric_limits<size_t>::max());
		assert(depth[i] < width.size());
		width[depth[i]] += 1;
		maxWidth = std::max(maxWidth, width[depth[i]]);
	}
	return maxWidth;
}

void recurseReachable(std::vector<bool>& reachable, const Automaton& graph, size_t node)
{
	if (reachable[node]) return;
	reachable[node] = true;
	for (auto transition : graph.transitions[node])
	{
		recurseReachable(reachable, graph, transition.first);
	}
}

Automaton filterReachable(const Automaton& graph, char startchar)
{
	std::vector<bool> reachable;
	reachable.resize(graph.transitions.size(), false);
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		for (auto transition : graph.transitions[i])
		{
			if (transition.second == startchar) recurseReachable(reachable, graph, i);
		}
	}
	std::vector<size_t> newIndex;
	newIndex.resize(graph.transitions.size(), std::numeric_limits<size_t>::max());
	size_t skipped = 0;
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		newIndex[i] = i - skipped;
		if (!reachable[i]) skipped++;
	}
	Automaton result;
	result.transitions.resize(graph.transitions.size() - skipped);
	for (size_t i = 0; i < graph.transitions.size(); i++)
	{
		if (!reachable[i]) continue;
		result.transitions[newIndex[i]] = graph.transitions[i];
	}
	for (size_t i = 0; i < result.transitions.size(); i++)
	{
		for (size_t j = 0; j < result.transitions[i].size(); j++)
		{
			assert(reachable[result.transitions[i][j].first]);
			result.transitions[i][j].first = newIndex[result.transitions[i][j].first];
		}
	}
	return result;
}

size_t recurseEquivalencePowersets(const Automaton& NFA, Automaton& result, const std::set<size_t>& currentSet, std::map<std::tuple<size_t, size_t, size_t, size_t, size_t>, size_t>& equivalencePowersets, std::map<std::set<size_t>, size_t>& powersetCache, size_t level, std::vector<size_t>& checksPerLevel, std::vector<size_t>& traversalsPerLevel, std::vector<size_t>& insertionsPerLevel, std::vector<size_t>& maxChecks, std::vector<size_t>& maxTraversals, std::vector<size_t>& maxInsertions, std::map<std::pair<size_t, size_t>, size_t>& maxEdgeTraversals)
{
	assert(currentSet.size() > 0);
	checksPerLevel[level]++;
	for (auto node : currentSet)
	{
		maxChecks[node]++;
	}
	if (powersetCache.count(currentSet) == 1) return powersetCache[currentSet];
	traversalsPerLevel[level]++;
	for (auto node : currentSet)
	{
		maxTraversals[node]++;
	}
	std::set<size_t> A, C, G, T, e;
	for (auto index : currentSet)
	{
		for (auto transition : NFA.transitions[index])
		{
			maxEdgeTraversals[std::make_pair(index, transition.first)]++;
			if (transition.second == 'A') A.insert(transition.first);
			if (transition.second == 'C') C.insert(transition.first);
			if (transition.second == 'G') G.insert(transition.first);
			if (transition.second == 'T') T.insert(transition.first);
			if (transition.second == 'e') e.insert(transition.first);
			if (transition.second >= 256)
			{
				size_t eqClassName = result.transitions.size();
				powersetCache[currentSet] = eqClassName;
				result.transitions.emplace_back();
				insertionsPerLevel[level]++;
				for (auto node : currentSet)
				{
					maxInsertions[node]++;
				}
				return eqClassName;
			}
		}
	}
	size_t classA = std::numeric_limits<size_t>::max(), classC = std::numeric_limits<size_t>::max(), classG = std::numeric_limits<size_t>::max(), classT = std::numeric_limits<size_t>::max(), classe = std::numeric_limits<size_t>::max();
	if (A.size() > 0)
	{
		classA = recurseEquivalencePowersets(NFA, result, A, equivalencePowersets, powersetCache, level+1, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
		// assert(classA != 0);
	}
	if (C.size() > 0)
	{
		classC = recurseEquivalencePowersets(NFA, result, C, equivalencePowersets, powersetCache, level+1, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
		// assert(classC != 0);
	}
	if (G.size() > 0)
	{
		classG = recurseEquivalencePowersets(NFA, result, G, equivalencePowersets, powersetCache, level+1, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
		// assert(classG != 0);
	}
	if (T.size() > 0)
	{
		classT = recurseEquivalencePowersets(NFA, result, T, equivalencePowersets, powersetCache, level+1, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
		// assert(classT != 0);
	}
	if (e.size() > 0)
	{
		classe = recurseEquivalencePowersets(NFA, result, e, equivalencePowersets, powersetCache, level+1, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
		// assert(classe == 0);
	}
	std::tuple<size_t, size_t, size_t, size_t, size_t> currentEqClass { classA, classC, classG, classT, classe };
	size_t eqClassName;
	if (equivalencePowersets.count(currentEqClass) == 1)
	{
		powersetCache[currentSet] = equivalencePowersets[currentEqClass];
		return equivalencePowersets[currentEqClass];
	}
	insertionsPerLevel[level]++;
	for (auto node : currentSet)
	{
		maxInsertions[node]++;
	}
	eqClassName = result.transitions.size();
	result.transitions.emplace_back();
	equivalencePowersets[currentEqClass] = eqClassName;
	powersetCache[currentSet] = eqClassName;
	if (A.size() > 0) result.transitions[eqClassName].emplace_back(classA, 'A');
	if (C.size() > 0) result.transitions[eqClassName].emplace_back(classC, 'C');
	if (G.size() > 0) result.transitions[eqClassName].emplace_back(classG, 'G');
	if (T.size() > 0) result.transitions[eqClassName].emplace_back(classT, 'T');
	if (e.size() > 0) result.transitions[eqClassName].emplace_back(classe, 'e');
	return eqClassName;
}

template <typename T>
T vecMax(const std::vector<T>& vec)
{
	T result = vec[0];
	for (auto item : vec)
	{
		result = std::max(result, item);
	}
	return result;
}

template <typename T, typename U>
U mapMax(const std::map<T, U>& map)
{
	U result {};
	for (auto pair : map)
	{
		result = std::max(result, pair.second);
	}
	return result;
}

Automaton powersetDFA(const Automaton& NFA)
{
	Automaton result;
	std::vector<size_t> maxChecks;
	std::vector<size_t> maxTraversals;
	std::vector<size_t> maxInsertions;
	std::vector<size_t> checksPerLevel;
	std::vector<size_t> traversalsPerLevel;
	std::vector<size_t> insertionsPerLevel;
	std::map<std::pair<size_t, size_t>, size_t> maxEdgeTraversals;
	maxChecks.resize(NFA.transitions.size(), 0);
	maxTraversals.resize(NFA.transitions.size(), 0);
	maxInsertions.resize(NFA.transitions.size(), 0);
	checksPerLevel.resize(NFA.transitions.size(), 0);
	traversalsPerLevel.resize(NFA.transitions.size(), 0);
	insertionsPerLevel.resize(NFA.transitions.size(), 0);
	std::map<std::tuple<size_t, size_t, size_t, size_t, size_t>, size_t> equivalencePowersets;
	std::map<std::set<size_t>, size_t> powersetCache;
	std::set<size_t> starts;
	for (size_t i = 0; i < NFA.transitions.size(); i++)
	{
		for (auto transition : NFA.transitions[i])
		{
			if (transition.second == 's')
			{
				starts.insert(transition.first);
			}
		}
	}
	size_t start = recurseEquivalencePowersets(NFA, result, starts, equivalencePowersets, powersetCache, 0, checksPerLevel, traversalsPerLevel, insertionsPerLevel, maxChecks, maxTraversals, maxInsertions, maxEdgeTraversals);
	std::cerr << "check per level " << vecMax(checksPerLevel) << std::endl;
	std::cerr << "traversal per level " << vecMax(traversalsPerLevel) << std::endl;
	std::cerr << "insertion per level " << vecMax(insertionsPerLevel) << std::endl;
	std::cerr << "max check " << vecMax(maxChecks) << std::endl;
	std::cerr << "max traversal " << vecMax(maxTraversals) << std::endl;
	std::cerr << "max insertion " << vecMax(maxInsertions) << std::endl;
	std::cerr << "max edge traversals " << mapMax(maxEdgeTraversals) << std::endl;
	result.transitions.emplace_back();
	result.transitions.back().emplace_back(start, 's');
	return result;
}

void testOneGraph(std::string filename, int numColumns)
{
	auto graph = GfaGraph::LoadFromFile(filename, true);
	auto alnGraph = DirectedGraph::BuildFromGFA(graph, false);
	auto currentGraph = toNxM(alnGraph, numColumns);
	std::cerr << "input graph nodes: " << graphNodeSize << std::endl;
	std::cerr << "input graph edges: " << graphEdgeSize << std::endl;
	std::cerr << "input graph is deterministic: " << (originalGraphDeterministic ? "yes" : "no") << std::endl;
	std::cerr << "grid size: " << currentGraph.transitions.size() << std::endl;
	auto noneqGrid = toNxMWithoutEquivalence(alnGraph, numColumns);
	std::cerr << "non-equivalent grid size: " << noneqGrid.transitions.size() << std::endl;
	currentGraph = reverse(filterReachable(reverse(filterReachable(currentGraph, 's')), 'e'));
	std::cerr << "reachable size: " << currentGraph.transitions.size() << std::endl;
	noneqGrid = reverse(filterReachable(reverse(filterReachable(noneqGrid, 's')), 'e'));
	std::cerr << "non-equivalent reachable size: " << currentGraph.transitions.size() << std::endl;
	auto DFA = powersetDFA(currentGraph);
	std::cerr << "DFA size: " << DFA.transitions.size() << std::endl;
	std::cerr << "DFA number of forward non-deterministic forks: " << nonDeterministicForks(DFA) << std::endl;
	// auto minDFA = reduce(DFA);
	// std::cerr << "minDFA size: " << minDFA.transitions.size() << std::endl;
	auto noneqDFA = powersetDFA(noneqGrid);
	std::cerr << "non-equivalent DFA size: " << noneqDFA.transitions.size() << std::endl;
	std::cerr << "Final size Q=" << DFA.transitions.size() << ", m|E|+3=" << (graphEdgeSize * numColumns + 3) << ", Q/(m|E|+3)=" << ((double)DFA.transitions.size() / (double)(graphEdgeSize * numColumns + 3)) << std::endl;
	std::cerr << "Final non-equivalent size Q=" << noneqDFA.transitions.size() << ", m|E|+3=" << (graphEdgeSize * numColumns + 3) << ", Q/(m|E|+3)=" << ((double)noneqDFA.transitions.size() / (double)(graphEdgeSize * numColumns + 3)) << std::endl;
}

template <typename F>
void enumerateDeterministicGraphs(size_t V, F f)
{
	Digraph graph;
	graph.seq.resize(V, 'A');
	std::vector<size_t> nodesA, nodesC, nodesG, nodesT;
	std::vector<size_t> neighborA, neighborC, neighborG, neighborT;
	for (size_t i = 0; i < V; i++)
	{
		nodesA.push_back(i);
	}
	neighborA.resize(V, 0);
	neighborC.resize(V, 0);
	neighborG.resize(V, 0);
	neighborT.resize(V, 0);
	graph.outNeighbors.clear();
	graph.outNeighbors.resize(V);
	for (size_t i = 0; i < V; i++)
	{
		assert(neighborA[i] <= nodesA.size());
		assert(neighborC[i] <= nodesC.size());
		assert(neighborG[i] <= nodesG.size());
		assert(neighborT[i] <= nodesT.size());
		if (neighborA[i] != nodesA.size()) graph.outNeighbors[i].push_back(nodesA[neighborA[i]]);
		if (neighborC[i] != nodesC.size()) graph.outNeighbors[i].push_back(nodesC[neighborC[i]]);
		if (neighborG[i] != nodesG.size()) graph.outNeighbors[i].push_back(nodesG[neighborG[i]]);
		if (neighborT[i] != nodesT.size()) graph.outNeighbors[i].push_back(nodesT[neighborT[i]]);
	}
	f(graph);
	while (true)
	{
		bool nextChar = false;
		for (size_t i = 0; i < V; i++)
		{
			bool cont = false;
			neighborT[i]++;
			if (neighborT[i] > nodesT.size())
			{
				neighborT[i] = 0;
				neighborG[i]++;
				if (neighborG[i] > nodesG.size())
				{
					neighborG[i] = 0;
					neighborC[i]++;
					if (neighborC[i] > nodesC.size())
					{
						neighborC[i] = 0;
						neighborA[i]++;
						if (neighborA[i] > nodesA.size())
						{
							neighborA[i] = 0;
							cont = true;
						}
					}
				}
			}
			if (cont && i == V-1) nextChar = true;
			if (!cont) break;
		}
		if (nextChar)
		{
			nodesA.clear();
			nodesC.clear();
			nodesG.clear();
			nodesT.clear();
			bool finished = false;
			for (size_t i = 0; i < V; i++)
			{
				if (graph.seq[i] == 'A')
				{
					graph.seq[i] = 'C';
					break;
				}
				if (graph.seq[i] == 'C')
				{
					graph.seq[i] = 'G';
					break;
				}
				if (graph.seq[i] == 'G')
				{
					graph.seq[i] = 'T';
					break;
				}
				if (graph.seq[i] == 'T')
				{
					graph.seq[i] = 'A';
					if (i == V-1) finished = true;
				}
			}
			if (finished) break;
			for (size_t i = 0; i < V; i++)
			{
				if (graph.seq[i] == 'A') nodesA.push_back(i);
				if (graph.seq[i] == 'C') nodesC.push_back(i);
				if (graph.seq[i] == 'G') nodesG.push_back(i);
				if (graph.seq[i] == 'T') nodesT.push_back(i);
			}
			neighborA.resize(V, 0);
			neighborC.resize(V, 0);
			neighborG.resize(V, 0);
			neighborT.resize(V, 0);
		}
		graph.outNeighbors.clear();
		graph.outNeighbors.resize(V);
		for (size_t i = 0; i < V; i++)
		{
			assert(neighborA[i] <= nodesA.size());
			assert(neighborC[i] <= nodesC.size());
			assert(neighborG[i] <= nodesG.size());
			assert(neighborT[i] <= nodesT.size());
			if (neighborA[i] != nodesA.size()) graph.outNeighbors[i].push_back(nodesA[neighborA[i]]);
			if (neighborC[i] != nodesC.size()) graph.outNeighbors[i].push_back(nodesC[neighborC[i]]);
			if (neighborG[i] != nodesG.size()) graph.outNeighbors[i].push_back(nodesG[neighborG[i]]);
			if (neighborT[i] != nodesT.size()) graph.outNeighbors[i].push_back(nodesT[neighborT[i]]);
		}
		f(graph);
	}
}

void testAllGraphs(size_t maxV, size_t maxColumns)
{
	double maxQdivME = 0;
	Digraph maxGraph;
	size_t maxLength = 0;
	for (size_t i = 2; i <= maxV; i++)
	{
		enumerateDeterministicGraphs(i, [&maxQdivME, &maxGraph, &maxLength, maxColumns](const Digraph& graph){
			size_t numEdges = graph.numEdges();
			for (size_t i = 2; i <= maxColumns; i++)
			{
				auto reachable = reverse(filterReachable(reverse(filterReachable(toNxM(graph, i), 's')), 'e'));
				if (reachable.transitions.size() == 0) return;
				auto DFA = powersetDFA(reachable);
				double fraction = DFA.transitions.size() / (maxColumns * numEdges + 3);
				if (fraction > maxQdivME)
				{
					maxQdivME = fraction;
					maxGraph = graph;
					maxLength = i;
				}
			}
		});
	}
	std::cerr << "max fraction " << maxQdivME << std::endl;
	std::cerr << "with length " << maxLength << std::endl;
	std::cerr << maxGraph.seq << std::endl;
	for (size_t i = 0; i < maxGraph.outNeighbors.size(); i++)
	{
		for (auto node : maxGraph.outNeighbors[i])
		{
			std::cerr << node << " ";
		}
		std::cerr << std::endl;
	}
}

int main(int argc, char** argv)
{
	testOneGraph(argv[1], std::stoi(argv[2]));
	// testAllGraphs(std::stoi(argv[1]), std::stoi(argv[2]));
}
