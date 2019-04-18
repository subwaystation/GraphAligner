#include <unordered_map>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "CommonUtils.h"
#include "stream.hpp"

struct Mapping
{
	std::unordered_map<int, NodePos> nodeMapping;
	std::unordered_set<std::pair<NodePos, NodePos>> keptEdges;
};

vg::Alignment translate(const vg::Alignment& input, const Mapping& mapping)
{
	vg::Alignment result;
	result.set_name(input.name());
	auto vgmapping = result.mutable_path()->add_mapping();
	int start = 0;
	NodePos translated;
	while (start < input.path().mapping_size() && mapping.nodeMapping.count(input.path().mapping(start).position().node_id()) == 0) start++;
	if (start == input.path().mapping_size()) return result;
	translated = mapping.nodeMapping.at(input.path().mapping(start).position().node_id());
	if (input.path().mapping(start).position().is_reverse()) translated = translated.Reverse();
	vgmapping->mutable_position()->set_node_id(translated.id);
	vgmapping->mutable_position()->set_is_reverse(!translated.end);
	for (int i = start+1; i < input.path().mapping_size(); i++)
	{
		NodePos oldPos, newPos;
		oldPos.id = input.path().mapping(i-1).position().node_id();
		oldPos.end = !input.path().mapping(i-1).position().is_reverse();
		newPos.id = input.path().mapping(i).position().node_id();
		newPos.end = !input.path().mapping(i).position().is_reverse();
		if (mapping.nodeMapping.count(input.path().mapping(i).position().node_id()) == 0) continue;
		auto newTranslated = mapping.nodeMapping.at(input.path().mapping(i).position().node_id());
		if (input.path().mapping(i).position().is_reverse()) newTranslated = newTranslated.Reverse();
		if (mapping.keptEdges.count(std::make_pair(oldPos, newPos)) == 1 || mapping.keptEdges.count(std::make_pair(newPos, oldPos)) == 1)
		{
			vgmapping = result.mutable_path()->add_mapping();
			vgmapping->mutable_position()->set_node_id(newTranslated.id);
			vgmapping->mutable_position()->set_is_reverse(!newTranslated.end);
		}
		else
		{
			assert(newTranslated == translated);
		}
		translated = newTranslated;
	}
	return result;
}

Mapping loadMapping(std::string filename)
{
	std::ifstream file { filename };
	Mapping result;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] != 'N' && line[0] != 'E') continue;
		if (line[0] == 'N')
		{
			std::string dummy, orientation;
			int fromid, toid;
			std::stringstream sstr {line};
			sstr >> dummy >> fromid >> toid >> orientation;
			assert(dummy == "N");
			assert(orientation == "+" || orientation == "-");
			result.nodeMapping[fromid] = NodePos { toid, orientation == "+" };
		}
		if (line[0] == 'E')
		{
			std::string dummy, fromorientation, toorientation;
			int fromid, toid;
			std::stringstream sstr {line};
			sstr >> dummy >> fromid >> fromorientation >> toid >> toorientation;
			assert(dummy == "E");
			assert(fromorientation == "+" || fromorientation == "-");
			assert(toorientation == "+" || toorientation == "-");
			result.keptEdges.emplace(NodePos { fromid, fromorientation == "+" }, NodePos { toid, toorientation == "+" });
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string inputAlns { argv[1] };
	std::string inputMapping { argv[2] };
	std::string outputAlns { argv[3] };
	
	auto mapping = loadMapping(inputMapping);

	std::ofstream alignmentOut { outputAlns, std::ios::out | std::ios::binary };

	std::ifstream input { inputAlns, std::ios::in | std::ios::binary };
	std::function<void(vg::Alignment&)> lambda = [&alignmentOut, &mapping](vg::Alignment& g) {
		std::vector<vg::Alignment> translated;
		translated.push_back(translate(g, mapping));
		if (translated.back().path().mapping_size() == 0) translated.pop_back();
		stream::write_buffered(alignmentOut, translated, 0);
	};
	stream::for_each(input, lambda);
}