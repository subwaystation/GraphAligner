#include <unordered_map>
#include <string>
#include <vector>
#include "CommonUtils.h"
#include "fastqloader.h"
#include "vg.pb.h"
#include "stream.hpp"

const int UnalignedNodeId = -1;

std::unordered_map<std::string, std::vector<vg::Alignment>> splitPerRead(const std::vector<vg::Alignment>& partials)
{
	std::unordered_map<std::string, std::vector<vg::Alignment>> result;
	for (auto aln : partials)
	{
		result[aln.name()].push_back(aln);
	}
	return result;
}

vg::Alignment combine(std::vector<vg::Alignment> partials, const std::string& seq)
{
	std::sort(partials.begin(), partials.end(), [](const vg::Alignment& left, const vg::Alignment& right) { return left.query_position() < right.query_position(); });
	vg::Alignment result;
	size_t readpos = 0;
	result.set_name(partials[0].name());
	for (size_t i = 0; i < partials.size(); i++)
	{
		if (partials[i].query_position() > readpos)
		{
			auto mapping = result.mutable_path()->add_mapping();
			mapping->mutable_position()->set_node_id(UnalignedNodeId);
			auto edit = mapping->add_edit();
			edit->set_sequence(seq.substr(readpos, partials[i].query_position() - readpos));
			edit->set_to_length(partials[i].query_position()-readpos);
			readpos = partials[i].query_position();
		}
		for (int j = 0; j < partials[i].path().mapping_size(); j++)
		{
			auto mapping = result.mutable_path()->add_mapping();
			*mapping = partials[i].path().mapping(j);
			readpos += partials[i].path().mapping(j).edit(0).to_length();
		}
	}
	if (readpos < seq.size())
	{
		auto mapping = result.mutable_path()->add_mapping();
		mapping->mutable_position()->set_node_id(UnalignedNodeId);
		auto edit = mapping->add_edit();
		edit->set_sequence(seq.substr(readpos));
		edit->set_to_length(seq.size()-readpos);
	}
	return result;
}

std::vector<vg::Alignment> combinePerRead(const std::vector<vg::Alignment>& partials, const std::vector<FastQ>& reads)
{
	auto perRead = splitPerRead(partials);
	std::vector<vg::Alignment> result;
	for (auto read : reads)
	{
		auto vec = perRead[read.seq_id];
		if (vec.size() == 0) continue;
		auto readresult = combine(vec, read.sequence);
		result.push_back(readresult);
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string inputAlns { argv[1] };
	std::string readfile { argv[2] };
	std::string outputAlns { argv[3] };

	auto alns = CommonUtils::LoadVGAlignments(inputAlns);
	auto reads = loadFastqFromFile(readfile);
	auto result = combinePerRead(alns, reads);

	std::ofstream resultFile { outputAlns, std::ios::out | std::ios::binary };
	stream::write_buffered(resultFile, result, 0);
}