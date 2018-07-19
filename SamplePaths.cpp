#include <fstream>
#include "CommonUtils.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	std::string inAlignments { argv[1] };
	std::string outAlignments { argv[2] };
	int minLen = std::stoi(argv[3]);

	auto aln = CommonUtils::LoadVGAlignment(inAlignments);
	std::vector<vg::Alignment> resultAlns;
	for (int i = 0; i < aln.path().mapping_size(); i++)
	{
		int currentLen = 1;
		resultAlns.emplace_back();
		for (int j = i; j < aln.path().mapping_size() && currentLen < minLen; j++)
		{
			auto mapping = resultAlns.back().mutable_path()->add_mapping();
			mapping->mutable_position()->set_node_id(aln.path().mapping(j).position().node_id());
			mapping->mutable_position()->set_is_reverse(aln.path().mapping(j).position().is_reverse());
			assert(aln.path().mapping(j).edit_size() == 1);
			assert(aln.path().mapping(j).edit(0).from_length() > 0);
			if (j > i) currentLen += aln.path().mapping(j).edit(0).from_length();
		}
	}

	std::ofstream outfile { outAlignments };
	stream::write_buffered(outfile, resultAlns, 0);
}
