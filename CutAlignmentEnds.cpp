#include <fstream>
#include "stream.hpp"
#include "CommonUtils.h"

int main(int argc, char** argv)
{
	std::string inAlignments { argv[1] };
	std::string outAlignments { argv[2] };
	int minTrustableSize = std::stoi(argv[3]);

	auto alns = CommonUtils::LoadVGAlignments(inAlignments);
	std::vector<vg::Alignment> result;
	for (auto aln : alns)
	{
		int firstTrustable = aln.path().mapping_size();
		int lastTrustable = 0;
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			if (aln.path().mapping(i).edit(0).to_length() >= minTrustableSize)
			{
				firstTrustable = std::min(firstTrustable, i);
				lastTrustable = std::max(lastTrustable, i);
			}
		}
		if (lastTrustable > firstTrustable)
		{
			result.emplace_back();
			for (int i = firstTrustable; i <= lastTrustable; i++)
			{
				auto mapping = result.back().mutable_path()->add_mapping();
				mapping->mutable_position()->set_node_id(aln.path().mapping(i).position().node_id());
				mapping->mutable_position()->set_is_reverse(aln.path().mapping(i).position().is_reverse());
			}
		}
	}

	std::ofstream outfile { outAlignments };
	stream::write_buffered(outfile, result, 0);
}