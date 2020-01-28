#pragma once

#include "SusyAnaTools/Tools/NTupleReader.h"
#include <vector>

namespace plotterFunctions
{
	class DownsizeBootstrap
	{
		private:
			void downsize(NTupleReader &tr)
			{
				if(tr.checkBranch("bootstrapWeight"))
				{
					const auto& original_weights = tr.getVec<int>("bootstrapWeight");
					std::vector<unsigned char> *new_weights = new std::vector<unsigned char>;

					for(const auto& weight : original_weights)
						new_weights->push_back((unsigned char) weight);

					tr.registerDerivedVec("bootstrapWeight_smaller", new_weights);
				}
			}

		public:

			DownsizeBootstrap() {}
			~DownsizeBootstrap() {}

			void operator()(NTupleReader &tr)
			{
				downsize(tr);
			}
	};
}
