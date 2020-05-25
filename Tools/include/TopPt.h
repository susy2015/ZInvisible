#pragma once

#include "SusyAnaTools/Tools/NTupleReader.h"
#include <vector>

namespace plotterFunctions
{
	class TopPt
	{
		private:

			void toppt(NTupleReader& tr)
			{
				if(tr.checkBranch("nGenPart"))
				{
					const auto& nGenPart = tr.getVar<unsigned int>("nGenPart");
					const auto& Id = tr.getVec<int>("GenPart_pdgId");
					const auto& pt = tr.getVec<float>("GenPart_pt");
					const auto& status = tr.getVec<int>("GenPart_statusFlags");

					float topPt = 0;
					float tbarPt = 0;

					int i;
					for(i = 0; i < nGenPart; ++i)
					{
						if((status[i] & (1<<13)) == (1<<13))
						{
							if(Id[i] == 6)
								topPt = pt[i];
							if(Id[i] == -6)
								tbarPt = pt[i];
						}
					}

					tr.registerDerivedVar("TopPt", topPt);
					tr.registerDerivedVar("TbarPt", tbarPt);
				}
			}

		public:

			TopPt() {}
			~TopPt(){}

			void operator()(NTupleReader& tr)
			{
				toppt(tr);
			}
	};
}
