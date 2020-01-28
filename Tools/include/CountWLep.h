#pragma once

#include "SusyAnaTools/Tools/NTupleReader.h"
#include <vector>

namespace plotterFunctions
{
	class CountWLep
	{
		private:

			void countWLep(NTupleReader& tr)
			{
				if(tr.checkBranch("nGenPart"))
				{
					const auto& nGenPart = tr.getVar<unsigned int>("nGenPart");
					const auto& Id = tr.getVec<int>("GenPart_pdgId");
					const auto& Mother = tr.getVec<int>("GenPart_genPartIdxMother");
					const auto& status = tr.getVec<int>("GenPart_status");

					int i;
					int nWplusLep = 0;
					int nWminusLep = 0;
					for(i = 0; i < nGenPart; ++i)
					{
						if((Mother[i] != -1) && ((status[i] == 1) || (status[i] == 2) || (status[i] == 23)))
						{
							if(Id[Mother[i]] == 24)
							{
								if(((Id[i] > 10) && (Id[i] < 20)) || ((-Id[i] > 10) && (-Id[i] < 20)))
								{
									nWplusLep += 1;
								}
							}
							else if(Id[Mother[i]] == -24)
							{
								if(((Id[i] > 10) && (Id[i] < 20)) || ((-Id[i] > 10) && (-Id[i] < 20)))
								{
									nWminusLep += 1;
								}
							}
						}
					}

					tr.registerDerivedVar("nWplusLep", nWplusLep);
					tr.registerDerivedVar("nWminusLep", nWminusLep);
				}
			}

		public:

			CountWLep() {}
			~CountWLep(){}

			void operator()(NTupleReader& tr)
			{
				countWLep(tr);
			}
	};
}
