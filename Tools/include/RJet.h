#pragma once

#include "SusyAnaTools/Tools/NTupleReader.h"
#include <vector>
#include <map>
#include <cmath>

namespace plotterFunctions
{
	class RJet
	{
		private:
			std::map<std::string, std::string> var_name_map;

			std::string get_name(std::string name)
			{
				auto search = var_name_map.find(name);
				if(search != var_name_map.end())
					return search->second;
				else
					return name;
			}

			void do_rjet(NTupleReader &tr)
			{
				const auto& nJet = tr.getVar<unsigned int>(get_name("nJet"));
				const auto& jet_phi = tr.getVec<float>(get_name("Jet_phi"));
				const auto& jet_eta = tr.getVec<float>(get_name("Jet_eta"));
				const auto& jet_pt  = tr.getVec<float>(get_name("Jet_pt"));
				const auto& jet_OK  = tr.getVec<unsigned char>(get_name("Jet_Stop0l"));
				const auto& MET_pt  = tr.getVar<float>(get_name("MET_pt"));
				//const auto& MET_phi = tr.getVar<float>(get_name("MET_phi"));
				const auto& jet_dPhi_MET = tr.getVec<float>(get_name("Jet_dPhiMET"));

				float best_dphi = 1000;
				int best_index = -1;
				float rpseudo = -999;
				float best_dR = 1000;
				int best_genindex = -1;
				float rjet = -999;
				for(int i = 0; i < nJet; ++i)
				{
					if(jet_OK[i] == 1)
					{
						if(jet_dPhi_MET[i] < best_dphi)
						{
							best_dphi = jet_dPhi_MET[i];
							best_index = i;
						}
					}
				}

				if((best_index < 0) || (best_index >= nJet))
				{
					tr.registerDerivedVar("best_dphi", best_dphi);
					tr.registerDerivedVar("best_index", best_index);
					tr.registerDerivedVar("rpseudo", rpseudo);

					tr.registerDerivedVar("best_dR", best_dR);
					tr.registerDerivedVar("best_genindex", best_genindex);
					tr.registerDerivedVar("rjet", rjet);
				}
				else
				{
					rpseudo = jet_pt[best_index] / (jet_pt[best_index] + MET_pt);

					tr.registerDerivedVar("best_dphi", best_dphi);
					tr.registerDerivedVar("best_index", best_index);
					tr.registerDerivedVar("rpseudo", rpseudo);

					if(tr.checkBranch("GenJet_pt"))
					{
						const auto& ngen = tr.getVar<unsigned int>("nGenJet");
						const auto& genJet_phi = tr.getVec<float>("GenJet_phi");
						const auto& genJet_eta = tr.getVec<float>("GenJet_eta");
						const auto& genJet_pt  = tr.getVec<float>("GenJet_pt");

						double dphi;
						double dR;

						for(int j = 0; j < ngen; ++j)
						{
							dphi = fmod(jet_phi[best_index] - genJet_phi[j], 6.283185307179586);
							if (dphi < -3.141592653589793)
								dphi += 6.283185307179586;
							if (dphi >= 3.141592653589793)
								dphi -= 6.283185307179586;

							dR = sqrt(dphi*dphi + (jet_eta[best_index] - genJet_eta[j])*(jet_eta[best_index] - genJet_eta[j]));
							if(dR < best_dR)
							{
								best_dR = dR;
								best_genindex = j;
							}
						}

						float rjet;
						if((best_genindex < 0) || (best_genindex >= ngen))
						{
							tr.registerDerivedVar("best_dR", best_dR);
							tr.registerDerivedVar("best_genindex", best_genindex);
							tr.registerDerivedVar("rjet", rjet);
						}
						else
						{
							rjet = jet_pt[best_index] / genJet_pt[best_genindex];

							tr.registerDerivedVar("best_dR", best_dR);
							tr.registerDerivedVar("best_genindex", best_genindex);
							tr.registerDerivedVar("rjet", rjet);
						}
					}
				}
			}

		public:
			RJet(std::map<std::string, std::string> var_name_map) : var_name_map(var_name_map) {}
			~RJet() {}

			void operator()(NTupleReader &tr)
			{
				do_rjet(tr);
			}
	};
}
