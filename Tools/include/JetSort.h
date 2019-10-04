#pragma once

#include "SusyAnaTools/Tools/NTupleReader.h"

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort


/*
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
*/

namespace plotterFunctions
{
	class JetSort
	{
		private:
			void jetSort(NTupleReader& tr)
			{
				const auto& jet_pt = tr.getVec<float>("Jet_pt");
				const auto& jet_eta = tr.getVec<float>("Jet_eta");
				std::vector<size_t> *idx = new std::vector<size_t>(jet_pt.size());
				std::iota(idx->begin(), idx->end(), 0); // Set up a vector of indices

				auto sort_comparator = [&jet_pt, &jet_eta](size_t i1, size_t i2)
				{
					if(jet_pt[i1] == jet_pt[i2])
					{
						return fabs(jet_eta[i1]) < fabs(jet_eta[i2]);
					}
					else
					{
						return jet_pt[i1] > jet_pt[i2];
					}
				};
				std::sort(idx->begin(), idx->end(), sort_comparator);
				tr.registerDerivedVec("Jet_order", idx);
			}

		public:
			JetSort(){}
			~JetSort(){}
			void operator()(NTupleReader& tr)
			{
				jetSort(tr);
			}
	};
}
