// -*- C++ -*-
//
//
// QDP data parallel interface
//
// prototypes used throughout the QDP code

#ifndef QDP_UTIL_INCLUDE
#define QDP_UTIL_INCLUDE

namespace QDP {

//! Decompose a lexicographic site into coordinates
multi1d<int> crtesn(int ipos, const multi1d<int>& latt_size);

//! Calculates the lexicographic site index from the coordinate of a site
int local_site(const multi1d<int>& coord, const multi1d<int>& latt_size);

//! Unique-ify a list
multi1d<int> uniquify_list(const multi1d<int>& ll);

//! Initializer for subsets
void initDefaultSets();

//! Initializer for maps
void initDefaultMaps();

} // namespace QDP

#endif
