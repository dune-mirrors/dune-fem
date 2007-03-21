#ifndef DUNE_GRAPETUPLES_HH
#define DUNE_GRAPETUPLES_HH

#include <dune/fem/io/file/iotuple.hh>

#warning "GrapeTuple is deprecated, use IOTuple instead, see <dune/fem/io/file/iotuple.hh> "

namespace Dune {

template <class TupType>
struct GrapeTuple : public IOTuple<TupType>
{
};
  
} // end namespace Dune 
#endif
