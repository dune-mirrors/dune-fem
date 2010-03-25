#ifndef DUNE_GRAPEDATAIO_HH
#define DUNE_GRAPEDATAIO_HH
#warning "Deprecated header and class use binarydataio.hh and BindaryDataIO instead!" 

#include <dune/fem/io/file/binarydataio.hh>

namespace Dune {

template <class GridImp>
class GrapeDataIO : public BinaryDataIO< GridImp >
{
public:
   GrapeDataIO () DUNE_DEPRECATED {};
};

}
#endif
