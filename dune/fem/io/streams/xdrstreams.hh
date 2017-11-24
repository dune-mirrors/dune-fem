#ifndef DUNE_FEM_XDRSTREAMS_HH
#define DUNE_FEM_XDRSTREAMS_HH
#warning "Deprecated header, use #include <dune/fem/io/streams/binarystreams.hh> instead!"

#include <dune/fem/io/streams/binarystreams.hh>

namespace Dune
{
  namespace Fem
  {

    using XDRFileOutStream = BinaryFileOutStream;
    using XDRFileInStream  = BinaryFileInStream;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_XDRSTREAMS_HH
