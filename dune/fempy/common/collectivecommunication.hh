#ifndef DUNE_FEMPY_COMMON_COLLECTIVECOMMUNICATION_HH
#define DUNE_FEMPY_COMMON_COLLECTIVECOMMUNICATION_HH

#include <array>
#include <vector>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpicommunication.hh>

namespace Dune
{

  namespace FemPy
  {

    template< class BinaryFunction, class C, class T >
    inline static int scan ( const Dune::Communication< C > &comm, T *inout, std::size_t size )
    {
      return 0;
    }

#if HAVE_MPI
    template< class BinaryFunction, class T >
    inline static int scan ( const Dune::Communication< MPI_Comm > &comm, T *inout, std::size_t size )
    {
      return MPI_Scan( MPI_IN_PLACE, inout, size, MPITraits< T >::getType(), (Generic_MPI_Op< T, BinaryFunction >::get()), static_cast< MPI_Comm >( comm ) );
    }
#endif // #if HAVE_MPI

    template< class BinaryFunction, class C, class T >
    inline static T scan ( const Dune::Communication< C > &comm, const T &in )
    {
      T out( in );
      scan< BinaryFunction >( &out, 1u );
      return out;
    }

    template< class BinaryFunction, class C, class T, std::size_t size >
    inline static int scan ( const Dune::Communication< C > &comm, std::array< T, size > &inout )
    {
      return scan( comm, inout.data(), size );
    }

    template< class BinaryFunction, class C, class T, class A >
    inline static int scan ( const Dune::Communication< C > &comm, std::vector< T, A > &inout )
    {
      return scan( comm, inout.data(), inout.size() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_COMMON_COLLECTIVECOMMUNICATION_HH
