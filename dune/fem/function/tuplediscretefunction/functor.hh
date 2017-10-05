#ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_FUNCTOR_HH
#define DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_FUNCTOR_HH

#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/tuplediscretefunction/dofvector.hh>


namespace Dune
{

  namespace Fem
  {

    // DofBlockFunctor
    // ---------------

    template< class ... DofVectors, class Functor >
    struct DofBlockFunctor< TupleDofVector< DofVectors ... >, Functor >
    {
      typedef TupleDofVector< DofVectors ... > DofVector;
      DofBlockFunctor ( DofVector &dofVector, Functor functor )
        : dofVector_( dofVector ), functor_( std::move( functor ) ) {}

      template< class GlobalKey >
      void operator() ( std::size_t local, const GlobalKey &globalKey ) const
      {
        const int localBlockSize
          = std::decay< decltype( std::get< GlobalKey::component() >( dofVector_ ) ) >::type::blockSize;

        const int index = globalKey.index();
        functor_( local,
                  std::get< GlobalKey::component() >( dofVector_ ) [ index / localBlockSize ][ index % localBlockSize ] );
      }

    private:
      DofVector &dofVector_;
      Functor functor_;
    };


    // DofBlockFunctor
    // ---------------

    template< class ... DofVectors, class Functor >
    struct DofBlockFunctor< const TupleDofVector< DofVectors ... >, Functor >
    {
      typedef TupleDofVector< DofVectors ... > DofVector;
      DofBlockFunctor ( const DofVector &dofVector, Functor functor )
        : dofVector_( dofVector ), functor_( std::move( functor ) ) {}

      template< class GlobalKey >
      void operator() ( std::size_t local, const GlobalKey &globalKey ) const
      {
        const int localBlockSize
          = std::decay< decltype( std::get< GlobalKey::component() >( dofVector_ ) ) >::type::blockSize;

        const int index = globalKey.index();
        functor_( local,
                  std::get< GlobalKey::component() >( dofVector_ ) [ index / localBlockSize ][ index % localBlockSize ] );
      }

    private:
      const DofVector &dofVector_;
      Functor functor_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_FUNCTOR_HH
