#ifndef DUNE_FEM_SPACE_GENERICDISCRETE_DOFMAPPERCODE_HH
#define DUNE_FEM_SPACE_GENERICDISCRETE_DOFMAPPERCODE_HH

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/space/dofmapper/code.hh>
#include <dune/fem/space/dofmapper/compile.hh>
#include <dune/fem/space/dofmapper/indexsetdofmapper.hh>


namespace Dune
{

  namespace Fem
  {
  
    // GenericBlockMapperSingletonKey
    // -------------------------------------

    template< class GridPart >
    class GenericBlockMapperSingletonKey
    {
      typedef GenericBlockMapperSingletonKey< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      GenericBlockMapperSingletonKey ( const GridPartType &gridPart )
      : gridPart_( gridPart )
      {}

      const GridPartType &gridPart () const { return gridPart_; }
    
      bool operator== ( const ThisType &other ) const
      {
        return ( &gridPart_ == &(other.gridPart_) );
      }
      
      bool operator!= ( const ThisType &other ) const
      {
        return !( *this == other );
      }

    private:
      const GridPartType &gridPart_;
    };



    // GenericDofMapperCodeFactory
    // ----------------------------------

    template< class LocalCoefficients >
    struct GenericDofMapperCodeFactory
    {
      typedef LocalCoefficients LocalCoefficientsType;

      GenericDofMapperCodeFactory ( const LocalCoefficients &localCoefficients = LocalCoefficients() )
      : localCoefficients_( localCoefficients )
      {}

      template< class Field, int dim >
      Dune::Fem::DofMapperCode operator() ( const Dune::ReferenceElement< Field, dim > &refElement ) const
      {
        return Dune::Fem::compile( refElement, localCoefficients() );
      }

      const LocalCoefficientsType &localCoefficients () const
      {
        return localCoefficients_;
      }

    private:
      LocalCoefficientsType localCoefficients_;
    };



    // GenericBlockMapperFactory
    // --------------------------------

    template< class GridPart, class LocalCoefficients >
    struct GenericBlockMapperFactory
    {
      typedef GridPart GridPartType;
      typedef LocalCoefficients LocalCoefficientsType;

      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      static BlockMapperType *createObject ( const GenericBlockMapperSingletonKey< GridPartType > &key )
      {
        return createObject( key.gridPart() );
      }

      static BlockMapperType *createObject ( const GridPart &gridPart )
      {
        GenericDofMapperCodeFactory< LocalCoefficientsType > codeFactory;
        return new BlockMapperType( gridPart, codeFactory );
      }

      static void deleteObject ( BlockMapperType *blockMapper )
      {
        delete blockMapper;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_GENERICDISCRETE_DOFMAPPERCODE_HH
