#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_DOFMAPPERCODE_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_DOFMAPPERCODE_HH

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/compile.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>


namespace Dune
{

  namespace Fem
  {

    // BDMBlockMapperSingletonKey
    // --------------------------

    template< class GridPart >
    class BDMBlockMapperSingletonKey
    {
      typedef BDMBlockMapperSingletonKey< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      BDMBlockMapperSingletonKey ( const GridPartType &gridPart )
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



    // BDMDofMapperCodeFactory
    // -----------------------

    template< class GridPart, class LocalCoefficients >
    struct BDMDofMapperCodeFactory
    {
      typedef LocalCoefficients LocalCoefficientsType;

      BDMDofMapperCodeFactory ( const LocalCoefficients &localCoefficients = LocalCoefficients() )
      : localCoefficients_( localCoefficients ), emptyCode_(0,0)
      {}

      template< class Field, int dim >
      Dune::Fem::DofMapperCode operator() ( const Dune::ReferenceElement< Field, dim > &refElement ) const
      {
        Dune::GeometryType type( GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId, dim );

        if( type == refElement.type() )
//        if( GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId == refElement.type().id() )
          return Dune::Fem::compile( refElement, localCoefficients() );
        else
          return emptyCode_;
      }

      const LocalCoefficientsType &localCoefficients () const
      {
        return localCoefficients_;
      }

    private:
      LocalCoefficientsType localCoefficients_;
      DofMapperCodeWriter emptyCode_;
    };



    // BDMBlockMapperFactory
    // ---------------------

    template< class GridPart, class LocalCoefficients >
    struct BDMBlockMapperFactory
    {
      typedef GridPart GridPartType;
      typedef LocalCoefficients LocalCoefficientsType;

      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      static BlockMapperType *createObject ( const BDMBlockMapperSingletonKey< GridPartType > &key )
      {
        return createObject( key.gridPart() );
      }

      static BlockMapperType *createObject ( const GridPart &gridPart )
      {
        BDMDofMapperCodeFactory< GridPart, LocalCoefficientsType > codeFactory;
        return new BlockMapperType( gridPart, codeFactory );
      }

      static void deleteObject ( BlockMapperType *blockMapper )
      {
        delete blockMapper;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_DOFMAPPERCODE_HH
