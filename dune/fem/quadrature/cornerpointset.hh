#ifndef DUNE_FEM_CORNERPOINTSET_HH
#define DUNE_FEM_CORNERPOINTSET_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/quadrature/cachingpointlist.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridPart >
    class CornerPointSet;

    template< class ct, Dune::GeometryType::Id geometryId >
    class CornerPointList;



    // CornerPointSetTraits
    // --------------------

    template< class GridPart >
    class CornerPointSetTraits
    {
      template< class ct, int dim >
      class PointListTraits;

    public:
      typedef IntegrationPointList< typename GridPart::ctype, GridPart::dimension, PointListTraits >
        IntegrationPointListType;
      typedef typename IntegrationPointListType::CoordinateType CoordinateType;
    };



    // CornerPointSetTraits::PointListTraits
    // -------------------------------------

    template< class GridPart >
    template< class ct, int dim >
    class CornerPointSetTraits< GridPart >::PointListTraits
    {
      static const int pdim = (dim > 0 ? dim : 1);

    public:
      typedef IntegrationPointListImp< ct, dim > IntegrationPointListType;

      static constexpr Dune::GeometryType::Id simplexId = Dune::GeometryTypes::simplex(dim);
      static constexpr Dune::GeometryType::Id cubeId    = Dune::GeometryTypes::cube(dim);
      static constexpr Dune::GeometryType::Id prismId   = Dune::GeometryTypes::prism ;
      static constexpr Dune::GeometryType::Id pyramidId = Dune::GeometryTypes::pyramid;

      typedef CornerPointList< ct, simplexId >   SimplexQuadratureType;
      typedef CornerPointList< ct, cubeId    >   CubeQuadratureType;
      typedef CornerPointList< ct, prismId   >   PrismQuadratureType;
      typedef CornerPointList< ct, pyramidId >   PyramidQuadratureType;

      typedef SimplexQuadratureType PointQuadratureType;
      typedef SimplexQuadratureType LineQuadratureType;

      typedef int QuadratureKeyType ;
    };



    // CornerPointSet
    // --------------

    template< class GridPart >
    class CornerPointSet
    : public CachingPointList< GridPart, 0, CornerPointSetTraits< GridPart > >
    {
      typedef CachingPointList< GridPart, 0, CornerPointSetTraits< GridPart > > BaseType;

    public:
      CornerPointSet ( const GeometryType &type )
      : BaseType( type, 1 )
      {
        // assert( ! type.isNone() );
      }

      CornerPointSet ( const typename GridPart::template Codim< 0 >::EntityType &entity )
      : CornerPointSet( entity.type() )
      {
      }
    };



    // CornerPointList
    // ---------------

    template< class ct, Dune::GeometryType::Id geometryId >
    class CornerPointList
    : public IntegrationPointListImp< ct, Dune::GeometryType(geometryId).dim() >
    {
      typedef IntegrationPointListImp< ct, Dune::GeometryType(geometryId).dim() > BaseType;

    public:
      typedef typename BaseType::CoordinateType CoordinateType;

      explicit CornerPointList ( const size_t id );
      CornerPointList ( const GeometryType &type, const int order, const size_t id );

      int order () const { return 1; }

      static unsigned int maxOrder () { return 1; }

      GeometryType geometryType () const { return GeometryType( geometryId ); }

    protected:
      using BaseType::addIntegrationPoint;

    private:
      void initialize ();
    };



    // Implementation of CornerPointList
    // ---------------------------------

    template< class ct, Dune::GeometryType::Id geometryId >
    inline CornerPointList< ct, geometryId >::CornerPointList ( const size_t id )
    : BaseType( id )
    {
      initialize();
    }


    template< class ct, Dune::GeometryType::Id geometryId >
    inline CornerPointList< ct, geometryId >
      ::CornerPointList ( const GeometryType &type, const int order, const size_t id )
    : BaseType( id )
    {
      initialize();
    }


    template< class ct, Dune::GeometryType::Id geometryId >
    inline void CornerPointList< ct, geometryId >::initialize ()
    {
      static constexpr GeometryType gt( geometryId );
      const auto &refElement = Dune::ReferenceElements< ct, gt.dim() >::general( gt );
      const unsigned int size = refElement.size( gt.dim() );
      for( unsigned int i = 0; i < size; ++i )
        addIntegrationPoint( refElement.position( i, gt.dim() ) );
    }

  } //namespace Fem

} //namespace Dune

#endif // #ifndef DUNE_FEM_CORNERPOINTSET_HH
