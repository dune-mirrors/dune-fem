#ifndef DUNE_FEM_CORNERPOINTSET_HH
#define DUNE_FEM_CORNERPOINTSET_HH

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referencedomain.hh>

#include <dune/fem/quadrature/cachingpointlist.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridPart >
  class CornerPointSet;

  template< class ct, class Topology >
  class CornerPointList;



  // CornerPointSetTraits
  // --------------------

  template< class GridPart >
  class CornerPointSetTraits
  {
    template< class ct, int dim >
    struct PointListTraits;

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

    typedef CornerPointList< ct, typename GenericGeometry::SimplexTopology< dim >::type > SimplexQuadratureType;
    typedef CornerPointList< ct, typename GenericGeometry::CubeTopology< dim >::type > CubeQuadratureType;
    typedef CornerPointList< ct, typename GenericGeometry::PrismTopology< pdim >::type > PrismQuadratureType;
    typedef CornerPointList< ct, typename GenericGeometry::PyramidTopology< pdim >::type > PyramidQuadratureType;

    typedef SimplexQuadratureType PointQuadratureType;
    typedef SimplexQuadratureType LineQuadratureType;
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
    {}

    CornerPointSet ( const typename GridPart::template Codim< 0 >::EntityType &entity )
    : BaseType( entity.type(), 1 )
    {}
  };



  // CornerPointList
  // ---------------

  template< class ct, class Topology >
  class CornerPointList
  : public IntegrationPointListImp< ct, Topology::dimension >
  {
    typedef IntegrationPointListImp< ct, Topology::dimension > BaseType;

    typedef GenericGeometry::ReferenceDomain< Topology > ReferenceDomain;
    typedef GenericGeometry::DuneGeometryTypeProvider< Topology::dimension, GeometryType::simplex >
      GeometryTypeProvider;

  public:
    typedef typename BaseType::CoordinateType CoordinateType;

    explicit CornerPointList ( const size_t id );
    CornerPointList ( const GeometryType &type, const int order, const size_t id );

    int order () const { return 1; }

    static unsigned int maxOrder () { return 1; }

    GeometryType geometry () const { return GeometryTypeProvider::type( Topology::id ); }

  protected:
    using BaseType::addIntegrationPoint;

  private:
    void initialize ();
  };



  // Implementation of CornerPointList
  // ---------------------------------

  template< class ct, class Topology >
  inline CornerPointList< ct, Topology >::CornerPointList ( const size_t id )
  : BaseType( id )
  {
    initialize();
  }


  template< class ct, class Topology >
  inline CornerPointList< ct, Topology >
    ::CornerPointList ( const GeometryType &type, const int order, const size_t id )
  : BaseType( id )
  {
    initialize();
  }


  template< class ct, class Topology >
  inline void CornerPointList< ct, Topology >::initialize ()
  {
    for( unsigned int i = 0; i < ReferenceDomain::numCorners; ++i )
    {
      CoordinateType pt;
      ReferenceDomain::corner( i, pt );
      addIntegrationPoint( pt );
    }
  }

}

#endif // #ifndef DUNE_FEM_CORNERPOINTSET_HH
