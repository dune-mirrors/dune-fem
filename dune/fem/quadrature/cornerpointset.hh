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

    template< class ct, class Topology >
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

      typedef CornerPointList< ct, typename Dune::Impl::SimplexTopology< dim >::type > SimplexQuadratureType;
      typedef CornerPointList< ct, typename Dune::Impl::CubeTopology< dim >::type > CubeQuadratureType;
      typedef CornerPointList< ct, typename Dune::Impl::PrismTopology< pdim >::type > PrismQuadratureType;
      typedef CornerPointList< ct, typename Dune::Impl::PyramidTopology< pdim >::type > PyramidQuadratureType;

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
        assert( ! type.isNone() );
      }

      CornerPointSet ( const typename GridPart::template Codim< 0 >::EntityType &entity )
      : CornerPointSet( entity.type() )
      {
      }
    };



    // CornerPointList
    // ---------------

    template< class ct, class Topology >
    class CornerPointList
    : public IntegrationPointListImp< ct, Topology::dimension >
    {
      typedef IntegrationPointListImp< ct, Topology::dimension > BaseType;

    public:
      typedef typename BaseType::CoordinateType CoordinateType;

      explicit CornerPointList ( const size_t id );
      CornerPointList ( const GeometryType &type, const int order, const size_t id );

      int order () const { return 1; }

      static unsigned int maxOrder () { return 1; }

      GeometryType geometryType () const { return GeometryType( Topology() ); }

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
      GeometryType gt( Topology::id, Topology::dimension );
      const auto &refElement = Dune::ReferenceElements< ct, Topology::dimension >::general( gt );
      const unsigned int size = refElement.size( Topology::dimension );
      for( unsigned int i = 0; i < size; ++i )
        addIntegrationPoint( refElement.position( i, Topology::dimension ) );
    }

  } //namespace Fem

} //namespace Dune

#endif // #ifndef DUNE_FEM_CORNERPOINTSET_HH
