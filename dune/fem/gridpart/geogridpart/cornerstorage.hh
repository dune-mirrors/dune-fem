#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH

#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/geometry.hh>

#include <dune/grid/geometrygrid/hostcorners.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>

#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  class IsDiscreteFunction;

  template< class, class, int, template< class > class >
  class LagrangeDiscreteFunctionSpace;



  namespace Fem
  {

    // GeoDiscreteCoordFunctionCaller
    // ------------------------------

    template< int codim, class CoordFunction, class DFSpace = typename CoordFunction::DiscreteFunctionSpaceType >
    class GeoDiscreteCoordFunctionCaller;

    template< int codim, class CoordFunction, class FunctionSpace, class GridPart, template< class > class Storage >
    class GeoDiscreteCoordFunctionCaller< codim, CoordFunction, LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1, Storage > >
    {
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1, Storage > DFSpace;
      dune_static_assert( (Conversion< DFSpace, typename CoordFunction::DiscreteFunctionSpaceType >::sameType), "Invalid use of template argument DFSpace." );

    public:
      typedef CoordFunction CoordFunctionType;

      typedef typename CoordFunctionType::GridPartType::template Codim< codim >::EntityType HostEntityType;
      typedef typename CoordFunctionType::RangeType RangeType;

      static const int dimRange = CoordFunctionType::FunctionSpaceType::dimRange;
      static const int dimension = HostEntityType::dimension;

      GeoDiscreteCoordFunctionCaller ( const CoordFunction &coordFunction,
                                       const HostEntityType &hostEntity )
      : coordFunction_( coordFunction ),
        hostEntity_( hostEntity )
      {}

      void evaluate ( unsigned int i, RangeType &y ) const
      {
        const int index = coordFunction_.gridPart().indexSet().subIndex( hostEntity_, i, dimension );
        assert( (index >= 0) && (index < (int)(coordFunction_.space().blockMapper().size())) );

        typedef typename CoordFunctionType::ConstDofBlockPtrType ConstDofBlockPtrType;
        ConstDofBlockPtrType block = coordFunction_.block( index );

        for( int k = 0; k < dimRange; ++k )
          y[ k ] = (*block)[ k ];
      }

      GeometryType type () const
      {
        return hostEntity_.type();
      }

      unsigned int numCorners () const
      {
        return GenericReferenceElements< typename CoordFunctionType::GridPartType::GridType::ctype, dimension >::general( type() ).size(dimension);
      }

    private:
      const CoordFunctionType &coordFunction_;
      const HostEntityType &hostEntity_;
    };



    // GeoCoordFunctionCaller
    // ----------------------

    template< int codim, class CoordFunction, bool discrete = Conversion< CoordFunction, Dune::IsDiscreteFunction >::exists >
    class GeoCoordFunctionCaller;

    template< int codim, class CoordFunction >
    class GeoCoordFunctionCaller< codim, CoordFunction, false >
    {
    public:
      typedef CoordFunction CoordFunctionType;

      typedef typename CoordFunctionType::GridPartType::template Codim< codim >::EntityType HostEntityType;
      typedef typename CoordFunctionType::RangeType RangeType;

      GeoCoordFunctionCaller ( const CoordFunction &coordFunction,
                               const HostEntityType &hostEntity )
      : coordFunction_( coordFunction ),
        hostCorners_( hostEntity )
      {}

      void evaluate ( unsigned int i, RangeType &y ) const
      {
        coordFunction_.evaluate( hostCorners_.corner( i ), y );
      }

      GeometryType type () const
      {
        return hostCorners_.type();
      }

      unsigned int numCorners () const
      {
        return hostCorners_.numCorners();
      }

    private:
      const CoordFunction &coordFunction_;
      GeoGrid::HostCorners< HostEntityType > hostCorners_;
    };

    template< int codim, class CoordFunction >
    class GeoCoordFunctionCaller< codim, CoordFunction, true >
    : public GeoDiscreteCoordFunctionCaller< codim, CoordFunction >
    {
      typedef GeoDiscreteCoordFunctionCaller< codim, CoordFunction > BaseType;

    public:
      typedef typename BaseType::CoordFunctionType CoordFunctionType;
      typedef typename BaseType::HostEntityType HostEntityType;

      GeoCoordFunctionCaller ( const CoordFunctionType &coordFunction,
                               const HostEntityType &hostEntity )
      : BaseType( coordFunction, hostEntity )
      {}
    };



    // GeoCoordVector
    // --------------

    template< int mydim, class GridFamily >
    class GeoCoordVector
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename Traits::CoordFunctionType CoordFunctionType;

      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      typedef GeoCoordFunctionCaller< codimension, CoordFunctionType > CoordFunctionCallerType;

    public:
      GeoCoordVector ( const CoordFunctionType &coordFunction,
                       const HostEntityType &hostEntity )
      : coordFunctionCaller_( coordFunction, hostEntity )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners == coordFunctionCaller_.numCorners() );
        for( unsigned int i = 0; i < numCorners; ++i )
          coordFunctionCaller_.evaluate( i, corners[ i ] );
      }

    private:
      const CoordFunctionCallerType coordFunctionCaller_;
    };



    // GeoLocalCoordVector
    // -------------------

    template< int mydim, class GridFamily, class LCFTraits >
    class GeoLocalCoordVector
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

    public:
      typedef LocalFunction< LCFTraits > LocalCoordFunctionType;

      explicit GeoLocalCoordVector ( const LocalCoordFunctionType &localCoordFunction )
      : localCoordFunction_( localCoordFunction )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners * dimensionworld == localCoordFunction_.numDofs() );
        for( unsigned int i = 0; i < numCorners; ++i )
        {
          for( int k = 0; k < dimensionworld; ++k )
            corners[ i ][ k ] = localCoordFunction_[ i*dimensionworld + k ];
        }
      }

    private:
      dune_static_assert( LocalCoordFunctionType::dimRange == dimensionworld, "Invalid local coordinate function." );

      const LocalCoordFunctionType &localCoordFunction_;
    };



    // IntersectionCoordVector
    // -----------------------

    template< class GridFamily >
    class GeoIntersectionCoordVector
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int codimension = 1;
      static const int mydimension = dimension-codimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometryType;
      typedef typename Traits::template Codim< codimension >::LocalGeometry HostLocalGeometryType;

      typedef typename ElementGeometryType::Implementation ElementGeometryImplType;

    public:
      GeoIntersectionCoordVector ( const ElementGeometryImplType &elementGeometry,
                                   const HostLocalGeometryType &hostLocalGeometry )
      : elementGeometry_( elementGeometry ),
        hostLocalGeometry_( hostLocalGeometry )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners == hostLocalGeometry_.corners() );
        for( unsigned int i = 0; i < numCorners; ++i )
          corners[ i ] = elementGeometry_.global( hostLocalGeometry_.corner( i ) );
      }

    private:
      const ElementGeometryImplType &elementGeometry_;
      HostLocalGeometryType hostLocalGeometry_;
    };



    // GeoCornerStorage
    // ----------------

    template< class Topology, class GridFamily >
    class GeoCornerStorage
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int mydimension = Topology::dimension;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

    public:
      static const unsigned int size = Topology::numCorners;

      template< class SubTopology >
      struct SubStorage
      {
        typedef GeoCornerStorage< SubTopology, GridFamily > type;
      };

      explicit GeoCornerStorage ( const GeoCoordVector< mydimension, GridFamily > &coords )
      {
        coords.calculate( coords_ );
      }

      template< class LCFTraits >
      explicit GeoCornerStorage ( const GeoLocalCoordVector< mydimension, GridFamily, LCFTraits > &coords )
      {
        coords.calculate( coords_ );
      }

      explicit GeoCornerStorage ( const GeoIntersectionCoordVector< GridFamily > &coords )
      {
        coords.calculate( coords_ );
      }

      template< class Mapping, unsigned int codim >
      explicit
      GeoCornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = coords[ i ];
      }

      const Coordinate &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

    private:
      Coordinate coords_[ size ];
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH
