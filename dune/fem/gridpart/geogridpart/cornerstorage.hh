#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH

#include <array>
#include <type_traits>

#include <dune/grid/geometrygrid/hostcorners.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>


namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    class IsDiscreteFunction;

    template< class, class, int, class >
    class LagrangeDiscreteFunctionSpace;



    // GeoDiscreteCoordFunctionCaller
    // ------------------------------

    template< int codim, class CoordFunction, class DFSpace = typename CoordFunction::DiscreteFunctionSpaceType >
    class GeoDiscreteCoordFunctionCaller;

    template< int codim, class CoordFunction, class FunctionSpace, class GridPart, class Storage >
    class GeoDiscreteCoordFunctionCaller< codim, CoordFunction, LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1, Storage > >
    {
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1, Storage > DFSpace;
      static_assert( (std::is_same< DFSpace, typename CoordFunction::DiscreteFunctionSpaceType >::value), "Invalid use of template argument DFSpace." );

    public:
      typedef CoordFunction CoordFunctionType;

      typedef typename CoordFunctionType::GridPartType::template Codim< codim >::EntityType HostEntityType;
      typedef typename CoordFunctionType::RangeType RangeType;

      static const int dimRange = CoordFunctionType::FunctionSpaceType::dimRange;
      static const int dimension = HostEntityType::dimension;
      static const int mydimension = HostEntityType::mydimension;

      GeoDiscreteCoordFunctionCaller ( const CoordFunction &coordFunction,
                                       const HostEntityType &hostEntity )
      : coordFunction_( coordFunction ),
        hostEntity_( hostEntity )
      {}

      void evaluate ( unsigned int i, RangeType &y ) const
      {
        const int index = coordFunction_.gridPart().indexSet().subIndex( hostEntity_, i, dimension );
        assert( (index >= 0) && (index < (int)(coordFunction_.space().blockMapper().size())) );

        for( int k = 0; k < dimRange; ++k )
          y[ k ] = coordFunction_.dofVector()[ index ][ k ];
      }

      GeometryType type () const
      {
        return hostEntity_.type();
      }

      std::size_t numCorners () const
      {
        return ReferenceElements< typename CoordFunctionType::GridPartType::GridType::ctype, mydimension >::general( type() ).size( mydimension );
      }

    private:
      const CoordFunctionType &coordFunction_;
      const HostEntityType &hostEntity_;
    };



    // GeoCoordFunctionCaller
    // ----------------------

    template< int codim, class CoordFunction, bool discrete = std::is_convertible< CoordFunction, IsDiscreteFunction >::value >
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
        coordFunction_.evaluate( hostCorners_[ i ], y );
      }

      GeometryType type () const
      {
        return hostCorners_.type();
      }

      std::size_t numCorners () const
      {
        return hostCorners_.size();
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
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

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

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > &corners ) const
      {
        const std::size_t numCorners = coordFunctionCaller_.numCorners();
        for( std::size_t i = 0; i < numCorners; ++i )
          coordFunctionCaller_.evaluate( i, corners[ i ] );
      }

    private:
      const CoordFunctionCallerType coordFunctionCaller_;
    };



    // GeoLocalCoordVector
    // -------------------

    template< int mydim, class GridFamily, class LocalFunction >
    class GeoLocalCoordVector
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

    public:
      typedef LocalFunction LocalCoordFunctionType;

      explicit GeoLocalCoordVector ( const LocalCoordFunctionType &localCoordFunction )
      : localCoordFunction_( localCoordFunction )
      {}

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > &corners ) const
      {
        assert( (localCoordFunction_.numDofs() % dimensionworld) == 0 );
        const std::size_t numCorners = localCoordFunction_.numDofs() / dimensionworld;
        assert( size >= numCorners );
        for( std::size_t i = 0; i < numCorners; ++i )
        {
          for( int k = 0; k < dimensionworld; ++k )
            corners[ i ][ k ] = localCoordFunction_[ i*dimensionworld + k ];
        }
      }

    private:
      static_assert( LocalCoordFunctionType::dimRange == dimensionworld, "Invalid local coordinate function." );

      const LocalCoordFunctionType &localCoordFunction_;
    };



    // IntersectionCoordVector
    // -----------------------

    template< class GridFamily >
    class GeoIntersectionCoordVector
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int codimension = 1;
      static const int mydimension = dimension-codimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

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

      template< std::size_t size >
      void calculate ( std::array< Coordinate, size > &corners ) const
      {
        const std::size_t numCorners = hostLocalGeometry_.corners();
        assert( size >= numCorners );
        for( std::size_t i = 0; i < numCorners; ++i )
          corners[ i ] = elementGeometry_.global( hostLocalGeometry_.corner( i ) );
      }

    private:
      const ElementGeometryImplType &elementGeometry_;
      HostLocalGeometryType hostLocalGeometry_;
    };



    // GeoCornerStorage
    // ----------------

    template< int mydim, int cdim, class GridFamily >
    class GeoCornerStorage
    {
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      typedef FieldVector< ctype, cdim > Coordinate;

      typedef std::array< Coordinate, (1 << mydim) > Coords;

    public:
      typedef typename Coords::const_iterator const_iterator;

      explicit GeoCornerStorage ( const GeoCoordVector< mydim, GridFamily > &coords )
      {
        coords.calculate( coords_ );
      }

      template< class LCFTraits >
      explicit GeoCornerStorage ( const GeoLocalCoordVector< mydim, GridFamily, LCFTraits > &coords )
      {
        coords.calculate( coords_ );
      }

      explicit GeoCornerStorage ( const GeoIntersectionCoordVector< GridFamily > &coords )
      {
        coords.calculate( coords_ );
      }

      const Coordinate &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

      const_iterator begin () const { return coords_.begin(); }
      const_iterator end () const { return coords_.end(); }

    private:
      Coords coords_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CORNERSTORAGE_HH
