#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH

#include <type_traits>

#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/geometrygrid/geometry.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/geogridpart/cornerstorage.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoGeometryTraits
    // -----------------

    template< class GridFamily >
    class GeoGeometryTraits
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      typedef Dune::Impl::FieldMatrixHelper< ctype > MatrixHelper;

      static ctype tolerance () { return 16 * std::numeric_limits< ctype >::epsilon(); }

      template< int mydim, int cdim >
      struct CornerStorage
      {
        typedef GeoCornerStorage< mydim, cdim, GridFamily > Type;
      };

      template< int mydim >
      struct hasSingleGeometryType
      : public Dune::GeoGrid::InferHasSingleGeometryType< GridPartCapabilities::hasSingleGeometryType< HostGridPartType >, dimension, mydim >
      {};
    };



    // GeoGeometry
    // -----------

    template< int mydim, int cdim, class GridFamily >
    class GeoGeometry
    {
      typedef GeoGeometry< mydim, cdim, GridFamily > ThisType;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef CachedMultiLinearGeometry< ctype, mydimension, coorddimension, GeoGeometryTraits< GridFamily > > Mapping;

    public:
      typedef typename Mapping::LocalCoordinate LocalCoordinate;
      typedef typename Mapping::GlobalCoordinate GlobalCoordinate;

      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      typedef typename Mapping::JacobianInverseTransposed JacobianInverseTransposed;

      GeoGeometry ()
      : mapping_( nullptr )
      {}

      template< class CoordVector >
      GeoGeometry ( const GeometryType &type, const CoordVector &coords )
      {
        assert( int( type.dim() ) == mydimension );
        mapping_ = new( mappingStorage_ ) Mapping( type, coords );
      }

      GeoGeometry ( const ThisType &other )
      {
        if( other.mapping_ )
          mapping_ = new( mappingStorage_ ) Mapping( *other.mapping_ );
        else
          mapping_ = nullptr;
      }

      ~GeoGeometry ()
      {
        if( mapping_ )
          mapping_->~Mapping();
      }

      const ThisType &operator= ( const ThisType &other )
      {
        if( mapping_ )
          mapping_->~Mapping();
        if( other.mapping_ )
          mapping_ = new( mappingStorage_ ) Mapping( *other.mapping_ );
        else
          mapping_ = nullptr;
        return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool affine () const { return mapping_->affine(); }
      GeometryType type () const { return mapping_->type(); }

      int corners () const { return mapping_->corners(); }
      GlobalCoordinate corner ( const int i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center () const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_->local( global ); }

      ctype integrationElement ( const LocalCoordinate &local ) const { return mapping_->integrationElement( local ); }
      ctype volume () const { return mapping_->volume(); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianInverseTransposed( local ); }

    private:
      Mapping* mapping_;
      char mappingStorage_[ sizeof( Mapping ) ];
    };

  } // namespace Fem


} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH
