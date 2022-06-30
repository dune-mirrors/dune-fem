#ifndef DUNE_FEM_GRIDPART_COMMON_LOCALFUNCTIONGEOMETRY_HH
#define DUNE_FEM_GRIDPART_COMMON_LOCALFUNCTIONGEOMETRY_HH

#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/fem/common/fmatrixcol.hh>

#include <dune/fem/gridpart/common/simplegeometry.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionBasicGeometry
    // --------------------------

    template< class LocalFunction >
    struct LocalFunctionBasicGeometry
    {
      typedef typename LocalFunction::EntityType Entity;
      typedef typename LocalFunction::FunctionSpaceType::RangeFieldType ctype;

      static const int mydimension = Entity::mydimension;
      static const int coorddimension = LocalFunction::FunctionSpaceType::dimRange;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

      template< class... Args, std::enable_if_t< std::is_constructible< LocalFunction, Args &&... >::value, int > = 0 >
      LocalFunctionBasicGeometry ( Args &&... args )
        : localFunction_( std::forward< Args >( args )... )
      {}

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        GlobalCoordinate ret;
        localFunction().evaluate( local, ret );
        return ret;
      }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
      {
        const auto gradFT = localFunction().entity().geometry().jacobianTransposed( local );

        typename LocalFunction::FunctionSpaceType::JacobianRangeType gradPhi;

        localFunction().jacobian( local, gradPhi );

        JacobianTransposed jacTransposed( 0 );
        for( int i = 0; i < coorddimension; ++i )
        {
          FieldMatrixColumn< JacobianTransposed > col( jacTransposed, i );
          gradFT.mv( gradPhi[ i ], col );
        }
        return jacTransposed;
      }

      const QuadratureRule< ctype, mydimension > &quadrature ( int order ) const
      {
        return QuadratureRules< ctype, mydimension >::rule( type(), order + (2*localFunction().order() + 1) );
      }

      GeometryType type () const { return localFunction().entity().type(); }

      void bind ( const Entity &entity ) { localFunction_.bind( entity ); }
      void init ( const Entity &entity ) { bind( entity ); }

      const LocalFunction &localFunction () const { return localFunction_; }

    private:
      LocalFunction localFunction_;
    };



    // LocalFunctionGeometry
    // ---------------------

    template< class LocalFunction >
    using LocalFunctionGeometry = SimpleGeometry< LocalFunctionBasicGeometry< LocalFunction > >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_LOCALFUNCTIONGEOMETRY_HH
