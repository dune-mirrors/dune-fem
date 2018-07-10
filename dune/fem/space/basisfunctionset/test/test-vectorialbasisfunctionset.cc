#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/testgrid.hh>



// Evaluate
// --------

template< class FunctionSpace >
struct Evaluate
{
  typedef FunctionSpace FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeType Vector;

  template< class BasisFunctionSet, class Point, class DofVector >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, const DofVector &dofs, Vector &vector )
  {
    basisFunctionSet.evaluateAll( x, dofs, vector );
  }

  template< class BasisFunctionSet, class Point, class Array >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, Array &array )
  {
    basisFunctionSet.evaluateAll( x, array );
  }
};



// Jacobian
// --------

template< class FunctionSpace >
struct Jacobian
{
  typedef FunctionSpace FunctionSpaceType;
  typedef typename FunctionSpaceType::JacobianRangeType Vector;

  template< class BasisFunctionSet, class Point, class DofVector >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, const DofVector &dofs, Vector &vector )
  {
    basisFunctionSet.jacobianAll( x, dofs, vector );
  }

  template< class BasisFunctionSet, class Point, class Array >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, Array &array )
  {
    basisFunctionSet.jacobianAll( x, array );
  }
};



// Hessian
// -------

template< class FunctionSpace >
struct Hessian
{
  typedef FunctionSpace FunctionSpaceType;
  typedef typename FunctionSpaceType::HessianRangeType Vector;

  template< class BasisFunctionSet, class Point, class DofVector >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, const DofVector &dofs, Vector &vector )
  {
    basisFunctionSet.hessianAll( x, dofs, vector );
  }

  template< class BasisFunctionSet, class Point, class Array >
  static void apply ( const BasisFunctionSet &basisFunctionSet,
                      const Point &x, Array &array )
  {
    basisFunctionSet.hessianAll( x, array );
  }
};



// DofConversion
// -------------

template< class DofAlignmentType >
struct DofConversion;

template< class ShapeFunctionSet, class Range >
struct DofConversion< Dune::Fem::HorizontalDofAlignment< ShapeFunctionSet, Range > >
{
  typedef Dune::Fem::HorizontalDofAlignment< ShapeFunctionSet, Range > DofAlignmentType;
  static std::size_t apply ( const DofAlignmentType &dofAlignment, std::size_t globalDof )
  {
    Dune::Fem::VerticalDofAlignment< ShapeFunctionSet, Range > conversion;
    return dofAlignment.globalDof( conversion.localDof( globalDof ) );
  }
};

template< class ShapeFunctionSet, class Range >
struct DofConversion< Dune::Fem::VerticalDofAlignment< ShapeFunctionSet, Range > >
{
  typedef Dune::Fem::VerticalDofAlignment< ShapeFunctionSet, Range > DofAlignmentType;

  static std::size_t apply ( const DofAlignmentType &, std::size_t globalDof )
  {
    return globalDof;
  }
};



// Error
// -----

template< class Vector >
struct Error;

template< class RangeFieldType, int dim >
struct Error< Dune::FieldVector< RangeFieldType, dim > >
{
  typedef Dune::FieldVector< RangeFieldType, dim > Vector;

  static RangeFieldType apply ( const Vector &a, const Vector &b )
  {
    return (a - b).two_norm();
  }
};

template< class RangeFieldType, int dimDomain, int dimRange >
struct Error< Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > >
{
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > Vector;

  static RangeFieldType apply ( const Vector &a, const Vector &b )
  {
    Vector c = a;
    c -= b;
    return c.frobenius_norm();
  }
};

template< class RangeFieldType, int dimDomain, int dimRange >
struct Error< Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >, dimRange > >
{
  typedef Dune::FieldVector< Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >, dimRange > Vector;

  static RangeFieldType apply ( const Vector &a, const Vector &b )
  {
    Dune::FieldVector< RangeFieldType, dimRange > error;
    for( int r = 0; r < dimRange; ++r )
      error[ r ] = Error< Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain > >::apply( a[ r ], b[ r ] );
    return error.infinity_norm();
  }
};



// random
// ------

template< class FieldType >
FieldType random ( FieldType a, FieldType b )
{
  FieldType delta = FieldType( rand() )/FieldType( RAND_MAX );
  return ( a + delta*( b - a ) );
}



// evaluateAll
// -----------

template< class Evaluate,
          class BasisFunctionSet, class VectorialBasisFunctionSet,
          class Point, class DofVector >
typename Evaluate::FunctionSpaceType::RangeFieldType
  evaluateAll ( const BasisFunctionSet &basisFunctionSet,
                const VectorialBasisFunctionSet &vectorialBasisFunctionSet,
                const Point &x, const DofVector &dofs )
{
  typedef typename Evaluate::FunctionSpaceType FunctionSpaceType;

  typedef typename VectorialBasisFunctionSet::DofAlignmentType DofAlignmentType;
  DofAlignmentType dofAlignment = vectorialBasisFunctionSet.dofAlignment();

  const std::size_t size = dofs.size();
  std::vector< typename FunctionSpaceType::RangeFieldType > convertedDofs( size );
  for( std::size_t i = 0; i < size; ++i )
  {
    const std::size_t j = DofConversion< DofAlignmentType >::apply( dofAlignment, i );
    convertedDofs.at( j ) = dofs[ i ];
  }

  typename Evaluate::Vector a, b;
  Evaluate::apply( basisFunctionSet, x, dofs, a );
  Evaluate::apply( vectorialBasisFunctionSet, x, convertedDofs, b );

  return Error< typename Evaluate::Vector >::apply( a, b );
}



// evaluateAll
// -----------

template< class Evaluate,
          class BasisFunctionSet, class VectorialBasisFunctionSet,
          class Point >
typename Evaluate::FunctionSpaceType::RangeFieldType
  evaluateAll ( const BasisFunctionSet &basisFunctionSet,
                const VectorialBasisFunctionSet &vectorialBasisFunctionSet,
                const Point &x )
{
  typedef typename Evaluate::Vector Vector;
  const std::size_t size = basisFunctionSet.size();
  std::vector< Vector > a( size ), b( size );

  Evaluate::apply( basisFunctionSet, x, a );
  Evaluate::apply( vectorialBasisFunctionSet, x, b );

  typedef typename VectorialBasisFunctionSet::DofAlignmentType DofAlignmentType;
  DofAlignmentType dofAlignment = vectorialBasisFunctionSet.dofAlignment();

  typedef typename Evaluate::FunctionSpaceType::RangeFieldType RangeFieldType;
  RangeFieldType error( 0 );
  for( std::size_t i = 0; i < size; ++i )
  {
    std::size_t j = DofConversion< DofAlignmentType >::apply( dofAlignment, i );
    error = std::max( error, Error< Vector >::apply( a[ i ], b[ j ] ) );
  }
  return error;
}


// compareQuadrature
// -----------------

template< class BasisFunctionSetType, class VectorialBasisFunctionSetType, class QuadratureType >
Dune::FieldVector< typename BasisFunctionSetType::RangeType::value_type, 3 >
  compareQuadrature ( const BasisFunctionSetType &basisFunctionSet,
                      const VectorialBasisFunctionSetType &vectorialBasisFunctionSet,
                      const QuadratureType &quadrature )
{
  // check order of basis function set
  if( vectorialBasisFunctionSet.order() != basisFunctionSet.order() )
    DUNE_THROW( Dune::InvalidStateException, "order() does not match" );

  // check size of basis function set
  const std::size_t size = basisFunctionSet.size();
  if( vectorialBasisFunctionSet.size() != basisFunctionSet.size() )
    DUNE_THROW( Dune::InvalidStateException, "size() does not match" );

  // get types
  typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  // init random dof vectors
  std::vector< RangeFieldType > dofs( size );
  for( std::size_t i = 0; i < size; ++i )
    dofs[ i ] = random< RangeFieldType >( 1e-3, 1e3 );

  // return value
  Dune::FieldVector< RangeFieldType, 3 > ret;

  const std::size_t nop = quadrature.nop();
  for( std::size_t qp = 0; qp < nop; ++qp )
  {
    // check evaluate methods
    {
      RangeFieldType error = 0;
      error += evaluateAll< Evaluate< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ], dofs );
      error += evaluateAll< Evaluate< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ] );
      ret[ 0 ] = std::max( ret[ 0 ], error );
    }

    // check jacobian methods
    {
      RangeFieldType error = 0;
      error += evaluateAll< Jacobian< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ], dofs );
      error += evaluateAll< Jacobian< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ] );
      ret[ 1 ] = std::max( ret[ 1 ], error );
    }

    // check hessian methods
    {
      RangeFieldType error = 0;
      error += evaluateAll< Hessian< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ], dofs );
      error += evaluateAll< Hessian< FunctionSpaceType > >( basisFunctionSet, vectorialBasisFunctionSet, quadrature[ qp ] );
      ret[ 2 ] = std::max( ret[ 1 ], error );
    }
  }

  return ret;
}


template< class GridPartType >
void traverse ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimRange = DIMRANGE;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, dimRange > FunctionSpaceType;

  // create discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType space( gridPart );

  // create scalar discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< typename FunctionSpaceType::ScalarFunctionSpaceType, GridPartType, POLORDER > ScalarDiscreteFunctionSpaceType;
  ScalarDiscreteFunctionSpaceType scalarSpace( gridPart );

  typedef Dune::FieldVector< typename FunctionSpaceType::RangeFieldType, 3 > ErrorType;
  ErrorType error( 0 );

  for( const auto& entity : space )
  {
    // get basis function set
    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    const BasisFunctionSetType basisFunctionSet = space.basisFunctionSet( entity );

    // get scalar basis function set
    typedef typename ScalarDiscreteFunctionSpaceType::BasisFunctionSetType ScalarBasisFunctionSetType;
    const ScalarBasisFunctionSetType scalarBasisFunctionSet = scalarSpace.basisFunctionSet( entity );

    // create vectorial basis function set
#if USE_VERTICAL_DOF_ALIGNMENT == 1
    typedef Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType, typename BasisFunctionSetType::RangeType, Dune::Fem::VerticalDofAlignment > VectorialBasisFunctionSetType;
#else
    typedef Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType, typename BasisFunctionSetType::RangeType, Dune::Fem::HorizontalDofAlignment > VectorialBasisFunctionSetType;
#endif
    VectorialBasisFunctionSetType vectorialBasisFunctionSet( scalarBasisFunctionSet );

    // create quadrature
    assert( basisFunctionSet.order() == vectorialBasisFunctionSet.order() );
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
    QuadratureType quadrature( basisFunctionSet.entity(), basisFunctionSet.order() );

    // compare basis function sets in quadrature points
    ErrorType localError = compareQuadrature( basisFunctionSet, vectorialBasisFunctionSet, quadrature );

    // update error
    for( int i = 0; i < 3; ++i )
      error[ i ] = std::max( error[ i ], localError[ i ] );
  }


  const double eps = 1e-8;
  // print error
  if( error.two_norm() > eps )
  {
    std::cerr << "Errors( evaluateAll, jacobianAll, hessianAll ): " << error << std::endl;
#if USE_VERTICAL_DOF_ALIGNMENT == 1
    DUNE_THROW( Dune::InvalidStateException, "VectorialBasisFunctionSet< VerticalDofAlignment > test failed." );
#else
    DUNE_THROW( Dune::InvalidStateException, "VectorialBasisFunctionSet< HorizontalDofAlignment > test failed." );
#endif
  }
}


int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( argc >= 2 ? argv[ 1 ] : "parameter" );

  srand( time( NULL ) );

  typedef Dune::GridSelector::GridType GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();

  int refCount = Dune::Fem::Parameter::getValue< int >( "startLevel", 0 );
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( grid, refCount*refineStepsForHalf );

  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );
  traverse( gridPart );

  return 0;
}
