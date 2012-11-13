#include <config.h>

// C++ includes
#include <vector>

// dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/exceptions.hh>

// dune-geometry includes
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/anisotropicdgspace/multiindexset.hh>
#include <dune/fem/space/anisotropicdgspace/shapefunctionset.hh>
#include <dune/fem/quadrature/quadrature.hh>



// CheckOrthonormalShapeFunctionSet
// --------------------------------

template< class ShapeFunctionSet >
class CheckOrthonormalShapeFunctionSet
{
  typedef CheckOrthonormalShapeFunctionSet< ShapeFunctionSet > ThisType;


public:
  // shape function set type
  typedef ShapeFunctionSet ShapeFunctionSetType;

private:
  // function space type
  typedef typename ShapeFunctionSetType::FunctionSpaceType FunctionSpaceType;
  // range type
  typedef typename FunctionSpaceType::RangeType RangeType;
  // range field type
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  // dimension
  static const int dimension = FunctionSpaceType::dimDomain;
  // quadrature type
  typedef Dune::Fem::Quadrature< RangeFieldType, dimension > QuadratureType;

  // Evaluate
  // --------

  struct Evaluate
  {
    Evaluate ( std::vector< RangeType > &phi )
    : phi_( phi )
    {}

    void operator() ( const int shapeFunction, const RangeType &value )
    {
      phi_[ shapeFunction ] = value;
    }

  private:
    std::vector< RangeType > &phi_;
  };

public:
  // mass matrix type
  typedef Dune::DynamicMatrix< RangeFieldType > MassMatrixType;

  CheckOrthonormalShapeFunctionSet ( const Dune::GeometryType &type,
                                     const ShapeFunctionSetType &shapeFunctionSet,
                                     int order )
  : shapeFunctionSet_( shapeFunctionSet )
  {
    QuadratureType quad( type, order );
    setupMassMatrix( quad );
    std::cout << massMatrix_ << std::endl;
  }

  // return mass matrix
  const MassMatrixType &massMatrix () const { return massMatrix_; }

private:
  void setupMassMatrix ( const QuadratureType &quad )
  {
    // get size of shape function set
    const std::size_t size = shapeFunctionSet().size();

    // resize mass matrix
    massMatrix_.resize( size, size, RangeFieldType( 0 ) );
    // create storage
    std::vector< RangeType > phi( size );

    const int nop = quad.nop();
    for( int qp = 0; qp < nop; ++qp )
    {
      // get quadrature weight
      const RangeFieldType weight = quad.weight(qp);

      // evaluate shape function set
      Evaluate functor( phi );
      shapeFunctionSet().evaluateEach( quad[ qp ], functor );

      // set up mass matrix
      for( std::size_t i = 0; i < size; ++i )
      {
        const RangeType &phi_i = phi[ i ];
        const RangeFieldType value = weight * ( phi_i * phi_i );
        massMatrix_[ i ][ i ] += value;

        for( std::size_t j = i+1; j < size; ++j )
        {
          const RangeFieldType value = weight * ( phi_i * phi[ j ] );
          massMatrix_[ i ][ j ] += value;
          massMatrix_[ j ][ i ] += value;
        }
      }
    }
  }

private:
  // return shape function set
  const ShapeFunctionSetType shapeFunctionSet () const
  {
    return shapeFunctionSet_;
  }

  ShapeFunctionSetType shapeFunctionSet_;
  MassMatrixType massMatrix_;
};



// ShapeFunctionSet
// ----------------

template< int dimension, int maxOrder >
class ShapeFunctionSet
{
  typedef Dune::Fem::FunctionSpace< double, double, dimension, 1 > FunctionSpaceType;

public:
  typedef AnisotropicDG::ShapeFunctionSet< FunctionSpaceType, maxOrder > Type;
};



int main ( int argc, char **argv )
try
{
  const int dimension = GRIDDIM;
  const int maxOrder = POLORDER;

  // get geometry type
  Dune::GenericGeometry::CubeTopology< dimension >::type topology; 
  Dune::GeometryType type( topology );

  // shape function set type
  typedef ShapeFunctionSet< dimension, maxOrder >::Type ShapeFunctionSetType;

  typedef AnisotropicDG::MultiIndexSet< dimension, maxOrder > MultiIndexSetType;
  typedef MultiIndexSetType::MultiIndexType MultiIndexType;
  typedef MultiIndexSetType::IteratorType IteratorType;

  // traverse all multi indices
  const IteratorType end = MultiIndexSetType::end();
  for( IteratorType it = MultiIndexSetType::begin(); it != end; ++it )
  {
    // get multi index
    const MultiIndexType &multiIndex = *it;

    // create shape function set
    ShapeFunctionSetType shapeFunctionSet( multiIndex );

    // desired order of quadrature rule
    int order = 2*maxOrder + 1;

    // check shape function set
    CheckOrthonormalShapeFunctionSet< ShapeFunctionSetType > check( type, shapeFunctionSet, order );
  }

  return 0;
}
catch( Dune::Exception &e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
catch( std::exception &e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
