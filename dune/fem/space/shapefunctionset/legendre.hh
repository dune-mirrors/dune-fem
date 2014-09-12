#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH

// C++ includes
#include <algorithm>
#include <cstddef>

// dune-common includes
#include <dune/common/array.hh>

// dune-fem includes
#include <dune/fem/space/shapefunctionset/legendrepolynomials.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Provides an implementation of Dune::Fem::ShapeFunctionSet consisting of Legendre polynomials
*/


namespace Dune 
{

  namespace Fem
  {

    // LegendreShapeFunction
    // ---------------------

    /**
     * \brief Class representing a single scalar-valued Legendre shape function
     *
     * \tparam  FunctionSpace  Scalar function space
     */
    template< class FunctionSpace >
    class LegendreShapeFunction
    {
      static_assert( (FunctionSpace::dimRange == 1),
                     "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      typedef LegendreShapeFunction< FunctionSpace > ThisType;

    public:
      typedef FunctionSpace FunctionSpaceType;

      static const int dimension = FunctionSpaceType::dimDomain;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef Dune::array< int, dimension > MultiIndexType;

      explicit LegendreShapeFunction ( const MultiIndexType &multiIndex )
      : multiIndex_( multiIndex )
      {}

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        value[ 0 ] = RangeFieldType( 1 );
        for( int i = 0; i < dimension; ++i )
          value[ 0 ] *= LegendrePolynomials::evaluate( multiIndex_[ i ], x[ i ] );
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 1 );
        for( int k = 0; k < dimension; ++k )
        {
          const RangeFieldType phi = LegendrePolynomials::evaluate( multiIndex_[ k ], x[ k ] );
          const RangeFieldType dphi = LegendrePolynomials::jacobian( multiIndex_[ k ], x[ k ]);
          for( int i = 0; i < dimension; ++i )
            jacobian[ 0 ][ i ] *= ( k == i ) ? dphi : phi;
        }
      }

      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 1 ) );
        for( int k = 0; k < dimension; ++k )
        {
          const RangeFieldType phi = LegendrePolynomials::evaluate( multiIndex_[ k ], x[ k ] );
          const RangeFieldType dphi = LegendrePolynomials::jacobian( multiIndex_[ k ], x[ k ] );
          for( int i = 0; i < dimension; ++i )
          {
            hessian[ i ][ i ] *= ( k == i ) ? LegendrePolynomials::hessian( multiIndex_[ i ], x[ i ]) : phi;
            for( int j = i+1; j < dimension; ++j )
            {
              RangeFieldType tmp = ( k == i || k == j ) ? dphi : phi;
              hessian[ i ][ j ] *= tmp;
              hessian[ j ][ i ] *= tmp;
            }
          }
        }
      }
      
      const ThisType *clone () const { return new ThisType( *this ); }

    private:
      MultiIndexType multiIndex_;
    };



    // LegendreShapeFunctionSetFactory
    // -------------------------------

    /**
     * \brief Factory to be used together with Dune::Fem::SimpleShapeFunctionSet
     *
     * \tparam  FunctionSpace  Scalar function space
     */
    template< class FunctionSpace >
    class LegendreShapeFunctionSetFactory
    {
      static const int dimension = FunctionSpace::dimDomain;

      static_assert( (FunctionSpace::dimRange == 1),
                     "FunctionSpace must be scalar (i.e., dimRange = 1)." );
    public:
      typedef FunctionSpace FunctionSpaceType;
      typedef LegendreShapeFunction< FunctionSpaceType > ShapeFunctionType;

      explicit LegendreShapeFunctionSetFactory ( int order )
      : order_( order )
      {}

      int order () const { return order_; }

      std::size_t numShapeFunctions () const
      {
        std::size_t size = 1;
        for( int i = 0; i < dimension; ++i )
          size *= (order()+1);
        return size;
      }

      ShapeFunctionType *createShapeFunction ( std::size_t shapeFunction ) const
      {
        typename ShapeFunctionType::MultiIndexType multiIndex;
        for( int i = 0; i < dimension; ++i )
        {
          multiIndex[ i ] = shapeFunction % (order()+1);
          shapeFunction /= (order()+1);
        }
        return new ShapeFunctionType( multiIndex );
      }

    private:
      int order_;
    };



    // LegendreShapeFunctionSet
    // ------------------------

    /**
     * \brief Implementation of Dune::Fem::ShapeFunctionSet based on Legendre polynomials
     *
     * \tparam  FunctionSpace  Scalar function space
     *
     * \note This shape function set can only be used with cubic reference elements
     */
    template< class FunctionSpace >
    class LegendreShapeFunctionSet
    : public SimpleShapeFunctionSet< typename LegendreShapeFunctionSetFactory< FunctionSpace >::ShapeFunctionType >
    {
      static_assert( (FunctionSpace::dimRange == 1),
                     "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      typedef LegendreShapeFunctionSetFactory< FunctionSpace > FactoryType;
      typedef typename LegendreShapeFunctionSetFactory< FunctionSpace >::ShapeFunctionType ShapeFunctionType;
      typedef SimpleShapeFunctionSet< ShapeFunctionType > BaseType;

    public:
      explicit LegendreShapeFunctionSet ( const int order )
      : BaseType( FactoryType( order ) )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH
