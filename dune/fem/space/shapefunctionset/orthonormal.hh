#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_ORTHONORMAL_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_ORTHONORMAL_HH

#include <cassert>
#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/orthonormal/orthonormalbase_1d.hh>
#include <dune/fem/space/shapefunctionset/orthonormal/orthonormalbase_2d.hh>
#include <dune/fem/space/shapefunctionset/orthonormal/orthonormalbase_3d.hh>

namespace Dune
{

  namespace Fem
  {

#ifndef DOXYGEN

    // OrthonormalShapeFunctions
    // -------------------------

    template< int dimension >
    struct OrthonormalShapeFunctions;

    template<>
    struct OrthonormalShapeFunctions< 1 >
    {
      static constexpr std::size_t size ( int order )
      {
        return static_cast< std::size_t >( order+1 );
      }
    };

    template<>
    struct OrthonormalShapeFunctions< 2 >
    {
      static constexpr std::size_t size ( int order )
      {
        return static_cast< std::size_t >( (order+2)*(order+1)/2 );
      }
    };

    template<>
    struct OrthonormalShapeFunctions< 3 >
    {
      static constexpr std::size_t size ( int order )
      {
        return static_cast< std::size_t >( ((order+1)*(order+2)*(2*order+3)/6+(order+1)*(order+2)/2)/2 );
      }
    };

#endif // #ifndef DOXYGEN



    // OrthonormalShapeFunctionSet
    // ---------------------------

    template< class FunctionSpace >
    class OrthonormalShapeFunctionSet
    {
      typedef OrthonormalShapeFunctionSet< FunctionSpace > ThisType;

      static_assert( FunctionSpace::dimDomain <= 3, "OrthonormalShapeFunctionSet only implemented up to domain dimension 3" );
      static_assert( FunctionSpace::dimRange == 1, "OrthonormalShapeFunctionSet only implemented for scalar function spaces" );

    public:
      /** \copydoc Dune::Fem::ShapeFunctionSet::FunctionSpaceType */
      typedef FunctionSpace FunctionSpaceType;

    private:
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef RangeFieldType (*Evaluate) ( const int, const DomainFieldType * );
      typedef void (*Jacobian) ( const int i, const DomainFieldType *, RangeFieldType * );

      // hessian is only available for 2d-simplices (aka triangles)
      typedef RangeFieldType Array[ 3 ];
      typedef void (*Hessian) ( const int i, const DomainFieldType *, Array & );

      typedef OrthonormalBase_1D< DomainFieldType, RangeFieldType > OrthonormalBase1d;
      typedef OrthonormalBase_2D< DomainFieldType, RangeFieldType > OrthonormalBase2d;
      typedef OrthonormalBase_3D< DomainFieldType, RangeFieldType > OrthonormalBase3d;

      static void setFunctionPointers(const Dune::GeometryType& geomType,
                                      Evaluate &evaluate, Jacobian &jacobian )
      {
        if( geomType.isLine() )
        {
          evaluate = OrthonormalBase1d::eval_line;
          jacobian = OrthonormalBase1d::grad_line;
          return ;
        }
        else if( geomType.isQuadrilateral() )
        {
          evaluate = OrthonormalBase2d::eval_quadrilateral_2d;
          jacobian = OrthonormalBase2d::grad_quadrilateral_2d;
          return ;
        }
        else if( geomType.isTriangle() )
        {
          evaluate = OrthonormalBase2d::eval_triangle_2d;
          jacobian = OrthonormalBase2d::grad_triangle_2d;
          return ;
        }
        else if( geomType.isHexahedron() )
        {
          evaluate = OrthonormalBase3d::eval_hexahedron_3d;
          jacobian = OrthonormalBase3d::grad_hexahedron_3d;
          return ;
        }
        else if ( geomType.isTetrahedron() )
        {
          evaluate = OrthonormalBase3d::eval_tetrahedron_3d;
          jacobian = OrthonormalBase3d::grad_tetrahedron_3d;
          return ;
        }
        else if( geomType.isPrism() )
        {
          evaluate = OrthonormalBase3d::eval_prism_3d;
          jacobian = OrthonormalBase3d::grad_prism_3d;
          return ;
        }
        else if ( geomType.isPyramid() )
        {
          evaluate = OrthonormalBase3d::eval_pyramid_3d;
          jacobian = OrthonormalBase3d::grad_pyramid_3d;
          return ;
        }

        DUNE_THROW(InvalidStateException,"Invalid geometry type " << geomType );
      }

    public:
      /** \copydoc Dune::Fem::ShapeFunctionSet::dimension */
      static const int dimension = FunctionSpaceType::dimDomain;

      /** \copydoc Dune::Fem::ShapeFunctionSet::Domaintype */
      typedef typename FunctionSpaceType::DomainType DomainType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::RangeType */
      typedef typename FunctionSpaceType::RangeType RangeType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::JacobianRangeType */
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::HessianRangeType */
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      /** \name Construction
       *  \{
       */

      OrthonormalShapeFunctionSet () = default;

      OrthonormalShapeFunctionSet ( GeometryType type, int order )
        : order_( order ),
          evaluate_( nullptr ),
          jacobian_( nullptr ),
          hessian_( nullptr )
      {
        // for type none simply create cube basis function set.
        if( type.isNone() )
          type = Dune::GeometryTypes::cube(type.dim());

        // set functions pointers for evaluate and jacobian
        // depending on geometry type
        setFunctionPointers( type, evaluate_, jacobian_ );
        assert( evaluate_ );
        assert( jacobian_ );

        if( dimension == 2 )
        {
          if( type.isTriangle() )
            hessian_ = OrthonormalBase2d::hess_triangle_2d;
          else if( type.isQuadrilateral() )
            hessian_ = OrthonormalBase2d::hess_quadrilateral_2d;
        }
      }

      OrthonormalShapeFunctionSet ( const ThisType & ) = default;

      OrthonormalShapeFunctionSet ( ThisType && ) = default;

      /** \} */

      /** \name Assignment
       *  \{
       */

      OrthonormalShapeFunctionSet &operator= ( const ThisType & ) = default;

      OrthonormalShapeFunctionSet &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::ShapeFunctionSet::order */
      int order () const { return order_; }

      /** \copydoc Dune::Fem::ShapeFunctionSet::size */
      std::size_t constexpr size () const { return size( order() ); }

      /** \brief please doc me */
      static std::size_t constexpr size ( int order )
      {
        return OrthonormalShapeFunctions< dimension >::size( order );
      }

      /** \copydoc Dune::Fem::ShapeFunctionSet::evaluateEach */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        const DomainType y = Dune::Fem::coordinate( x );
        RangeType value;
        const std::size_t size = this->size();
        for( std::size_t i = 0; i < size; ++i )
        {
          evaluate( i, y, value );
          functor( i, value );
        }
      }

      /** \copydoc Dune::Fem::ShapeFunctionSet::jacobianEach */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        const DomainType y = Dune::Fem::coordinate( x );
        JacobianRangeType jacobian;
        const std::size_t size = this->size();
        for( std::size_t i = 0; i < size; ++i )
        {
          this->jacobian( i, y, jacobian );
          functor( i, jacobian );
        }
      }

      /** \copydoc Dune::Fem::ShapeFunctionSet::hessianEach */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        const DomainType y = Dune::Fem::coordinate( x );
        HessianRangeType hessian;
        const std::size_t size = this->size();
        for( std::size_t i = 0; i < size; ++i )
        {
          this->hessian( i, y, hessian );
          functor( i, hessian );
        }
      }

      /** \} */

    private:
      void evaluate ( std::size_t i, const DomainType &x, RangeType &value ) const
      {
        value[ 0 ] = evaluate_( i, &x[ 0 ] );
      }

      void jacobian ( std::size_t i, const DomainType &x, JacobianRangeType &jacobian ) const
      {
        jacobian_( i, &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }

      void hessian ( std::size_t i, const DomainType &x, HessianRangeType &hessian ) const
      {
        assert( hessian_ );

        RangeFieldType values[] = { 0, 0, 0 };
        hessian_( i , &x[ 0 ], values );
        for( unsigned int j = 0; j < FunctionSpaceType::dimDomain;  ++j )
          for( unsigned int k = 0; k < FunctionSpaceType::dimDomain; ++k )
            hessian[ 0 ][ j ][ k ] = values[ j + k ];
      }

      int order_;
      Evaluate evaluate_;
      Jacobian jacobian_;
      Hessian hessian_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_ORTHONORMAL_HH
