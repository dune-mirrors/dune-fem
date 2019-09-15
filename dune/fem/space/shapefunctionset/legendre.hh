#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/space/shapefunctionset/legendrepolynomials.hh>

namespace Dune
{

  namespace Fem
  {

    // LegendreShapeFunction
    // ---------------------

    /** \brief implementation of a single scalar-valued Legendre shape function
     *
     *  \note The range field type used in the evaluation is fixed to `double`.
     *
     *  \tparam  FunctionSpace  (scalar) function space
     */
    template< class FunctionSpace >
    class LegendreShapeFunction
    {
      static_assert( FunctionSpace::dimRange == 1, "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    public:
      /** \copydoc Dune::Fem::Function::FunctionSpaceType */
      typedef FunctionSpace FunctionSpaceType;

      /** \copydoc Dune::Fem::Function::DomainFieldType */
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      /** \copydoc Dune::Fem::Function::RangeFieldType */
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      /** \copydoc Dune::Fem::Function::DomainType */
      typedef typename FunctionSpaceType::DomainType DomainType;
      /** \copydoc Dune::Fem::Function::RangeType */
      typedef typename FunctionSpaceType::RangeType RangeType;
      /** \copydoc Dune::Fem::Function::JacobianRangeType */
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      /** \copydoc Dune::Fem::Function::HessianRangeType */
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    private:
      static const int dimDomain = FunctionSpaceType::dimDomain;

    public:
      /** \name Construction
       *  \{
       */

      LegendreShapeFunction () = default;

      template< class MultiIndex >
      explicit LegendreShapeFunction ( const MultiIndex &multiIndex )
      {
        using std::begin;
        using std::end;

        auto first = begin( multiIndex );
        auto last = end( multiIndex );
        assert( std::distance( first, last ) == dimDomain );
        std::copy( first, last, multiIndex_.begin() );
      }

      /** \} */

      /** \brief return polynomial order of this function */
      int order () const noexcept
      {
        return *std::max_element( multiIndex_.begin(), multiIndex_.end() );
      }

      /** \brief return monomial orders of this function */
      const std::array< int, FunctionSpaceType::dimDomain > &orders () const noexcept
      {
        return multiIndex_;
      }

      /** \copydoc Dune::Fem::Function::evaluate */
      void evaluate ( const DomainType &x, RangeType &value ) const noexcept
      {
        value[ 0 ] = RangeFieldType( 1 );
        for( int i = 0; i < dimDomain; ++i )
          value[ 0 ] *= LegendrePolynomials::evaluate( multiIndex_[ i ], x[ i ] );
      }

      /** \copydoc Dune::Fem::Function::jacobian */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const noexcept
      {
        jacobian = JacobianRangeType( 1 );
        for( int k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = LegendrePolynomials::evaluate( multiIndex_[ k ], x[ k ] );
          const RangeFieldType dphi = LegendrePolynomials::jacobian( multiIndex_[ k ], x[ k ]);
          for( int i = 0; i < dimDomain; ++i )
            jacobian[ 0 ][ i ] *= ( k == i ) ? dphi : phi;
        }
      }

      /** \copydoc Dune::Fem::Function::hessian */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const noexcept
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 1 ) );
        for( int k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = LegendrePolynomials::evaluate( multiIndex_[ k ], x[ k ] );
          const RangeFieldType dphi = LegendrePolynomials::jacobian( multiIndex_[ k ], x[ k ] );
          for( int i = 0; i < dimDomain; ++i )
          {
            hessian[ 0 ][ i ][ i ] *= ( k == i ) ? LegendrePolynomials::hessian( multiIndex_[ i ], x[ i ]) : phi;
            for( int j = i+1; j < dimDomain; ++j )
            {
              RangeFieldType tmp = ( k == i || k == j ) ? dphi : phi;
              hessian[ 0 ][ i ][ j ] *= tmp;
              hessian[ 0 ][ j ][ i ] *= tmp;
            }
          }
        }
      }

    private:
      std::array< int, dimDomain > multiIndex_;
    };



#ifndef DOXYGEN

    namespace __LegendreShapeFunctionSet
    {

      // DefaultFactory
      // --------------

      template< class FunctionSpace >
      class DefaultFactory
      {
        typedef LegendreShapeFunction< FunctionSpace > ShapeFunctionType;

        static const int dimDomain = FunctionSpace::dimDomain;
        typedef std::array< int, dimDomain > MultiIndex;

      public:
        explicit DefaultFactory ( int order )
          : order_( order )
        {}

        int order () const noexcept { return order_; }

        std::size_t size () const noexcept
        {
          std::size_t size = 1;
          for( int i = 0; i < dimDomain; ++i )
            size *= order() + 1;
          return size;
        }

        template< class InputIterator >
        void operator() ( InputIterator first ) const noexcept
        {
          auto function = [&first]( const MultiIndex &multiIndex )
          {
            *first = ShapeFunctionType( multiIndex );
            ++first;
          };

          MultiIndex multiIndex;
          fill( multiIndex, function, &multiIndex[ 0 ], dimDomain, order() );
        }

      private:
        template< class Function >
        static void fill ( MultiIndex &multiIndex, Function function,
                           int *begin, std::size_t n, int order )
        {
          if( n > 0u )
          {
            for( *begin = 0; *begin <= order; *begin += 1 )
              fill( multiIndex, function, begin+1, n-1, order );
          }
          else
            function( multiIndex );
        }

        int order_;
      };

    } // namespace __LegendreShapeFunctionSet

#endif // #ifndef DOXYGEN



    // LegendreShapeFunctionSet
    // ------------------------

    /** \brief a Dune::Fem::ShapeFunctionSet of Legendre ansatz polynomials
     *
     *  \note The range field type used in the evaluation is fixed to `double`.

     *  \note This shape function set can only be used with cubic reference elements.
     *
     *  \tparam  FunctionSpace        (scalar) function space
     *  \tparam  hierarchicalOrdering (bool) if true shape functions are ordered according to their polynomial order
     */
    template< class FunctionSpace, bool hierarchicalOrdering = false >
    class LegendreShapeFunctionSet
    {
    protected:
      typedef LegendreShapeFunction< FunctionSpace > ShapeFunctionType;

    public:
      /** \copydoc Dune::Fem::ShapeFunctionSet::FunctionSpaceType */
      typedef FunctionSpace FunctionSpaceType;

      /** \copydoc Dune::Fem::ShapeFunctionSet::DomainType */
      typedef typename FunctionSpaceType::DomainType DomainType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::RangeType */
      typedef typename FunctionSpaceType::RangeType RangeType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::JacobianRangeType */
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      /** \copydoc Dune::Fem::ShapeFunctionSet::HessianRangeType */
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    protected:
      struct Compare
      {
        bool operator() ( const ShapeFunctionType &lhs, const ShapeFunctionType &rhs ) noexcept
        {
          if( lhs.order() != rhs.order() )
            return lhs.order() < rhs.order();
          else
          {
            const auto &a = lhs.orders();
            const auto &b = rhs.orders();
            return std::lexicographical_compare( a.begin(), a.end(), b.begin(), b.end() );
          }
        }
      };

    public:
      /** \name Construction
       *  \{
       */

      /** \brief default constructor resulting in uninitialized shape function
       *         set */
      LegendreShapeFunctionSet () = default;

      /** \brief initialize with polynomial order
       *
       *  \param[in]  order
       */
      explicit LegendreShapeFunctionSet ( int order )
        : LegendreShapeFunctionSet( __LegendreShapeFunctionSet::DefaultFactory< FunctionSpaceType >( order ) )
      {}

      /** \brief initialize from user-defined factory object
       *
       *  \param[in]  factory  a factory, see description
       *
       *  \note The parameter factory must implement the following methods:
       *  \code
       *  struct Factory
       *  {
       *    // return number of shape functions
       *    std::size_t size () const noexcept;
       *
       *    // return maximum order
       *    int order () const noexcept;
       *
       *    // fill range from begin to begin+size() with unique shape functions
       *    void operator() ( InputIterator begin ) const noexcept;
       *  };
       *  \endcode
       */
      template< class Factory >
      LegendreShapeFunctionSet ( const Factory &factory )
        : shapeFunctions_( factory.size() ),
          order_( factory.order() )
      {
        factory( shapeFunctions_.begin() );

        // if hierarchical ordering was selected sort shape functions
        // according to polynomial order
        if( hierarchicalOrdering )
        {
          std::sort( shapeFunctions_.begin(), shapeFunctions_.end(), Compare() );
        }
      }

      /** \} */

      /** \copydoc Dune::Fem::ShapeFunctionSet::order */
      int order () const noexcept { return order_; }

      /** \copydoc Dune::Fem::ShapeFunctionSet::size */
      std::size_t size () const noexcept { return shapeFunctions_.size(); }

      /** \copydoc Dune::Fem::ShapeFunctionSet::evaluateEach */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const noexcept
      {
        RangeType value;
        std::size_t i = 0u;
        for( const ShapeFunctionType &shapeFunction : shapeFunctions_ )
        {
          shapeFunction.evaluate( coordinate( x ), value );
          functor( i++, value );
        }
      }

      /** \copydoc Dune::Fem::ShapeFunctionSet::jacobianEach */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const noexcept
      {
        JacobianRangeType jacobian;
        std::size_t i = 0u;
        for( const ShapeFunctionType &shapeFunction : shapeFunctions_ )
        {
          shapeFunction.jacobian( coordinate( x ), jacobian );
          functor( i++, jacobian );
        }
      }

      /** \copydoc Dune::Fem::ShapeFunctionSet::hessianEach */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const noexcept
      {
        HessianRangeType hessian;
        std::size_t i = 0u;
        for( const ShapeFunctionType &shapeFunction : shapeFunctions_ )
        {
          shapeFunction.hessian( coordinate( x ), hessian );
          functor( i++, hessian );
        }
      }

    protected:
      std::vector< ShapeFunctionType > shapeFunctions_;
      int order_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDRE_HH
