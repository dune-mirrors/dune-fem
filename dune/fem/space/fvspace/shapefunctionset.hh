#ifndef DUNE_FEM_SPACE_FVSPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_FVSPACE_SHAPEFUNCTIONSET_HH

// dune-common includes 
#include <dune/common/static_assert.hh>

// dune-fem includes 
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

/*
  \file
  \author Christoph Gersbacher
  \brief Shape function set for Finite Volume space
 */


namespace Dune
{

  namespace Fem 
  {

    // FVScalarShapeFunctionSet
    // ------------------------

    template< class FunctionSpace, int polOrder >
    struct FVScalarShapeFunctionSet;



    // Template specialization for polynomial order 0
    // ----------------------------------------------

    template< class FunctionSpace >
    class FVScalarShapeFunctionSet< FunctionSpace, 0 >
    : public ShapeFunctionSet< FunctionSpace, FVScalarShapeFunctionSet< FunctionSpace, 0 > >
    {
      typedef FVScalarShapeFunctionSet< FunctionSpace, 0 > ThisType;
      typedef ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    public:
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianType JacobianRangeType;
      typedef typename BaseType::HessianType HessianRangeType;

    private:
      typedef typename RangeType::value_type RangeFieldType;
      static const int dimRange = RangeFieldType::dimension;

    public:
      FVScalarShapeFunctionSet ( const GeometryType &type )
      : type_( type )
      {}

      GeometryType type () const { return type_; }

      std::size_t size () const { return 1; }

      template< class Point, class Functor >
      void evaluateEach ( const Point &, Functor functor ) const
      {
        functor( 0, RangeType( 1 ) );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &, Functor functor ) const
      {
        functor( 0, JacobianRangeType( 0 ) );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &, Functor functor ) const
      {
        functor( 0, HessianRangeType( 0 ) );
      }

    private:
      GeometryType type_;
    };



    // FVShapeFunctionSet
    // ------------------

    /*
     * \brief Implementation of Dune::Fem::ShapeFunctionSet for Finite Volume spaces 
     *
     * \tparam  FunctionSpace  Function space
     * \tparam  polOrder  Polynomial order
     *
     * \note This shape function set is only implemented for polynomial order 0.
     */
    template< class FunctionSpace, int polOrder >
    struct FVShapeFunctionSet
    : public ShapeFunctionSet< FunctionSpace, FVShapeFunctionSet< FunctionSpace, polOrder > >
    {
      dune_static_assert( (polOrder == 0),
                           "FVShapeFunctionSet only availabe for polynomial order 0." );

      typedef FVShapeFunctionSet< FunctionSpace, polOrder > ThisType;
      typedef ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

    public:
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::RangeType RangeType;

    protected:
      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef FVScalarShapeFunctionSet< ScalarFunctionSpaceType, polOrder > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > ImplementationType;

    public:
      FVShapeFunctionSet( const GeometryType &type )
      : implementation_( ScalarShapeFunctionSetType( type ) )
      {}

      GeometryType type () const { return implementation().type(); }

      std::size_t size () const { return implementation().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        implementation().evaluateEach( x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point & x, Functor functor ) const
      {
        implementation().jacobianEach( x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point & x, Functor functor ) const
      {
        implementation().hessianEach( x, functor );
      }

    protected:
      const ImplementationType &implementation () const
      {
        return implementation_;
      }

    private:
      ImplementationType implementation_;
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FVSPACE_SHAPEFUNCTIONSET_HH
