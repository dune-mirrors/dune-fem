#ifndef DUNE_FEM_SPACE_FiniteVolumeSPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_FiniteVolumeSPACE_SHAPEFUNCTIONSET_HH

// dune-common includes 
#include <dune/common/static_assert.hh>

// dune-fem includes 
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

    // FiniteVolumeScalarShapeFunctionSet
    // ----------------------------------

    template< class FunctionSpace, int polOrder >
    struct FiniteVolumeScalarShapeFunctionSet;



    // Template specialization for polynomial order 0
    // ----------------------------------------------

    template< class FunctionSpace >
    class FiniteVolumeScalarShapeFunctionSet< FunctionSpace, 0 >
    {
      typedef FiniteVolumeScalarShapeFunctionSet< FunctionSpace, 0 > ThisType;

      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    public:
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    private:
      typedef typename RangeType::value_type RangeFieldType;
      static const int dimRange = RangeFieldType::dimension;

    public:
      static std::size_t size () { return 1; }

      template< class Point, class Functor >
      static void evaluateEach ( const Point &, Functor functor )
      {
        functor( 0, RangeType( 1 ) );
      }

      template< class Point, class Functor >
      static void jacobianEach ( const Point &, Functor functor )
      {
        functor( 0, JacobianRangeType( 0 ) );
      }

      template< class Point, class Functor >
      static void hessianEach ( const Point &, Functor functor )
      {
        functor( 0, HessianRangeType( 0 ) );
      }
    };



    // FiniteVolumeShapeFunctionSet
    // ----------------------------

    /*
     * \brief Implementation of Dune::Fem::ShapeFunctionSet for Finite Volume spaces 
     *
     * \tparam  FunctionSpace  Function space
     * \tparam  polOrder  Polynomial order
     *
     * \note This shape function set is only implemented for polynomial order 0.
     */
    template< class FunctionSpace, int polOrder >
    struct FiniteVolumeShapeFunctionSet
    {
      dune_static_assert( (polOrder == 0),
                           "FiniteVolumeShapeFunctionSet only availabe for polynomial order 0." );

      typedef FiniteVolumeShapeFunctionSet< FunctionSpace, polOrder > ThisType;

    public:
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::RangeType JacobianRangeType;
      typedef typename FunctionSpaceType::RangeType HessianRangeType;

    protected:
      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef FiniteVolumeScalarShapeFunctionSet< ScalarFunctionSpaceType, polOrder > ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > VectorialShapeFunctionSetType;

    public:
      FiniteVolumeShapeFunctionSet () {}

      std::size_t size () const { return vectorialShapeFunctionSet().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        vectorialShapeFunctionSet().evaluateEach( x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point & x, Functor functor ) const
      {
        vectorialShapeFunctionSet().jacobianEach( x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point & x, Functor functor ) const
      {
        vectorialShapeFunctionSet().hessianEach( x, functor );
      }

    protected:
      const VectorialShapeFunctionSetType &vectorialShapeFunctionSet () const
      {
        return vectorialShapeFunctionSet_;
      }

    private:
      VectorialShapeFunctionSetType vectorialShapeFunctionSet_;
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FiniteVolumeSPACE_SHAPEFUNCTIONSET_HH
