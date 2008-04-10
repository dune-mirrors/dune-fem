#ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
#define DUNE_FEM_LAGRANGEINTERPOLATION_HH

#include <dune/common/typetraits.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>

namespace Dune 
{

  /** \class LagrangeInterpolation
   *  \brief Generates the Lagrange Interpolation of an analytic function
   */
  template< class DiscreteFunctionImp >
  class LagrangeInterpolation
  {
  public:
    //! type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of local functions
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
      LagrangePointSetType;
    //! type of vectors in function's domain
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of vectors in function's range
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  public:
    /** interpolate an analytical function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  globally) at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class FunctionType >
    inline static void interpolateFunction ( const FunctionType &function,
                                             DiscreteFunctionType &discreteFunction );

  private:
    template< class FunctionType, bool hasLocalFunction >
    struct callInterpolateDiscreteFunction;

  protected:
    /** interpolate a discrete function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  locally at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class FunctionType >
    inline static void interpolateDiscreteFunction
      ( const FunctionType &function,
        DiscreteFunctionType &discreteFunction );
  };

}

#include "lagrangeinterpolation_inline.hh"

#endif
