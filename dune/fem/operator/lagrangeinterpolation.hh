#ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
#define DUNE_FEM_LAGRANGEINTERPOLATION_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class LagrangeInterpolation
     *  \brief Generates the Lagrange Interpolation of an analytic function
     */
    template< class Function, class DiscreteFunction >
    class LagrangeInterpolation : public Operator< Function, DiscreteFunction >
    {
      typedef LagrangeInterpolation< Function, DiscreteFunction > ThisType;

    public:
      //! type of discrete functions
      typedef DiscreteFunction DiscreteFunctionType;

      //! type of discrete function space
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      //! type of local functions
      typedef typename DiscreteFunctionType::LocalFunctionType
        LocalFunctionType;

      //! type of grid partition
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      //! type of Lagrange point set
      typedef typename DiscreteFunctionSpaceType::LagrangePointSetType
        LagrangePointSetType;
      //! type of vectors in function's domain
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      //! type of vectors in function's range
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    public:
      //! empty contructor
      LagrangeInterpolation () {}

      //! virtual destructor because of inheritance from Operator
      virtual ~LagrangeInterpolation () {}

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
      void operator () ( const Function &function,
                         DiscreteFunctionType &discreteFunction ) const
      {
        interpolateFunction( function, discreteFunction );
      }

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
      static void interpolateFunction ( const Function &function, DiscreteFunctionType &discreteFunction )
      {
        const bool hasLocalFunction = Conversion< Function, HasLocalFunction >::exists;
        interpolateFunction( function, discreteFunction, integral_constant< bool, hasLocalFunction >() );
      }

      //! \copydoc interpolateFunction
      //- make interface equal to DGL2Projection
      static void apply ( const Function &function, DiscreteFunctionType &discreteFunction )
      {
        interpolateFunction( function, discreteFunction );
      }

    private:
      static void interpolateFunction ( const Function &function, DiscreteFunctionType &discreteFunction, integral_constant< bool, true > )
      {
        interpolateDiscreteFunction( function, discreteFunction );
      }

      static void interpolateFunction ( const Function &function, DiscreteFunctionType &discreteFunction, integral_constant< bool, false > )
      {
        typedef GridFunctionAdapter< Function, GridPartType > GridFunctionType;

        const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();
        GridFunctionType dfAdapter( "function", function, dfSpace.gridPart() );
        interpolateDiscreteFunction( dfAdapter, discreteFunction );
      }

      template< class GridFunction >
      static void interpolateDiscreteFunction ( const GridFunction &function, DiscreteFunctionType &discreteFunction );
    };



    template< class Function, class DiscreteFunction >
    template< class GridFunction >
    inline void LagrangeInterpolation< Function, DiscreteFunction >
      ::interpolateDiscreteFunction ( const GridFunction &function,
                                      DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionType::DofType DofType;
      typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
      static const int dimRange = DiscreteFunctionSpaceType::dimRange;

      typedef typename GridFunction::LocalFunctionType FunctionLocalFunctionType;

      // set all DoFs to infinity
      const DofIteratorType dend = discreteFunction.dend();
      for( DofIteratorType dit = discreteFunction.dbegin(); dit != dend; ++dit )
        *dit = std::numeric_limits< DofType >::infinity();

      const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();
      for( const typename DiscreteFunctionSpaceType::EntityType &entity : dfSpace )
      {
        const LagrangePointSetType &lagrangePointSet
          = dfSpace.lagrangePointSet( entity );

        FunctionLocalFunctionType f_local = function.localFunction( entity );
        LocalFunctionType df_local = discreteFunction.localFunction( entity );

        // assume point based local dofs
        const int nop = lagrangePointSet.nop();
        int k = 0;
        for( int qp = 0; qp < nop; ++qp )
        {
          // if the first DoF for this point is already valid, continue
          if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
          {
            // evaluate the function in the Lagrange point
            RangeType phi;
            f_local.evaluate( lagrangePointSet[ qp ], phi );

            // assign the appropriate values to the DoFs
            for( int i = 0; i < dimRange; ++i, ++k )
              df_local[ k ] = phi[ i ];
          }
          else
            k += dimRange;
        }
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
