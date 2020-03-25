#ifndef DUNE_FEM_DGL2PROJECTION_HH
#define DUNE_FEM_DGL2PROJECTION_HH
#warning "Deprecated header, use #include <dune/fem/space/common/interpolate.hh> instead!"

#include <type_traits>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include <dune/fem/space/common/interpolate.hh>

namespace Dune
{

  namespace Fem
  {

    // DGL2ProjectionImpl
    // ------------------

    // implementation of L2 projection for discontinuous spaces
    class DGL2ProjectionImpl
    {
      template <int dummy, bool hasLocalFunction>
      struct ProjectChooser
      {
        template <class FunctionImp, class FunctionSpace>
        class FunctionAdapter
        {
          const FunctionImp& function_;
        public:
          typedef FunctionSpace FunctionSpaceType;
          typedef typename FunctionSpaceType :: RangeType RangeType;
          typedef typename FunctionSpaceType :: DomainType DomainType;
          FunctionAdapter(const FunctionImp& f) : function_(f) {}

          void evaluate(const DomainType& local,
                        RangeType& ret) const
          {
            function_.evaluate( local , ret );
          }
        };

        template <class FunctionImp, class DiscreteFunctionImp>
        static void project(const FunctionImp& f,
                            DiscreteFunctionImp& discFunc,
                            const int polOrd)
        {
          // some typedefs
          typedef typename DiscreteFunctionImp :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
          typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
          typedef FunctionAdapter<FunctionImp, typename  DiscreteFunctionSpaceType :: FunctionSpaceType> FunctionAdapterType;
          // create function adapter in case of incorrect implementation
          FunctionAdapterType af( f );
          // create discrete function adapter
          GridFunctionAdapter< FunctionAdapterType, GridPartType> adapter(
              "DGL2projection::adapter" , f , discFunc.space().gridPart(), discFunc.space().order() );
          DGL2ProjectionImpl::projectFunction(adapter, discFunc, polOrd);
        }
      };

      template <int dummy>
      struct ProjectChooser<dummy,true>
      {
        template <class FunctionImp, class DiscreteFunctionImp>
        static void project(const FunctionImp& f,
                            DiscreteFunctionImp& discFunc,
                            const int polOrd )
        {
          DGL2ProjectionImpl::projectFunction(f, discFunc, polOrd);
        }
      };

    public:
      //! project function to discrete space
      template <class FunctionImp, class DiscreteFunctionImp>
      static void apply(const FunctionImp& f, DiscreteFunctionImp& discFunc )
      {
        project( f, discFunc );
      }

    protected:
      template <class FunctionImp, class DiscreteFunctionImp>
      static void projectFunction(const FunctionImp& func,
                                  DiscreteFunctionImp& discFunc,
                                  int polOrd = -1)
      {
        interpolate( func, discFunc );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DGL2PROJECTION_HH
