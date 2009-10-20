#ifndef DUNE_L2PROJECTION_HH
#define DUNE_L2PROJECTION_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include "dgl2projection.hh"

namespace Dune 
{

// deprecated class (use DGL2Projection instead)
struct L2ProjectionImpl : public DGL2ProjectionImpl
{
  //! project function onto discrete function space  
  template <class FunctionImp, class DiscreteFunctionImp>
  DUNE_DEPRECATED static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd = -1)
  {
    DGL2ProjectionImpl :: project( f, discFunc, polOrd );
  }
};

/*======================================================================*/
/*! @ingroup L2ProjectionOperator
 *  \class L2Projection
 *  \brief The L2Projection class provides methods for projecting a function
 *         unto a given discrete function space. Note that this implementation
 *         assumes orthorgonal base functions!
 */
/*======================================================================*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class L2Projection : public Operator<DFieldType, RFieldType,DType , RType> {
 public:
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  //! Constructor taking degree for quadrature rule
  //! if no argument given a default value is chosen depending on the order
  //! in the discrete function space
  L2Projection(int polOrd = -1) : polOrd_(polOrd) {}

  //! apply L2 projection
  virtual void operator() (const DomainType& f, RangeType& discFunc) const 
  {
    if( discFunc.space().continuous() )
      DUNE_THROW(NotImplemented,"L2-Projection not implemented for contiuous spaces!"); 
    else 
      DGL2ProjectionImpl::project(f,discFunc,polOrd_);
  }

protected:
  const int polOrd_;
};

}
#endif
