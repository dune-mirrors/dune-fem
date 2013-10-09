#ifndef DUNE_FEM_L2PROJECTION_HH
#define DUNE_FEM_L2PROJECTION_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include "dgl2projection.hh"

namespace Dune 
{

  namespace Fem 
  { 

    /*! @ingroup L2ProjectionOperator
     *  \class L2Projection
     *  \brief The L2Projection class provides methods for projecting a function
     *         unto a given discrete function space. Note that this implementation
     *         assumes orthorgonal base functions!
     */
    template < class DType, class RType>
    class L2Projection : public Operator<DType , RType> 
    {
    public:
      //! domain function type 
      typedef DType  DomainType;
      //! range function type 
      typedef RType  RangeType;

      /** \brief Constructor 
       *    \param  quadOrder      degree for quadrature rule (default = 2*space.order())
       *    \param  doCommunicate  apply communication for the result (default = true)
       */    
      explicit L2Projection( const int quadOrder = -1, 
                             const bool doCommunicate = true ) 
        : quadOrder_( quadOrder ), doCommunicate_( doCommunicate ) 
      {}

      /** \brief  calculates the L2 projection of a function onto the 
       *          discrete space discreteFunction belongs to.
       *  \param  function          function to be projected 
       *  \param  discreteFunction  discrete result of projection  */
      virtual void operator() ( const DomainType& function, RangeType& discreteFunction ) const 
      {
        if( discreteFunction.space().continuous() )
          DUNE_THROW(NotImplemented,"L2-Projection not implemented for contiuous spaces!"); 
        else 
          DGL2ProjectionImpl::project( function, discreteFunction, quadOrder_, doCommunicate_ );
      }

    protected:
      const int quadOrder_;          // order of quadrature  
      const bool doCommunicate_ ; // true if communication is applied for the result 
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_L2PROJECTION_HH
