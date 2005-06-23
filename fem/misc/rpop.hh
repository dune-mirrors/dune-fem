#ifndef __DUNE_RPOP_HH__
#define __DUNE_RPOP_HH__

#include <dune/common/array.hh>
#include <dune/fem/common/objpointer.hh>

#include <dune/quadrature/barycenter.hh>

namespace Dune{

//***********************************************************************

    /** \brief I suppose it does the restriction/prolongation for FV
     * 
     */
template <class DofManagerType, class DiscreteFunctionType>
class NewRPOperatorFV
{
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
public:  
    //! Constructor
  NewRPOperatorFV (DofManagerType& dm,  
                   DiscreteFunctionType & df , 
                   RangeFieldType weight) : 
    dm_ (dm),
    df_ (df) , 
    vati_ ( df_.newLocalFunction() ) , 
    sohn_ ( df_.newLocalFunction() ) , 
    weight_(weight),
    chunkSize_(0) {}

  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, 
                       EntityType &son,
                       int chunkSize,
                       bool initialize ) const
  {
    // if weight < 0.0 , weight has not been calced
    assert(weight_ > 0.0);
   
    dm_.resizeChunk((father.globalIndex()+1)*
                    DiscreteFunctionType::FunctionSpaceType::DimRange,
                    chunkSize*
                    DiscreteFunctionType::FunctionSpaceType::DimRange);
 
    df_.localFunction( father, vati_ );
    df_.localFunction( son   , sohn_ );

    if(initialize)
    {
      for(int i=0; i<vati_.numberOfDofs(); i++)
      {
        vati_[i] = weight_ * sohn_[i];
      }
    }
    else 
    {
      for(int i=0; i<vati_.numberOfDofs(); i++)
      {
        vati_[i] += weight_ * sohn_[i];
      }
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, 
                      EntityType &son, 
                      int chunkSize,
                      bool initialize ) const
  {
    dm_.resizeChunk((son.globalIndex()+1)* 
                    DiscreteFunctionType::FunctionSpaceType::DimRange,
                    chunkSize*
                    DiscreteFunctionType::FunctionSpaceType::DimRange);

    df_.localFunction( father, vati_ );
    df_.localFunction( son   , sohn_ );

    for(int i=0; i<vati_.numberOfDofs(); i++)
    {
      sohn_[i] = vati_[i];
    }
  }

private:
  mutable DofManagerType& dm_;
  mutable DiscreteFunctionType & df_;

  mutable LocalFunctionType vati_;
  mutable LocalFunctionType sohn_;

  const RangeFieldType weight_;
  int chunkSize_;
};


/** @} end documentation group */

}

#endif
