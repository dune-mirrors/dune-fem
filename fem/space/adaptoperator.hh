#ifndef DUNE_DGADAPTOPERATORIMP_HH
#define DUNE_DGADAPTOPERATORIMP_HH

#include <dune/fem/transfer/adaptoperator.hh>

#include "../quadrature/cachequad.hh"
#include "dgspace.hh"

namespace Dune{

//***********************************************************************

/** \brief This is a general restriction/prolongation operator 
 */
 template <class DiscreteFunctionType> class RestProlOperator;
/** \brief This is a restriction/prolongation operator for DG data. 
 */
 template <template <class> class DiscFunc,
  	  class FunctionSpaceImp, 
   	  class GridPartImp, 
   	  int polOrd, 
   	  template <class> class StorageImp> 
 class RestProlOperator<DiscFunc<
   DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp> 
 > > 
 {
   typedef DiscFunc<DiscontinuousGalerkinSpace<FunctionSpaceImp, 
					       GridPartImp, 
					      polOrd,
					      StorageImp> > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionType::GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef CachingQuadrature<GridType,0> QuadratureType;
public:  
  //! Constructor
  RestProlOperator ( DiscreteFunctionType & df , GeometryType eltype ) : 
    df_ (df) , quadord_(2*df.getFunctionSpace().polynomOrder()),
    weight_(-1.0)
  {
  }

  //! calculates the weight, i.e. (volume son)/(volume father)
  template <class EntityType>
  void calcFatherChildWeight (EntityType &father, EntityType &son) const
  {
  }
  
  //! the weight can also be seted
  void setFatherChildWeight (const RangeFieldType& val) const
  {
  }
  
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const
  {
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    assert( !father.isLeaf() );

    LocalFunctionType vati_ =df_.localFunction( father);
    LocalFunctionType sohn_ =df_.localFunction( son   );

    QuadratureType quad(son,quadord_);
    const typename FunctionSpaceType::BaseFunctionSetType & baseset =
      vati_.getBaseFunctionSet();
    const int nop=quad.nop();
    if(initialize) {
      for(int qP = 0; qP < nop; qP++) {
	sohn_.evaluate(son,quad,qP,ret);
	for(int i=0; i<sohn_.numDofs(); i++) {
	  baseset.eval(i,son.geometryInFather().global(quad.point(qP)),phi);
	  vati_[i] = quad.weight(qP) * (ret * phi) ;
	}
      }
    }
    else 
    {
      for(int i=0; i<vati_.numDofs(); i++)
      {
        vati_[i] += weight_ * sohn_[i];
      }
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    //assert( son.state() == REFINED );
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    LocalFunctionType vati_ = df_.localFunction( father);
    LocalFunctionType sohn_ = df_.localFunction( son   );
    for(int i=0; i<sohn_.numDofs(); i++) 
      sohn_[i] = 0.;

    QuadratureType quad(son,quadord_);
    const typename FunctionSpaceType::BaseFunctionSetType & baseset =
      sohn_.getBaseFunctionSet();
    const int nop=quad.nop();
    for(int qP = 0; qP < nop; qP++) {
      vati_.evaluateLocal(father,
			  son.geometryInFather().global(quad.point(qP)),ret);
      for(int i=0; i<sohn_.numDofs(); i++) {
	baseset.eval(i,quad,qP,phi);
	sohn_[i] += quad.weight(qP) * (ret * phi) ;
      }
    }
  }

private:
  mutable DiscreteFunctionType & df_;

  int quadord_;
  mutable RangeFieldType weight_;
};


/** \brief This is a restriction/prolongation operator for DG data of order zero. 
 */
 template <template <class> class DiscFunc,
  	  class FunctionSpaceImp, 
   	  class GridPartImp, 
   	  template <class> class StorageImp> 
 class RestProlOperator<DiscFunc<
   DiscontinuousGalerkinSpace<FunctionSpaceImp,GridPartImp,0,StorageImp> 
 > > 
 {
   typedef DiscFunc<DiscontinuousGalerkinSpace<FunctionSpaceImp, 
					       GridPartImp, 
					       0,
					       StorageImp> > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionType::GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef CachingQuadrature<GridType,0> QuadratureType;
public:  
  //! Constructor
  RestProlOperator ( DiscreteFunctionType & df , GeometryType eltype ) : 
    df_ (df),
    weight_(-1.0)
  {}

  //! calculates the weight, i.e. (volume son)/(volume father)
  template <class EntityType>
  void calcFatherChildWeight (EntityType &father, EntityType &son) const
  {
  }
  
  //! the weight can also be seted
  void setFatherChildWeight (const RangeFieldType& val) const
  {
  }
  
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const
  {
    assert( !father.isLeaf() );

    // if weight < 0.0 , weight has not been calculated
    assert(weight_ > 0.0);
    
    LocalFunctionType vati_ =df_.localFunction( father);
    LocalFunctionType sohn_ =df_.localFunction( son   );

    if(initialize)
    {
      for(int i=0; i<vati_.numDofs(); i++)
      {
        vati_[i] = weight_ * sohn_[i];
      }
    }
    else 
    {
      for(int i=0; i<vati_.numDofs(); i++)
      {
        vati_[i] += weight_ * sohn_[i];
      }
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    //assert( son.state() == REFINED );
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    LocalFunctionType vati_ = df_.localFunction( father);
    LocalFunctionType sohn_ = df_.localFunction( son   );
    for(int i=0; i<vati_.numDofs(); i++)
    {
      sohn_[i] = vati_[i];
    }
  }

private:
  mutable DiscreteFunctionType & df_;
  mutable RangeFieldType weight_;
};


/** @} end documentation group */

}

#endif
