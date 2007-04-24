#ifndef DUNE_DGADAPTMANAGERIMP_HH
#define DUNE_DGADAPTMANAGERIMP_HH

//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachequad.hh>

//- local includes 
#include "dgspace.hh"

/************************************************************
1) Gewichte zwischen Vater/Sohn default Implementieren auf Gitter
   (father.weight(son) oder so
2) Caching der Basisfunktionen fuer Vater-BasisFunktionen fuer
   Kinderquadraturen
*************************************************************/

namespace Dune{
/** @ingroup RestrictProlongInterface
**/

//***********************************************************************
/** \brief This is a restriction/prolongation operator for DG data. 
 */
 template <template <class> class DiscFunc,
  	  class FunctionSpaceImp, 
   	  class GridPartImp, 
   	  int polOrd, 
   	  template <class> class StorageImp> 
 class RestrictProlongDefault<DiscFunc<
   DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp> 
 > > : 
  public RestrictProlongInterface<RestrictProlongTraits<RestrictProlongDefault<DiscFunc<
    DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp
     > > > > >
 {
 public:
   typedef DiscFunc<DiscontinuousGalerkinSpace<FunctionSpaceImp, 
					       GridPartImp, 
					      polOrd,
					      StorageImp> > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType :: GridPartType GridPartType;
  typedef typename FunctionSpaceType :: GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef CachingQuadrature<GridPartType,0> QuadratureType;
  typedef typename GridType::template Codim<0>::Entity::Geometry Geometry;
 public:  
  //! Constructor
  RestrictProlongDefault ( DiscreteFunctionType & df ) : 
    df_ (df) , quadord_(2*df.space().order()),
    weight_(-1.0)
  {
    // make sure that index set is used that can handle adaptivity 
    assert( (Capabilities::IsUnstructured<GridType>::v) ? (df.space().indexSet().adaptive()) : true );
  }
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const
  {
    // volume of son / volume of father  
    weight_ = val; 
  }

  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const
  {
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    assert( !father.isLeaf() );
    const RangeFieldType weight = 
      (weight_ < 0.0) ? (this->calcWeight(father,son)) : weight_; 

    LocalFunctionType vati_ = df_.localFunction( father);
    LocalFunctionType sohn_ = df_.localFunction( son   );

    QuadratureType quad(son,quadord_);
    const typename FunctionSpaceType::BaseFunctionSetType & baseset =
      vati_.baseFunctionSet();
    const int nop=quad.nop();
    const Geometry& geometryInFather = son.geometryInFather();

    const int vati_numDofs = vati_.numDofs();
    if(initialize) 
    {
      for(int i=0; i<vati_numDofs; ++i) 
      {
      	vati_[i] = 0.0;
      }
    }
    
    for(int qP = 0; qP < nop; ++qP) 
    {
      sohn_.evaluate(quad,qP,ret);
      for(int i=0; i<vati_numDofs; ++i) 
      {
      	baseset.evaluate(i,geometryInFather.global(quad.point(qP)),phi);
      	vati_[i] += quad.weight(qP) * weight * (ret * phi) ;
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
    const int sohn_numDofs = sohn_.numDofs();
    for(int i=0; i<sohn_numDofs; ++i) sohn_[i] = 0.;

    QuadratureType quad(son,quadord_);
    const typename FunctionSpaceType::BaseFunctionSetType & baseset =
      sohn_.baseFunctionSet();
    const Geometry& geometryInFather = son.geometryInFather();
    const int nop=quad.nop();
    for(int qP = 0; qP < nop; ++qP) 
    {
      vati_.evaluate(geometryInFather.global(quad.point(qP)),ret);
      
      for(int i=0; i<sohn_numDofs; ++i) {
      	baseset.evaluate(i,quad,qP,phi);
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
 class RestrictProlongDefault<DiscFunc<
   DiscontinuousGalerkinSpace<FunctionSpaceImp,GridPartImp,0,StorageImp> 
 > > : 
  public RestrictProlongInterface<RestrictProlongTraits<RestrictProlongDefault<DiscFunc<
    DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, 0,StorageImp
     > > > > >
 
 {
   public:
   typedef DiscFunc<DiscontinuousGalerkinSpace<FunctionSpaceImp, 
					       GridPartImp, 
					       0,
					       StorageImp> > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType ::GridPartType GridPartType;
  typedef typename FunctionSpaceType ::GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef CachingQuadrature<GridPartType,0> QuadratureType;
  public:  
  //! Constructor
  RestrictProlongDefault ( DiscreteFunctionType & df ) : 
    df_ (df),
    weight_(-1.0)
  {
    // make sure that index set is used that can handle adaptivity 
    assert( (Capabilities::IsUnstructured<GridType>::v) ? (df.space().indexSet().adaptive()) : true );
  }
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const
  {
    // volume of son / volume of father  
    weight_ = val; 
  }
  
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const
  {
    assert( !father.isLeaf() );

    // if weight < 0.0 , weight has not been calculated
    const RangeFieldType weight = 
      (weight_ < 0.0) ? (this->calcWeight(father,son)) : weight_; 
    
    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );

    const int numDofs = vati.numDofs();
    assert( sohn.numDofs() == numDofs );
    if(initialize) {
      for(int i=0; i<numDofs; i++) {
        vati[i] = weight * sohn[i];
      }
    }
    else 
    {
      for(int i=0; i<numDofs; i++) {
        vati[i] += weight * sohn[i];
      }
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );
    const int numDofs = vati.numDofs();
    assert( sohn.numDofs() == numDofs );
    for(int i=0; i<numDofs; i++) {
      sohn[i] = vati[i];
    }
  }

private:
  mutable DiscreteFunctionType & df_;
  mutable RangeFieldType weight_;
};


}

#endif
