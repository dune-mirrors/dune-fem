#ifndef DUNE_PRODUCTFUNCTION_CC
#define DUNE_PRODUCTFUNCTION_CC



#include <algorithm>
// BUG: return values are initialized with *= 0. Problematic for nan
// values!
namespace Dune {

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
  inline ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
  ProductDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) :
    // DiscreteFunctionDefaultType ( f ), 
    name_ ("no name"),
    functionSpace_(f),
    functionSpace2_ (f2),
    dofVec_ ( f.size() * f2.size() )
{
   // std::cout << "Size old = " << dofVec_.size() << "  Size  f2 = " << functionSpace2_.size();
   // std::cout << " New size will be " << dofVec_.size() * functionSpace2_.size() << "\n";
}
// Desctructor 
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
~ProductDiscreteFunction()
{
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::clear ()
{
  const int size = dofVec_.size();
  DofStorageType &vec = dofVec_;
  for(int i=0; i<size; ++i) vec[i] = 0.0; 
  //set ( 0.0 ); 
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
print(std::ostream &s ) const
{
  RangeFieldType sum = 0.;
  ConstDofIteratorType enddof = dend ( );
  for(ConstDofIteratorType itdof = dbegin ( ); itdof != enddof; ++itdof) 
  {
    s << (*itdof) << " DofValue \n";
    sum += std::abs(*itdof);
  } 
  s << "sum = " << sum << "\n";
}
//*************************************************************************
//  Interface Methods 
//*************************************************************************
//return local function for global dof index dofIndex2 of space two
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*dofIndex2]));
 
}

//return local function for local dof index dofIndex2 of space two and given entity2
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
template < class Entity2Type>
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2,int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*space2().mapToGlobal( en2 ,dofIndex2)]));
 
}

//local function for given entity2, quadrature type and qudrature point number of second space and discrete function of first space
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class QuadratureType>
inline void ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2, const QuadratureType &quad2,  int pointNr, DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(i,quad2[pointNr],tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
	
}

//local function for given entity2, quadrature type and qudrature point number of second space and discrete function of first space
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class PointType>
inline void ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
localFunction(const Entity2Type &en2, const PointType &pt, DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(pt,tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( )
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend ()
{
  return dofVec_.end();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( ) const
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend () const 
{
  return dofVec_.end();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
addScaled( const ThisType &org, const RangeFieldType &c )
{
  dofVec_.axpy(org.dofVec_ , c);
}


} // end namespace 
#endif
