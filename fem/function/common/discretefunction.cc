#ifndef DUNE_DISCRETEFUNCTION_CC
#define DUNE_DISCRETEFUNCTION_CC

#include <fstream>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune 
{
 
//************************************************************
//  Default Implementations 
//************************************************************
template <class DiscreteFunctionTraits>
inline void DiscreteFunctionDefault<DiscreteFunctionTraits>::clear() 
{
  DofIteratorType endit = this->dend();
  for (DofIteratorType it = this->dbegin(); it != endit; ++it) 
  {
    (*it) = 0;
  }
}

template <class DiscreteFunctionTraits>
inline void DiscreteFunctionDefault<DiscreteFunctionTraits>::
addScaled(const DiscreteFunctionType& g, const RangeFieldType& c) {
  assert(this->size() == g.size());
  DofIteratorType endit = this->dend();
  ConstDofIteratorType oit = g.dbegin();
  for (DofIteratorType it = this->dbegin(); it != endit; ++it, ++oit) 
  {
    (*it) += (*oit) * c;
  }
}

// scalarProductDofs
template <class DiscreteFunctionTraits>
inline typename DiscreteFunctionTraits::DiscreteFunctionSpaceType::RangeFieldType 
DiscreteFunctionDefault<DiscreteFunctionTraits>::
scalarProductDofs(const DiscreteFunctionType& g) const
{
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType; 
  assert(this->size() == g.size());

  RangeFieldType skp = 0.;

  ConstDofIteratorType endit = this->dend ();
  ConstDofIteratorType git =  g.dbegin ();

  // multiply
  for(ConstDofIteratorType it = this->dbegin(); it != endit; ++it,++git)
  {   
    skp += (*it) * (*git);
  }
  
  return skp;
}

// operator=
template<class DiscreteFunctionTraits>
inline void 
DiscreteFunctionDefault<DiscreteFunctionTraits >::
assign(const DiscreteFunctionType& g) 
{
  assert(this->size() == g.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = g.dbegin ();
  for (DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) {
    *it = *git;
  }
}

// operator +=
/** \todo This operator can add a discretefunction defined on all levels to another
 * one defined only on one level.  We should somehow issue a warning in this case.
 */
template<class DiscreteFunctionTraits>
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator += ( const DiscreteFunctionType& g ) 
{
  assert(this->size() == g.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = g.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it += *git;
  }
  return asImp();
}

// operator -=
template<class DiscreteFunctionTraits>
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator -= ( const DiscreteFunctionType& g ) 
{
  assert(this->size() == g.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = g.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it -= *git;
  }
  return asImp();
}

// operator *=
template<class DiscreteFunctionTraits >
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator*=(const typename DiscreteFunctionDefault<DiscreteFunctionTraits>::RangeFieldType & scalar)
{
  DofIteratorType endit = this->dend ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it) 
  {
    *it *= scalar;
  }
  return asImp();
}

// operator /=
template<class DiscreteFunctionTraits>
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator/=(const typename DiscreteFunctionDefault<DiscreteFunctionTraits>::RangeFieldType & scalar)
{
  (*this) *= (RangeFieldType(1)/scalar);
  return asImp();
}


// print 
template<class DiscreteFunctionTraits >
inline void 
DiscreteFunctionDefault<DiscreteFunctionTraits >::
print(std::ostream & s) const 
{
  typedef DiscreteFunctionDefault<
    DiscreteFunctionTraits 
    > DiscreteFunctionDefaultType;
    
  s << this->name() << std::endl;
  ConstDofIteratorType endit = this->dend ();
  for(ConstDofIteratorType it = this->dbegin(); it != endit; ++it) 
  {
    s << (*it) << std::endl;
  }
}

} // end namespace Dune
#endif
