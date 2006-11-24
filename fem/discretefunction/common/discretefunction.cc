#ifndef DUNE_DISCRETEFUNCTION_CC
#define DUNE_DISCRETEFUNCTION_CC

#include <fstream>
#include "../../io/file/asciiparser.hh"

namespace Dune 
{
 
//************************************************************
//  Default Implementations 
//************************************************************
template <class DiscreteFunctionTraits>
void DiscreteFunctionDefault<DiscreteFunctionTraits>::
set(const RangeFieldType & value) 
{
  DofIteratorType endit = this->dend();
  for (DofIteratorType it = this->dbegin(); it != endit; ++it) 
  {
    (*it) = value;
  }
}

template <class DiscreteFunctionTraits>
void DiscreteFunctionDefault<DiscreteFunctionTraits>::clear() 
{
  this->set(RangeFieldType(0.0));
}

template <class DiscreteFunctionTraits>
void DiscreteFunctionDefault<DiscreteFunctionTraits>::
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
  DiscreteFunctionDefault<DiscreteFunctionTraits >&
  DiscreteFunctionDefault<DiscreteFunctionTraits >::
  assign(const MappingType& g) 
  {
    typedef DiscreteFunctionDefault<
      DiscreteFunctionTraits 
      > DiscreteFunctionDefaultType;
    
    const DiscreteFunctionDefaultType &gc = 
      static_cast<const DiscreteFunctionDefaultType &> ( g );
    
    assert(this->size() == gc.size());

    DofIteratorType endit = this->dend ();
    ConstDofIteratorType git = gc.dbegin ();
    for (DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) {
      *it = *git;
    }
    return *this;
  }


// operator +=
/** \todo This operator can add a discretefunction defined on all levels to another
 * one defined only on one level.  We should somehow issue a warning in this case.
 */
template<class DiscreteFunctionTraits>
DiscreteFunctionDefault<DiscreteFunctionTraits >&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator += (const MappingType& g) 
{
  typedef DiscreteFunctionDefault<
    DiscreteFunctionTraits 
    > DiscreteFunctionDefaultType;

  const DiscreteFunctionDefaultType &gc = 
    static_cast<const DiscreteFunctionDefaultType &> ( g );

  assert(this->size() == gc.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = gc.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it += *git;
  }
  return *this;
}


// operator -=
template<class DiscreteFunctionTraits>
DiscreteFunctionDefault<DiscreteFunctionTraits> &
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator -= ( const MappingType& g ) 
{
  typedef DiscreteFunctionDefault<
    DiscreteFunctionTraits 
    > DiscreteFunctionDefaultType;

  // cast to class discrete functions     
  const DiscreteFunctionDefaultType &gc = 
    static_cast<const DiscreteFunctionDefaultType &> ( g );

  assert(this->size() == gc.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = gc.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it -= *git;
  }
  return asImp();
}

// operator *=
template<class DiscreteFunctionTraits >
inline DiscreteFunctionDefault<DiscreteFunctionTraits>&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator*=(const typename DiscreteFunctionDefault<DiscreteFunctionTraits>::RangeFieldType & scalar)
{
  DofIteratorType endit = this->dend ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it) 
    *it *= scalar;

  return *this;
}

// operator /=
template<class DiscreteFunctionTraits>
inline DiscreteFunctionDefault<DiscreteFunctionTraits > &
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator/=(const typename DiscreteFunctionDefault<DiscreteFunctionTraits>::RangeFieldType & scalar)
{
  (*this) *= (1./scalar);
  return *this;
}


// add
template<class DiscreteFunctionTraits >
typename DiscreteFunctionTraits::DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
add(const DiscreteFunctionType& g, RangeFieldType scalar) 
{
  typedef DiscreteFunctionDefault<
    DiscreteFunctionTraits 
    > DiscreteFunctionDefaultType;
    
  const DiscreteFunctionDefaultType &gc = 
     static_cast<const DiscreteFunctionDefaultType &> ( g );

  assert(this->size() == gc.size());
    
  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = gc.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it += (*git) * scalar;
  }
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
