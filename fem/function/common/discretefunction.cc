#ifndef DUNE_DISCRETEFUNCTION_CC
#define DUNE_DISCRETEFUNCTION_CC

#include <fstream>

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/streams/streams.hh>

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

  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: addScaled ( const DiscreteFunctionType &g,
                   const RangeFieldType &s )
  {
    assert( this->size() == g.size() );
    const DofIteratorType end = this->dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = this->dbegin(); it != end; ++it, ++git )
      (*it) += s * (*git);
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

  template< class DiscreteFunctionTraits >
  inline bool DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: operator== ( const DiscreteFunctionType &g ) const
  {
    if( size() != g.size() )
      return false;
    
    const ConstDofIteratorType end = dend();

    ConstDofIteratorType fit = dbegin();
    ConstDofIteratorType git = g.dbegin();
    for( ; fit != end; ++fit, ++git )
      if( *fit != *git )
        return false;
    
    return true;
  }


  // print 
  template< class DiscreteFunctionTraits >
  inline void DiscreteFunctionDefault<DiscreteFunctionTraits >
    :: print ( std::ostream &out ) const
  {
    out << name() << std::endl;
    
    const ConstDofIteratorType end = dend();
    for( ConstDofIteratorType dit = dbegin(); dit != end; ++dit )
      out << (*dit) << std::endl;
  }



  template< class StreamTraits, class DiscreteFunctionTraits >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const DiscreteFunctionInterface< DiscreteFunctionTraits > &df )
  {
    typedef DiscreteFunctionInterface< DiscreteFunctionTraits > DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: ConstDofIteratorType IteratorType;
   
    out << df.size();

    const IteratorType end = df.dend();
    for( IteratorType it = df.dbegin(); it != end; ++it )
      out << *it;

    return out;
  }



  template< class StreamTraits, class DiscreteFunctionTraits >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 DiscreteFunctionInterface< DiscreteFunctionTraits > &df )
  {
    typedef DiscreteFunctionInterface< DiscreteFunctionTraits > DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DofIteratorType IteratorType;

    int size;
    in >> size;
    if( size != df.size() )
      DUNE_THROW( IOError, "Trying to read discrete function of different size." );

    const IteratorType end = df.dend();
    for( IteratorType it = df.dbegin(); it != end; ++it )
      in >> *it;

    return in;
  }

} // end namespace Dune
#endif
