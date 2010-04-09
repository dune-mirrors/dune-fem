#ifndef DUNE_BLOCKVECTORFUNCTION_INLINE_HH
#define DUNE_BLOCKVECTORFUNCTION_INLINE_HH

#include "blockvectorfunction.hh"

namespace Dune 
{

  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType &f ) 
  : BaseType( "no name ", f, lfFactory_ ), 
    lfFactory_( *this ),
    mapper_( f.blockMapper() ) ,
    memObject_( 0 ), 
    dofVec_( allocateDofStorage() ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}


  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const std::string &name,
                                     const DiscreteFunctionSpaceType &dfSpace )
  : BaseType( name, dfSpace, lfFactory_ ),
    lfFactory_( *this ),
    mapper_( dfSpace.blockMapper() ),
    memObject_( 0 ), 
    dofVec_( allocateDofStorage() ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &f,
                                     const DofStorageType &data )
  : BaseType( name, f , lfFactory_ ),
    lfFactory_( *this ),
    mapper_( f.blockMapper() ),
    memObject_( 0 ),
    dofVec_( const_cast<DofStorageType&> (data) ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const ThisType &other ) 
  : BaseType( other.name(), other.space(), lfFactory_ ),
    lfFactory_( *this ),
    mapper_( other.space().blockMapper() ),
    memObject_( 0 ),
    dofVec_( allocateDofStorage() ),
    leakPtr_( dofVec_ ),
    localFunc_( *this )
  {
    // copy values 
    dofVec_ = other.dofVec_;
  } 

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: ~BlockVectorDiscreteFunction ()
  {
    if( memObject_ ) delete memObject_; 
  }

template<class DiscreteFunctionSpaceType>
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::DofStorageType& 
BlockVectorDiscreteFunction< DiscreteFunctionSpaceType > :: allocateDofStorage()
{
  if( memObject_ != 0 ) 
    DUNE_THROW(InvalidStateException,"DofStorage already allocated!");
  
  std::pair< DofStorageInterface*, DofStorageType* > memPair
    = allocateManagedDofStorage( this->space().grid(),
                                 mapper_ ,
                                 this->name(),
                                 (DofStorageType *) 0 );
  // store memory 
  memObject_ = memPair.first;

  return *(memPair.second);
}

template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::clear ()
{
  const int size = dofVec_.size();
  for(int i=0; i<size; ++i) dofVec_[i] = 0.0; 
}

template< class DiscreteFunctionSpaceType >
inline void BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: print ( std::ostream &out ) const
{
  out << "BlockVectorDiscreteFunction '" << name() << "'" << std :: endl;
  
  const ConstDofIteratorType end = dend();
  for( ConstDofIteratorType dit = dbegin(); dit != end; ++dit )
    out << (*dit) << std :: endl;
}

//*************************************************************************
//  Interface Methods 
//*************************************************************************

#if 0
template<class DiscreteFunctionSpaceType> template <class EntityType>
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: LocalFunctionType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: localFunction ( const EntityType &en ) const
{
  return LocalFunctionType (en,*this);
}
#endif

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::DofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dbegin ()
{
  return DofIteratorType( dofVec_ , 0 );     
}


template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::DofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dend ( )
{
  return DofIteratorType( dofVec_  , dofVec_.size() );     
}

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::ConstDofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dbegin ( ) const
{
  DofIteratorType tmp ( dofVec_ , 0 );     
  return ConstDofIteratorType(tmp);
}

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::ConstDofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dend ( ) const 
{
  DofIteratorType tmp ( dofVec_ , dofVec_.size() );     
  return ConstDofIteratorType(tmp);
}

template< class DiscreteFunctionSpaceType >
inline void BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: addScaled ( const DiscreteFunctionType &g,
                 const RangeFieldType &s )
{
  int length = dofVec_.size();
  const DofStorageType &gvec = g.dofVec_;
  assert(length == gvec.size());
  dofVec_.axpy( s, gvec );
}

template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
enableDofCompression()
{
  if( memObject_ ) 
    memObject_->enableDofCompression();
}

} // end namespace Dune 
#endif
