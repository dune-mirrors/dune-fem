#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_INLINE_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_INLINE_HH

#include "blockvectorfunction.hh"

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunctionSpaceImp >
    inline ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    ::ISTLBlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( "no name ", dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType * ) ) * space().blockMapper().maxNumDofs() * localBlockSize ),
        blockMapper_( space().blockMapper() ),
        memObject_( 0 ),
        dofVec_( allocateDofStorage() ),
        leakPtr_( dofVec_ )
    {}


    template< class DiscreteFunctionSpaceImp >
    inline ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    ::ISTLBlockVectorDiscreteFunction ( const std::string &name,
                                        const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType * ) ) * space().blockMapper().maxNumDofs() * localBlockSize ),
        blockMapper_( space().blockMapper() ),
        memObject_( 0 ),
        dofVec_( allocateDofStorage() ),
        leakPtr_( dofVec_ )
    {}


    template< class DiscreteFunctionSpaceImp >
    inline ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    ::ISTLBlockVectorDiscreteFunction ( const std::string &name,
                                        const DiscreteFunctionSpaceType &dfSpace,
                                        const DofStorageType &data )
      : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( std::max( sizeof( DofType ), sizeof( DofType * ) ) * space().blockMapper().maxNumDofs() * localBlockSize ),
        blockMapper_( space().blockMapper() ),
        memObject_( 0 ),
        dofVec_( const_cast< DofStorageType & >(data) ),
        leakPtr_( dofVec_ )
    {}


    template< class DiscreteFunctionSpaceImp >
    inline ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    ::ISTLBlockVectorDiscreteFunction ( const ThisType &other )
      : BaseType( other.name(), other.space(), LocalDofVectorAllocatorType( &ldvStack_ ) ),
        ldvStack_( other.ldvStack_ ),
        blockMapper_( space().blockMapper() ),
        memObject_( 0 ),
        dofVec_( allocateDofStorage() ),
        leakPtr_( dofVec_ )
    {
      // copy values
      dofVec_ = other.dofVec_;
    }


    template< class DiscreteFunctionSpaceType >
    inline typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::DofStorageType&
    ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::allocateDofStorage()
    {
      if( memObject_ != 0 )
        DUNE_THROW( InvalidStateException, "DofStorage already allocated!" );

      std::pair< Fem::DofStorageInterface *, DofStorageType * > memPair
        = Fem::allocateManagedDofStorage( this->space().gridPart().grid(),
                                          blockMapper_,
                                          this->name(),
                                          (DofStorageType *) 0 );
      // store memory
      memObject_ = memPair.first;

      return *(memPair.second);
    }

    template< class DiscreteFunctionSpaceType >
    inline void ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::clear ()
    {
      const int size = dofVec_.size();
      for( int i = 0; i < size; ++i )
        dofVec_[ i ] = 0.0;
    }

    template< class DiscreteFunctionSpaceType >
    inline void ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
    ::print ( std::ostream &out ) const
    {
      const ConstDofIteratorType end = dend();
      for( ConstDofIteratorType dit = dbegin(); dit != end; ++dit )
        out << (*dit) << std::endl;
    }

    //*************************************************************************
    //  Interface Methods
    //*************************************************************************

    template< class DiscreteFunctionSpaceType >
    inline typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::DofIteratorType
    ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::dbegin ()
    {
      return DofIteratorType( dofVec_, 0 );
    }


    template< class DiscreteFunctionSpaceType >
    inline typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::DofIteratorType
    ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::dend ()
    {
      return DofIteratorType( dofVec_, dofVec_.size() );
    }

    template< class DiscreteFunctionSpaceType >
    inline typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::ConstDofIteratorType
    ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::dbegin () const
    {
      DofIteratorType tmp( dofVec_, 0 );
      return ConstDofIteratorType( tmp );
    }

    template< class DiscreteFunctionSpaceType >
    inline typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::ConstDofIteratorType
    ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::dend () const
    {
      DofIteratorType tmp( dofVec_, dofVec_.size() );
      return ConstDofIteratorType( tmp );
    }

    template< class DiscreteFunctionSpaceType >
    inline void ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
    ::axpy ( const RangeFieldType &s, const DiscreteFunctionType &g )
    {
      const DofStorageType &gvec = g.dofVec_;
      assert( dofVec_.size() == gvec.size() );
      dofVec_.axpy( s, gvec );
    }

    template< class DiscreteFunctionSpaceType >
    inline void ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >::
    enableDofCompression ()
    {
      if( memObject_ )
        memObject_->enableDofCompression();
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORFUNCTION_INLINE_HH
