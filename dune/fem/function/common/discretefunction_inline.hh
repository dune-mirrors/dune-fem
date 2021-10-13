#ifndef DUNE_FEM_DISCRETEFUNCTION_INLINE_HH
#define DUNE_FEM_DISCRETEFUNCTION_INLINE_HH

#include <fstream>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/common/bindguard.hh>
#include <dune/fem/common/memory.hh>
#include <dune/fem/gridpart/common/persistentindexset.hh>
#include <dune/fem/io/streams/streams.hh>

#include "discretefunction.hh"

namespace Dune
{

  namespace Fem
  {

    // DiscreteFunctionDefault
    // -----------------------

    template< class Impl >
    inline DiscreteFunctionDefault< Impl >
      :: DiscreteFunctionDefault ( const std::string &name,
                                   const DiscreteFunctionSpaceType &dfSpace )
    : dfSpace_( referenceToSharedPtr( dfSpace ) ),
      ldvStack_( std::max( std::max( sizeof( DofType ), sizeof( DofType* ) ),
                           sizeof(typename LocalDofVectorType::value_type) ) // for PetscDiscreteFunction
                 * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
      ldvAllocator_( &ldvStack_ ),
      localFunction_( space() ),
      name_( name ),
      scalarProduct_( dfSpace )
    {
    }


    template< class Impl >
    inline DiscreteFunctionDefault< Impl >::DiscreteFunctionDefault ( std::string name, std::shared_ptr< const DiscreteFunctionSpaceType > dfSpace )
      : dfSpace_( std::move( dfSpace ) ),
        ldvStack_( std::max( std::max( sizeof( DofType ), sizeof( DofType* ) ),
                             sizeof(typename LocalDofVectorType::value_type) ) // for PetscDiscreteFunction
                   * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
        ldvAllocator_( &ldvStack_ ),
        localFunction_( space() ),
        name_( std::move( name ) ),
        scalarProduct_( dfSpace )
    {}


    template< class Impl >
    inline DiscreteFunctionDefault< Impl >::DiscreteFunctionDefault ( const DiscreteFunctionDefault &other )
      : BaseType( static_cast< const BaseType & >( other ) ),
        dfSpace_( other.dfSpace_ ),
        ldvStack_( std::max( std::max( sizeof( DofType ), sizeof( DofType* ) ),
                             sizeof(typename LocalDofVectorType::value_type) ) // for PetscDiscreteFunction
                   * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
        ldvAllocator_( &ldvStack_ ),
        localFunction_( space() ),
        name_( other.name_ ),
        scalarProduct_( other.scalarProduct_ )
    {
      if( other.assembleOperation_ != std::type_index( typeid( void ) ) )
        DUNE_THROW( InvalidStateException, "Cannot copy discrete function during assembly" );
      assert( other.assembleCount_ == 0 );
    }


    template< class Impl >
    inline DiscreteFunctionDefault< Impl >
      :: DiscreteFunctionDefault ( DiscreteFunctionDefault && other )
    : BaseType( static_cast< BaseType&& >( other ) ),
      dfSpace_( std::move( other.dfSpace_ ) ),
      ldvStack_( std::move( other.ldvStack_ ) ),
      ldvAllocator_( &ldvStack_ ),
      localFunction_( space() ),
      name_( std::move( other.name_ ) ),
      scalarProduct_( std::move( other.scalarProduct_ ) )
    {
      if( other.assembleOperation_ != std::type_index( typeid( void ) ) )
        DUNE_THROW( InvalidStateException, "Cannot move discrete function during assembly" );
      assert( other.assembleCount_ == 0 );
    }


    template< class Impl >
    inline void DiscreteFunctionDefault<Impl >
      :: print ( std::ostream &out ) const
    {
      const auto end = BaseType :: dend();
      for( auto dit = BaseType :: dbegin(); dit != end; ++dit )
        out << (*dit) << std::endl;
    }


    template< class Impl >
    inline bool DiscreteFunctionDefault< Impl >
      :: dofsValid () const
    {
      const auto end = BaseType :: dend();
      for( auto it = BaseType :: dbegin(); it != end; ++it )
      {
        if( ! std::isfinite( *it ) )
          return false ;
      }

      return true;
    }


    template< class Impl >
    template< class DFType >
    inline void DiscreteFunctionDefault< Impl >
      ::axpy ( const RangeFieldType &s, const DiscreteFunctionInterface< DFType > &g )
    {
      if( BaseType::size() != g.size() )
        DUNE_THROW(InvalidStateException,"DiscreteFunctionDefault: sizes do not match in axpy");

      // apply axpy to all dofs from g
      const auto end = BaseType::dend();
      auto git = g.dbegin();
      for( auto it = BaseType::dbegin(); it != end; ++it, ++git )
        *it += s * (*git );
    }


    template< class Impl >
    template< class DFType >
    inline void DiscreteFunctionDefault< Impl >
      ::assign ( const DiscreteFunctionInterface< DFType > &g )
    {
      if( BaseType::size() != g.size() )
      {
        std::cout << BaseType::size() << " vs " << g.size() << std::endl;
        DUNE_THROW(InvalidStateException,"DiscreteFunctionDefault: sizes do not match in assign");
      }

      // copy all dofs from g to this
      const auto end = BaseType::dend();
      auto git = g.dbegin();
      for( auto it = BaseType::dbegin(); it != end; ++it, ++git )
        *it = *git;
    }


    template< class Impl >
    template< class Operation >
    inline typename DiscreteFunctionDefault< Impl >
      :: template CommDataHandle< Operation > :: Type
    DiscreteFunctionDefault< Impl > :: dataHandle ( const Operation &operation )
    {
      return BaseType :: space().createDataHandle( asImp(), operation );
    }


    template< class Impl >
    template< class Functor >
    inline void DiscreteFunctionDefault< Impl >
      ::evaluateGlobal ( const DomainType &x, Functor functor ) const
    {
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      EntitySearch< GridPartType, EntityType::codimension > entitySearch( BaseType::space().gridPart() );

      const auto& entity(entitySearch( x ));
      const auto geometry = entity.geometry();

      // bind local function to entity
      auto guard = bindGuard( localFunction_, entity );
      // obtain dofs (since it's a temporary local function)
      getLocalDofs( entity, localFunction_.localDofVector());
      // evaluate functor
      functor( geometry.local( x ), localFunction_ );
    }


    template< class Impl >
    template< class DFType >
    inline typename DiscreteFunctionDefault< Impl > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Impl >
      ::operator+= ( const DiscreteFunctionInterface< DFType > &g )
    {
      if( BaseType::size() != g.size() )
        DUNE_THROW(InvalidStateException,"DiscreteFunctionDefault: sizes do not match in operator +=");

      const auto end = BaseType::dend();
      auto git = g.dbegin();
      for( auto it = BaseType::dbegin(); it != end; ++it, ++git )
        *it += *git;
      return asImp();
    }


    template< class Impl >
    template< class DFType >
    inline typename DiscreteFunctionDefault< Impl > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Impl >
      ::operator-= ( const DiscreteFunctionInterface< DFType > &g )
    {
      if( BaseType::size() != g.size() )
        DUNE_THROW(InvalidStateException,"DiscreteFunctionDefault: sizes do not match in operator -=");

      const auto end = BaseType :: dend();
      auto git = g.dbegin();
      for( auto it = BaseType :: dbegin(); it != end; ++it, ++git )
        *it -= *git;
      return asImp();
    }


    template< class Impl >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Impl >
      :: read ( InStreamInterface< StreamTraits > &in )
    {
      auto versionId = in.readUnsignedInt();
      if( versionId < DUNE_VERSION_ID(0,9,1) )
        DUNE_THROW( IOError, "Trying to read outdated file." );
      else if( versionId > DUNE_MODULE_VERSION_ID(DUNE_FEM) )
        std :: cerr << "Warning: Reading discrete function from newer version: "
                    << versionId << std :: endl;

      // verify space id for files written with dune-fem version 1.5 or newer
      if( versionId >= DUNE_VERSION_ID(1,5,0) )
      {
        // make sure that space of discrete function matches the space
        // of the data that was written
        const auto spaceId = space().type();
        int mySpaceIdInt;
        in >> mySpaceIdInt;
        const auto mySpaceId = static_cast<DFSpaceIdentifier>(mySpaceIdInt);

        if( spaceId != mySpaceId )
          DUNE_THROW( IOError, "Trying to read discrete function from different space: DFSpace (" << spaceName( spaceId ) << ") != DataSpace (" << spaceName( mySpaceId ) << ")" );
      }

      // read name
      in >> name_;

      // read size as integer
      int32_t mysize;
      in >> mysize;

      // check size
      if( static_cast< int >( mysize ) != BaseType::size() )
        DUNE_THROW( IOError, "Trying to read discrete function of different size." );
      if( BaseType::size() != static_cast< int >( this->space().size() ) )
        DUNE_THROW( InvalidStateException, "Trying to read discrete function in uncompressed state." );

      // read all dofs
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::LocalBlockIndices LocalBlockIndices;
      typedef typename BlockMapperType::SizeType SizeType;

      const BlockMapperType &blockMapper = this->space().blockMapper();
      const SizeType nBlocks = blockMapper.size();
      for( SizeType i = 0; i < nBlocks; ++i )
      {
        auto &&block = this->dofVector()[ i ];
        Hybrid::forEach( LocalBlockIndices(), [ &in, &block ] ( auto j ) { in >> block[ j ]; } );
      }

      // a discrete function that is part of backup/restore
      // most likely needs to keep the dofs consistent
      asImp().enableDofCompression();
    }


    template< class Impl >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Impl >
      :: write ( OutStreamInterface< StreamTraits > &out ) const
    {
      unsigned int versionId = DUNE_MODULE_VERSION_ID(DUNE_FEM);
      out << versionId ;

      // write space id to for testing when function is read
      auto spaceId = space().type();
      out << spaceId ;

      // write name
      out << name_;

      // only allow write when vector is compressed
      if( BaseType::size() != static_cast< int >( this->space().size() ) )
        DUNE_THROW( InvalidStateException, "Trying to write DiscreteFunction in uncompressed state." );

      // write size as 32-bit integer
      const int32_t mysize = BaseType::size();
      out << mysize;

      // write all dofs
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::LocalBlockIndices LocalBlockIndices;
      typedef typename BlockMapperType::SizeType SizeType;

      const BlockMapperType &blockMapper = this->space().blockMapper();
      const SizeType nBlocks = blockMapper.size();
      for( SizeType i = 0; i < nBlocks; ++i )
      {
        const auto block = this->dofVector()[ i ];
        Hybrid::forEach( LocalBlockIndices(), [ &out, &block ] ( auto j ) { out << block[ j ]; } );
      }
    }


    template< class Impl >
    void DiscreteFunctionDefault< Impl >
      :: insertSubData()
    {
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      IndexSetType& indexSet = (IndexSetType&)space().indexSet();
      if( Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::v )
      {
        auto persistentIndexSet = Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::map( indexSet );

        // this marks the index set in the DofManager's list of index set as persistent
        if( persistentIndexSet )
          persistentIndexSet->addBackupRestore();
      }
    }

    template< class Impl >
    void DiscreteFunctionDefault< Impl >
      :: removeSubData()
    {
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      IndexSetType& indexSet = (IndexSetType&)space().indexSet();
      if( Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::v )
      {
        auto persistentIndexSet = Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::map( indexSet );

        // this unmarks the index set in the DofManager's list of index set as persistent
        if( persistentIndexSet )
          persistentIndexSet->removeBackupRestore();
      }
    }


    template< class Impl >
    template< class DFType >
    inline bool DiscreteFunctionDefault< Impl >
      :: compare ( const DiscreteFunctionInterface< DFType > &g ) const
    {
      if( BaseType :: size() != g.size() )
        return false;

      const auto end = BaseType :: dend();

      auto fit = BaseType :: dbegin();
      auto git = g.dbegin();
      for( ; fit != end; ++fit, ++git )
      {
        if( std::abs( *fit - *git ) > 1e-15 )
          return false;
      }

      return true;
    }


    // Stream Operators
    // ----------------

    /** \brief write a discrete function into an STL stream
     *  \relates DiscreteFunctionInterface
     *
     *  \param[in]  out  STL stream to write to
     *  \param[in]  df   discrete function to write
     *
     *  \returns the STL stream (for concatenation)
     */
    template< class Impl >
    inline std :: ostream &
      operator<< ( std :: ostream &out,
                   const DiscreteFunctionInterface< Impl > &df )
    {
      df.print( out );
      return out;
    }



    /** \brief write a discrete function into an output stream
     *  \relates DiscreteFunctionInterface
     *  \relatesalso OutStreamInterface
     *
     *  \param[in]  out  stream to write to
     *  \param[in]  df   discrete function to write
     *
     *  \returns the output stream (for concatenation)
     */
    template< class StreamTraits, class Impl >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out,
                   const DiscreteFunctionInterface< Impl > &df )
    {
      df.write( out );
      return out;
    }



    /** \brief read a discrete function from an input stream
     *  \relates DiscreteFunctionInterface
     *  \relatesalso InStreamInterface
     *
     *  \param[in]   in  stream to read from
     *  \param[out]  df  discrete function to read
     *
     *  \returns the input stream (for concatenation)
     */
    template< class StreamTraits, class Impl >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in,
                   DiscreteFunctionInterface< Impl > &df )
    {
      df.read( in );
      return in;
    }

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_DISCRETEFUNCTION_INLINE_HH
