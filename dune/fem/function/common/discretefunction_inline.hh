#ifndef DUNE_FEM_DISCRETEFUNCTION_INLINE_HH
#define DUNE_FEM_DISCRETEFUNCTION_INLINE_HH

#include <fstream>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/gridpart/common/persistentindexset.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

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
                                   const DiscreteFunctionSpaceType &dfSpace,
                                   const LocalDofVectorAllocatorType &ldvAllocator )
    : dfSpace_( dfSpace ),
      ldvAllocator_( ldvAllocator ),
      name_( name ),
      scalarProduct_( dfSpace )
    {
    }


    template< class Impl >
    inline void DiscreteFunctionDefault< Impl > :: clear ()
    {
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it = 0;
    }


    template< class Impl >
    inline typename DiscreteFunctionDefault< Impl > :: RangeFieldType *
    DiscreteFunctionDefault< Impl > :: allocDofPointer()
    {
      dofPointerLock_.lock();

      const unsigned int size = BaseType :: size();
      RangeFieldType *dofPointer = new RangeFieldType[ size ];

      unsigned int i = 0;
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        dofPointer[ i++ ] = *it;
      assert( i == size );

      return dofPointer;
    }


    template< class Impl >
    inline void DiscreteFunctionDefault< Impl >
      :: freeDofPointer( RangeFieldType *dofPointer )
    {
      unsigned int i = 0;
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it = dofPointer[ i++ ];
      assert( i == BaseType :: size() );

      delete[] dofPointer;
      dofPointerLock_.unlock();
    }


    template< class Impl >
    inline void DiscreteFunctionDefault< Impl >
      ::axpy ( const RangeFieldType &s, const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );
      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        (*it) += s * (*git);
    }


    template< class Impl >
    inline void DiscreteFunctionDefault<Impl >
      :: print ( std::ostream &out ) const
    {
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType dit = BaseType :: dbegin(); dit != end; ++dit )
        out << (*dit) << std::endl;
    }


    template< class Impl >
    inline bool DiscreteFunctionDefault< Impl >
      :: dofsValid () const
    {
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        if( *it != *it )
          return false;

      return true;
    }


    template< class Impl >
    inline void DiscreteFunctionDefault< Impl >
      ::assign ( const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );

      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        *it = *git;
    }


    template< class Impl >
    template< class Operation >
    inline typename DiscreteFunctionDefault< Impl >
      :: template CommDataHandle< Operation > :: Type
    DiscreteFunctionDefault< Impl > :: dataHandle ( const Operation *operation )
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

      const EntityType entity(entitySearch( x ));
      const typename EntityType::Geometry geometry = entity.geometry();
      functor( geometry.local( x ), BaseType::localFunction( entity ) );
    }


    template< class Impl >
    inline typename DiscreteFunctionDefault< Impl > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Impl >
      ::operator+= ( const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );

      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        *it += *git;
      return asImp();
    }


    template< class Impl >
    template< class DFType >
    inline typename DiscreteFunctionDefault< Impl > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Impl >
      :: operator-= ( const DFType &g )
    {
      assert( BaseType :: size() == g.size() );

      const DofIteratorType end = BaseType :: dend();
      typename DFType :: ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
        *it -= *git;
      return asImp();
    }


    template< class Impl >
    inline typename DiscreteFunctionDefault< Impl > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Impl >
      :: operator*= ( const RangeFieldType &scalar )
    {
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it *= scalar;
      return asImp();
    }


    template< class Impl >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Impl >
      :: read ( InStreamInterface< StreamTraits > &in )
    {
      unsigned int versionId = in.readUnsignedInt();
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
        const DFSpaceIdentifier spaceId = space().type();
        int mySpaceIdInt;
        in >> mySpaceIdInt;
        const DFSpaceIdentifier mySpaceId = (DFSpaceIdentifier) mySpaceIdInt;

        if( spaceId != mySpaceId )
          DUNE_THROW( IOError, "Trying to read discrete function from different space: DFSpace (" << spaceName( spaceId ) << ") != DataSpace (" << spaceName( mySpaceId ) << ")" );
      }

      // read name
      in >> name_;

      // read size as integer
      int mysize;
      in >> mysize ;

      // check size
      if( mysize != BaseType :: size() &&
          BaseType :: size() != this->space().size() ) // only read compressed vectors
      {
        DUNE_THROW( IOError, "Trying to read discrete function of different size." );
      }

      // read all dofs
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        in >> *it;
    }


    template< class Impl >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Impl >
      :: write ( OutStreamInterface< StreamTraits > &out ) const
    {
      unsigned int versionId = DUNE_MODULE_VERSION_ID(DUNE_FEM);
      out << versionId ;

      // write space id to for testing when function is read
      int spaceId = space().type();
      out << spaceId ;

      // write name
      out << name_;

      // only allow write when vector is compressed
      if( BaseType :: size() != this->space().size() )
        DUNE_THROW(InvalidStateException,"Writing DiscreteFunction in uncompressed state!");

      // write size as integer
      const int mysize = BaseType :: size();
      out << mysize;

      // write all dofs
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        out << *it;
    }


    template< class Impl >
    void DiscreteFunctionDefault< Impl >
      :: insertSubData()
    {
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      IndexSetType &indexSet = (IndexSetType &)space().indexSet();
      if( Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::v )
      {
        PersistentIndexSetInterface* persistentIndexSet
          = Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::map( indexSet );

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
      IndexSetType &indexSet = (IndexSetType &)space().indexSet();
      if( Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::v )
      {
        PersistentIndexSetInterface* persistentIndexSet
          = Dune::Fem::Capabilities::isPersistentIndexSet< IndexSetType >::map( indexSet );

        // this unmarks the index set in the DofManager's list of index set as persistent
        if( persistentIndexSet )
          persistentIndexSet->removeBackupRestore();
      }
    }


    template< class Impl >
    inline bool DiscreteFunctionDefault< Impl >
      :: operator== ( const DiscreteFunctionType &g ) const
    {
      if( BaseType :: size() != g.size() )
        return false;

      const ConstDofIteratorType end = BaseType :: dend();

      ConstDofIteratorType fit = BaseType :: dbegin();
      ConstDofIteratorType git = g.dbegin();
      for( ; fit != end; ++fit, ++git )
        if( *fit != *git )
          return false;

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
