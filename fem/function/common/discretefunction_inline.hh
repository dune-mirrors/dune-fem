#ifndef DUNE_DISCRETEFUNCTION_INLINE_HH
#define DUNE_DISCRETEFUNCTION_INLINE_HH

#include <fstream>

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/streams/streams.hh>

#include "discretefunction.hh"

namespace Dune 
{
 
  // DiscreteFunctionDefault 
  // -----------------------

  template< class Traits >
  inline DiscreteFunctionDefault< Traits >
    :: DiscreteFunctionDefault ( const std :: string &name,
                                 const DiscreteFunctionSpaceType &dfSpace,
                                 const LocalFunctionFactoryType &lfFactory )
  : DiscreteFunctionInterfaceType( dfSpace ),
    lfStorage_( lfFactory ),
    name_( name ),
    scalarProduct_( dfSpace )
  {}

  
  template< class Traits >
  inline DiscreteFunctionDefault< Traits >
    :: ~DiscreteFunctionDefault ()
  {
    assert( !dofPointerLock_ );
  }


  template< class Traits >
  inline const std :: string &DiscreteFunctionDefault< Traits > :: name () const
  {
    return name_;
  }


  template< class Traits >
  template< class EntityType >
  inline const typename DiscreteFunctionDefault< Traits > ::  LocalFunctionType
  DiscreteFunctionDefault< Traits >
    :: localFunction ( const EntityType &entity ) const
  {
    return lfStorage_.localFunction( entity );
  }


  template< class Traits >
  template< class EntityType >
  inline typename DiscreteFunctionDefault< Traits > ::  LocalFunctionType
  DiscreteFunctionDefault< Traits >
    :: localFunction ( const EntityType &entity )
  {
    return lfStorage_.localFunction( entity );
  }


  template< class Traits >
  inline void DiscreteFunctionDefault< Traits > :: clear ()
  {
    const DofIteratorType end = BaseType :: dend();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
      *it = 0;
  }

  
  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: RangeFieldType *
  DiscreteFunctionDefault< Traits > :: allocDofPointer()
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

  
  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
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


  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: addScaled ( const DiscreteFunctionType &g,
                   const RangeFieldType &s )
  {
    assert( BaseType :: size() == g.size() );
    const DofIteratorType end = BaseType :: dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
      (*it) += s * (*git);
  }


  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: RangeFieldType
  DiscreteFunctionDefault< Traits >
    :: scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const
  {
    return scalarProduct_.scalarProductDofs( *this, other );
  }

  
  template< class Traits >
  inline void DiscreteFunctionDefault<Traits >
    :: print ( std::ostream &out ) const
  {
    out << BaseType :: name() << std::endl;
    
    const ConstDofIteratorType end = BaseType :: dend();
    for( ConstDofIteratorType dit = BaseType :: dbegin(); dit != end; ++dit )
      out << (*dit) << std::endl;
  }


  template< class Traits >
  inline bool DiscreteFunctionDefault< Traits >
    :: dofsValid () const
  {
    const ConstDofIteratorType end = BaseType :: dend();
    for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
      if( *it != *it )
        return false;

    return true;
  }

  
  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: assign( const DiscreteFunctionType &g )
  {
    assert( BaseType :: size() == g.size() );

    const DofIteratorType end = BaseType :: dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
      *it = *git;
  }


  template< class Traits >
  template< class Operation >
  inline typename DiscreteFunctionDefault< Traits >
    :: template CommDataHandle< Operation > :: Type
  DiscreteFunctionDefault< Traits > :: dataHandle ( const Operation *operation )
  {
    return BaseType :: space().createDataHandle( asImp(), operation );
  }


  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: evaluate ( const DomainType &x,
                  RangeType &ret ) const
  {
    FieldVector< deriType, 0 > diffVariable;
    BaseType :: evaluate( diffVariable, x, ret );
  }


  template< class Traits >
  template< int diffOrder >
  inline void DiscreteFunctionDefault< Traits >
    :: evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                  const DomainType &x,
                  RangeType &ret ) const
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
    
    const DiscreteFunctionSpaceType &space = BaseType :: space();
    const IteratorType end = space.end();
    for( IteratorType it = space.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      const DomainType xlocal = geometry.local( x );
      if( geometry.checkInside( xlocal ) )
      {
        const LocalFunctionType localFunction
          = BaseType :: localFunction( entity );

        localFunction.evaluate( diffVariable, xlocal, ret );
        return;
      }
    }
    DUNE_THROW( RangeError,
                "DiscreteFunctionDefault :: evaluate: x is within domain." );
  }


  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
  DiscreteFunctionDefault< Traits >
    :: operator+= ( const DiscreteFunctionType &g )
  {
    assert( BaseType :: size() == g.size() );

    const DofIteratorType end = BaseType :: dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
      *it += *git;
    return asImp();
  }


  template< class Traits >
  template< class DFType >
  inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
  DiscreteFunctionDefault< Traits >
    :: operator-= ( const DFType &g )
  {
    assert( BaseType :: size() == g.size() );

    const DofIteratorType end = BaseType :: dend();
    typename DFType :: ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
      *it -= *git;
    return asImp();
  }


  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
  DiscreteFunctionDefault< Traits >
    :: operator*= ( const RangeFieldType &scalar )
  {
    const DofIteratorType end = BaseType :: dend();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
      *it *= scalar;
    return asImp();
  }


  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
  DiscreteFunctionDefault< Traits >
    :: operator/= ( const RangeFieldType &scalar )
  {
    return BaseType :: operator*=( RangeFieldType( 1 ) / scalar );
  }


  template< class Traits >
  template< class StreamTraits >
  inline void DiscreteFunctionDefault< Traits >
    :: read ( InStreamInterface< StreamTraits > &in )
  {
#if ! DUNE_FEM_COMPATIBILITY
    unsigned int versionId = in.readUnsignedInt();
    if( versionId < DUNE_VERSION_ID(0,9,1) )
      DUNE_THROW( IOError, "Trying to read outdated file." );
    else if( versionId > DUNE_MODULE_VERSION_ID(DUNE_FEM) )
      std :: cerr << "Warning: Reading discrete function from newer version: "
                  << versionId << std :: endl;

    in >> name_;
#endif
    
    if( in.readInt() != BaseType :: size() && 
        BaseType :: size() != this->space().size() ) // only read compressed vectors 
      DUNE_THROW( IOError, "Trying to read discrete function of different size." );

    const DofIteratorType end = BaseType :: dend();
    for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
      in >> *it;
  }


  template< class Traits >
  template< class StreamTraits >
  inline void DiscreteFunctionDefault< Traits >
    :: write ( OutStreamInterface< StreamTraits > &out ) const
  {
#if ! DUNE_FEM_COMPATIBILITY
    out << DUNE_MODULE_VERSION_ID(DUNE_FEM);
    out << name_;
#endif
  
    // only allow write when vector is compressed 
    if( BaseType :: size() != this->space().size() )
      DUNE_THROW(InvalidStateException,"Writing DiscreteFunction in uncompressed state!");
    
    out << BaseType :: size();

    const ConstDofIteratorType end = BaseType :: dend();
    for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
      out << *it;
  }

  template< class Traits >
  bool DiscreteFunctionDefault< Traits >
    :: read_xdr ( const std :: string filename )
  {
    try
    {
      XDRFileInStream in( filename );
      BaseType :: read( in );
      return true;
    }
    catch( Exception e )
    {
      return false;
    }
  }
 
  
  template< class Traits >
  bool DiscreteFunctionDefault< Traits >
    :: write_xdr ( const std :: string filename ) const
  {
    try
    {
      XDRFileOutStream out( filename );
      BaseType :: write( out );
      return true;
    }
    catch( Exception e )
    {
      return false;
    }
  }
  
  
  template< class Traits >
  bool DiscreteFunctionDefault< Traits >
    :: read_ascii ( const std :: string filename )
  {
    try
    {
      ASCIIInStream in( filename );
      BaseType :: read( in );
      return true;
    }
    catch( Exception e )
    {
      return false;
    }
  }
 
  
  template< class Traits >
  bool DiscreteFunctionDefault< Traits >
    :: write_ascii ( const std :: string filename ) const
  {
    try
    {
      ASCIIOutStream out( filename );
      BaseType :: write( out );
      return true;
    }
    catch( Exception e )
    {
      return false;
    }
  }

  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: enableDofCompression ()
  {}


  template< class Traits >
  inline bool DiscreteFunctionDefault< Traits >
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
  
 
  template< class Traits >
  inline bool DiscreteFunctionDefault< Traits >
    :: operator!= ( const DiscreteFunctionType &g ) const
  {
    return !(operator==( g ));
  }

  
  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: LocalFunctionStorageType &
  DiscreteFunctionDefault< Traits > :: localFunctionStorage () const
  {
    return lfStorage_;
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
  template< class Traits >
  inline std :: ostream &
    operator<< ( std :: ostream &out,
                 const DiscreteFunctionInterface< Traits > &df )
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
  template< class StreamTraits, class Traits >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const DiscreteFunctionInterface< Traits > &df )
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
  template< class StreamTraits, class Traits >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 DiscreteFunctionInterface< Traits > &df )
  {
    df.read( in );
    return in;
  }

} // end namespace Dune
#endif
