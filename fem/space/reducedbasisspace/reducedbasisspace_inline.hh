#ifndef DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_INLINE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_INLINE_HH

#include "reducedbasisspace.hh"

namespace Dune
{

  // ReducedBasisSpace
  // -----------------

  template< class BaseFunction >
  inline ReducedBasisSpace< BaseFunction >
    :: ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace )
  : BaseType( baseFunctionSpace.gridPart() ),
    baseFunctionSpace_( baseFunctionSpace ),
    baseFunctionList_(),
    mapper_( baseFunctionList_ )
  {
  }

  
  template< class BaseFunction >
  template< class StreamTraits >
  inline ReducedBasisSpace< BaseFunction >
    :: ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                           InStreamInterface< StreamTraits > &in )
  : BaseType( baseFunctionSpace.gridPart() ),
    baseFunctionSpace_( baseFunctionSpace ),
    baseFunctionList_(),
    mapper_( baseFunctionList_ )
  {
    read( in );
  }


  template< class BaseFunction >
  inline ReducedBasisSpace< BaseFunction > :: ~ReducedBasisSpace ()
  {
    clear();
  }


  template< class BaseFunction >
  inline void ReducedBasisSpace< BaseFunction >
    :: addBaseFunction ( const BaseFunctionType &baseFunction )
  {
    BaseFunctionType *f = new BaseFunctionType( baseFunction );
    assert( f != NULL );
    baseFunctionList_.append( f );
  }


  template< class BaseFunction >
  inline const typename ReducedBasisSpace< BaseFunction > :: BaseFunctionType &
  ReducedBasisSpace< BaseFunction > :: baseFunction ( unsigned int i ) const
  {
    return *(baseFunctionList_[ i ]);
  }


  template< class BaseFunction >
  inline const typename ReducedBasisSpace< BaseFunction > :: BaseFunctionSpaceType &
  ReducedBasisSpace< BaseFunction > :: baseFunctionSpace () const
  {
    return baseFunctionSpace_;
  }


  template< class BaseFunction >
  inline void ReducedBasisSpace< BaseFunction > :: clear ()
  {
    crop( 0 );
  }


  template< class BaseFunction >
  inline void ReducedBasisSpace< BaseFunction > :: crop ( unsigned int n )
  {
    const unsigned int size = baseFunctionList_.size();
    assert( n <= size );
    for( unsigned int i = n; i < size; ++i )
      delete baseFunctionList_[ i ];
    baseFunctionList_.resize( n );
  }
 

  template< class BaseFunction >
  inline unsigned int ReducedBasisSpace< BaseFunction > :: numBaseFunctions () const
  {
    return baseFunctionList_.size();
  }

  
  template< class BaseFunction >
  template< class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunction >
    :: project ( const DiscreteFunctionType &sourceFunction,
                 BaseFunctionType &destFunction ) const
  {
    typedef typename DiscreteFunctionType :: RangeFieldType DofType;

    const unsigned int size = baseFunctionList_.size();

    destFunction.clear();
    for( unsigned int i = 0; i < size; ++i )
    {
      const BaseFunctionType &baseFunction = *(baseFunctionList_[ i ]);
      const DofType &dof = sourceFunction.dof( i );

      destFunction.addScaled( baseFunction, dof );
    }
  }


  template< class BaseFunction >
  template< class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunction >
    :: restrictFunction ( const BaseFunctionType &sourceFunction,
                          DiscreteFunctionType &destFunction ) const
  {
    const unsigned int size = baseFunctionList_.size();
    assert( size == destFunction.size() );

    for( unsigned int i = 0; i < size; ++i )
    {
      const BaseFunctionType &baseFunction = *(baseFunctionList_[ i ]);
      destFunction.dof( i )
        = baseFunction.scalarProductDofs( sourceFunction );
    }
  }


  template< class BaseFunction > 
  template< class MatrixOfflineType, class MatrixOnlineType >
  inline void ReducedBasisSpace< BaseFunction >
    :: restrictMatrix ( const MatrixOfflineType &matrixOffline,
                        MatrixOnlineType &matrixOnline ) const
  {
    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (matrixOnline.cols() == size)
            && (matrixOnline.rows() == size) );
    assert( (matrixOffline.cols() == baseFunctionSpace_.size()) 
            && (matrixOffline.rows() == baseFunctionSpace_.size()) );

    BaseFunctionType Sphi( *(baseFunctionList_[ 0 ]) );
    for( unsigned int i = 0; i < size; ++i )
    {
      const BaseFunctionType &phi = *(baseFunctionList_[ i ]);
      matrixOffline( phi, Sphi );

      for( unsigned int j = 0; j < size; ++j)
      {
        const BaseFunctionType &psi = *(baseFunctionList_[ j ]);
        matrixOnline.set( i, j, Sphi.scalarProductDofs( psi ) );
      }
    }
  }


  template< class BaseFunction >
  template< class StreamTraits >
  inline void ReducedBasisSpace< BaseFunction >
    :: read ( InStreamInterface< StreamTraits > &in )
  {
    clear();
    
    unsigned int size;
    in >> size;
    baseFunctionList_.resize( size );
   
    for( unsigned int i = 0; i < size; ++i )
    {
      BaseFunctionType *baseFunction
        = new BaseFunctionType( "BaseFunction", baseFunctionSpace() );
      in >> *baseFunction;
      baseFunctionList_[ i ] = baseFunction;
    }
  }


  template< class BaseFunction >
  template< class StreamTraits >
  inline void ReducedBasisSpace< BaseFunction >
    :: write ( OutStreamInterface< StreamTraits > &out ) const
  {
    const unsigned int size = numBaseFunctions();

    out << size;
    for( unsigned int i = 0; i < size; ++i )
      out << baseFunction( i );
  }



  // Stream Operators for ReducedBasisSpace
  // --------------------------------------

  /** \brief write a ReducedBasisSpace into an output stream
   *  \relates ReducedBasisSpace
   *  \relatesalso OutStreamInterface
   *
   *  \param[in]  out    stream to write to
   *  \param[in]  space  ReducedBasisSpace to write
   *
   *  \returns the output stream (for concatenation)
   */
  template< class StreamTraits, class BaseFunctionType >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.write( out );
    return out;
  }


  /** \brief read a ReducedBasisSpace from an input stream
   *  \relates ReducedBasisSpace
   *  \relatesalso InStreamInterface
   *
   *  \param[in]   in     stream to read from
   *  \param[out]  space  ReducedBasisSpace to read
   *
   *  \returns the input stream (for concatenation)
   */
  template< class StreamTraits, class BaseFunctionType >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.read( in ); 
    return in;
  }
  
}

#endif
