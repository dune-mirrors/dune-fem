namespace Dune
{

  template< class BaseFunctionImp >
  template< class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunctionImp >
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



  template< class BaseFunctionImp >
  template< class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunctionImp >
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



  template< class BaseFunctionImp > 
  template< class MatrixOfflineType, class MatrixOnlineType >
  inline void ReducedBasisSpace< BaseFunctionImp >
    :: restrictMatrix ( const MatrixOfflineType &matrixOffline,
                        MatrixOnlineType &matrixOnline ) const
  {
    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (matrixOnline.cols() == size)
            && (matrixOnline.rows == size) );
    assert( (matrixOffline.cols() == baseFunctionSpace_.size()) 
            && (matrixOffline.rows() == baseFunctionSpace_.size()) );

    BaseFunctionType Sphi( baseFunctionList_[ 0 ] );
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



  template< class BaseFunctionImp >
  template< class StreamTraits >
  inline void ReducedBasisSpace< BaseFunctionImp >
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

  

  template< class BaseFunctionImp >
  template< class StreamTraits >
  inline void ReducedBasisSpace< BaseFunctionImp >
    :: write ( OutStreamInterface< StreamTraits > &out )
  {
    const unsigned int size = numBaseFunctions();

    out << size;
    for( unsigned int i = 0; i < size; ++i )
      out << baseFunction( i );
  }

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
