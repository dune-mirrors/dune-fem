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
    :: restrictVector ( const BaseFunctionType &sourceFunction,
                 DiscreteFunctionType &destFunction ) const
  {
    typedef typename DiscreteFunctionType :: RangeFieldType DofType;

    const unsigned int size = baseFunctionList_.size();

    destFunction.clear();
    for( unsigned int i = 0; i < size; ++i )
    {
      const BaseFunctionType &baseFunction = *(baseFunctionList_[ i ]);
      destFunction.addScaled( baseFunction, sourceFunction );
    }
  }

  
  template< class BaseFunctionImp > 
  template< class MatrixOnlineType , class MatrixOfflineType >
  inline void ReducedBasisSpace< BaseFunctionImp >
    :: restrictMatrix ( const MatrixOnlineType &matrixOnline,
                 MatrixOfflineType &matrixOffline ) const
  {

    const unsigned int size = baseFunctionList_.size();
    assert( size == matrixOffline.cols() );
    assert( matrixOnline.cols() == baseFunctionSpace.size());
    matrixOffline.clear();
    BaseFunctionType Sphi("Sphi", baseFunctionSpace_);

    for( unsigned int i = 0; i < size; ++i )
    {
      const BaseFunctionType &baseFunction = *(baseFunctionList_[ i ]);
      matrixOnline.mult(baseFunction, Sphi);

      for( unsigned int j = 0; j < size; ++j)
      {
      const BaseFunctionType &baseFunctionj = *(baseFunctionList_[ j ]);
      double val = Sphi * baseFunctionj ;
      matrixOffline.set(i,j,val);
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

  
  template< class StreamTraits, class BaseFunctionType >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.write( out );
    return out;
  }

  template< class StreamTraits, class BaseFunctionType >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.read( in ); 
    return in;
  }
  
}
