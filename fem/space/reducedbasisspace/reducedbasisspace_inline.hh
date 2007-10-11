namespace Dune
{

  template< class StreamTraits, class BaseFunctionType >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const ReducedBasisSpace< BaseFunctionType > &space )
  {
    const unsigned int size = space.numBaseFunctions();

    out << size;
    for( unsigned int i = 0; i < size; ++i )
      out << space.baseFunction( i );

    return out;
  }

  template< class StreamTraits, class BaseFunctionType >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.clear();
    
    unsigned int size;
    in >> size;
    
    BaseFunctionType baseFunction( "Base Function", space.baseFunctionSpace() );
    for( unsigned int i = 0; i < size; ++i )
    {
      in >> baseFunction;
      space.addBaseFunction( baseFunction );
    }
    
    return in;
  }
  
}
