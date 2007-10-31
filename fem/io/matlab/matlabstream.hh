#ifndef DUNE_FEM_MATLABSTREAM_HH
#define DUNE_FEM_MATLABSTREAM_HH

#include <dune/grid/common/gridpart.hh>

#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/function/common/discretefunction.hh>


namespace Dune
{

  class MatlabOutStream
  {
  protected:
    XDRFileOutStream stream_;

  public:
    inline MatlabOutStream ( const std :: string &filename )
    : stream_( filename )
    {
    }

    inline void write ( const double &value )
    {
      stream_ << value;
    }

    inline void write ( const float &value )
    {
      stream_ << value;
    }

    inline void write ( const int &value )
    {
      stream_ << value;
    }

    inline void write ( const unsigned int &value )
    {
      stream_ << value;
    }

    template< class Traits >
    inline void write ( const VectorInterface< Traits > &v )
    {
      const unsigned int size = v.size();

      write( size );
      for( unsigned int i = 0; i < size; ++i )
        write( v[ i ] );
    }

    template< class T >
    inline void write ( const DenseMatrix< T > &matrix )
    {
      const unsigned int rows = matrix.rows();
      const unsigned int cols = matrix.cols();

      write( rows );
      write( cols );
      for( unsigned int i = 0; i < rows; ++i )
      {
        for( unsigned int j = 0; j < cols; ++j )
          write( matrix[ i ][j] );
      }
    }

    template< class Traits >
    inline void write ( const GridPartInterface< Traits > &gridPart )
    {
      typedef GridPartInterface< Traits > GridPartType;
      typedef typename GridPartType :: IndexSetType IndexSetType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename Traits :: template Codim< 0 > :: IteratorType IteratorType;

      enum { dim = GridType :: dimension };
      enum { dimworld = GridType :: dimensionworld };

      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename GridType :: template Codim< dim > :: EntityPointer
        VertexPtrType;
      typedef typename GridType :: template Codim< dim > :: Entity VertexType;
      typedef FieldVector< typename GridType :: ctype, dimworld > PointType;

      const IndexSetType &indexSet = gridPart.indexSet();
      const unsigned int numEntities = indexSet.size( 0 );
      const unsigned int numVertices = indexSet.size( dim );

      DenseMatrix< typename GridType :: ctype > points( dimworld, numVertices );
      DenseMatrix< unsigned int > triangles( dim+1, numEntities );

      const IteratorType end = gridPart.template end< 0 >();
      for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
      {
        const EntityType &entity = *it;

        if( entity.template count< dim >() != dim+1 )
          DUNE_THROW( IOError, "MatlabOutStream: Cannot write non-simplex grid." );

        const unsigned int enindex = indexSet.index( entity );
        // run over all vertices of the simplex 
        for( unsigned int i = 0; i <= dim; ++i )
        {
          const VertexPtrType vertexptr = entity.template entity< dim >( i );
          const VertexType &vertex = *vertexptr;

          const PointType &point = vertex.geometry()[ 0 ];
          const unsigned int ptindex = indexSet.index( vertex );

          for( unsigned int j = 0; j < dimworld; ++j )
            points[ j ][ ptindex ] = point[ j ];
          
          triangles[ i ][ enindex ] = ptindex + 1;
        }
      }

      write( points );
      write( triangles );
    }

  };

  template< class T >
  inline MatlabOutStream &operator<< ( MatlabOutStream &out,
                                       const T &value )
  {
    out.write( value );
    return out;
  }

}

#endif
