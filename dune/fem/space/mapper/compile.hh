#ifndef DUNE_FEM_DOFMAPPER_COMPILE_HH
#define DUNE_FEM_DOFMAPPER_COMPILE_HH

#include <algorithm>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/localkey.hh>

namespace Dune
{

  namespace Fem
  {

    // generateCodimensionCode
    // -----------------------

    template< class Field, int dim >
    inline DofMapperCode
    generateCodimensionCode ( const Dune::ReferenceElement< Field, dim > &refElement, int codim )
    {
      unsigned int count = refElement.size( codim );
      DofMapperCodeWriter code( count, count );
      for( unsigned int i = 0; i < count; ++i )
      {
        code[ 4*i+0 ] = GlobalGeometryTypeIndex::index( refElement.type( i, codim ) );
        code[ 4*i+1 ] = i;
        code[ 4*i+2 ] = 1;
        code[ 4*i+3 ] = i;
      }
      return code;
    }



    // compile (for LocalCoefficients)
    // -------------------------------

    template< class Field, int dim, class LocalCoefficients >
    inline DofMapperCode
    compile ( const Dune::ReferenceElement< Field, dim > &refElement,
              const LocalCoefficients &localCoefficients )
    {
      const std::size_t numDofs = localCoefficients.size();

      // count number of keys per subentity

      unsigned int numSubEntities = 0;
      for( int codim = 0; codim <= dim; ++codim )
        numSubEntities += refElement.size( codim );
      assert( numSubEntities > 0 );

      unsigned int *count[ dim+1 ];
      count[ 0 ] = new unsigned int[ numSubEntities ];
      assert( count[ 0 ] );
      std::fill( count[ 0 ], count[ 0 ] + numSubEntities, (unsigned int)0 );
      for( int codim = 0; codim < dim; ++codim )
        count[ codim+1 ] = count[ codim ] + refElement.size( codim );

      unsigned int numBlocks = 0;
      for( std::size_t i = 0; i < numDofs; ++i )
      {
        const LocalKey &key = localCoefficients.localKey( i );

        const int codim = key.codim();
        const int subEntity = key.subEntity();

        assert( (codim >= 0) && (codim <= dim) );
        assert( (subEntity >= 0) && (subEntity < refElement.size( codim )) );

        if( count[ codim ][ subEntity ] == 0 )
          ++numBlocks;
        ++count[ codim ][ subEntity ];
      }

      // format the code into subentity blocks
      // result: count will hold the first local index in the block (0 = unused)

      DofMapperCodeWriter code( numBlocks, numDofs );

      unsigned int next = 0;
      for( int codim = 0; codim <= dim; ++codim )
      {
        for( int i = 0; i < refElement.size( codim ); ++i )
        {
          const unsigned int cnt = count[ codim ][ i ];
          if( cnt == 0 )
            continue;

          code[ next++ ] = GlobalGeometryTypeIndex::index( refElement.type( i, codim ) );
          code[ next++ ] = i;
          code[ next++ ] = cnt;

          count[ codim ][ i ] = next;
          next += cnt;
        }
      }

      // fill in the local indices

      for( std::size_t i = 0; i < numDofs; ++i )
      {
        const LocalKey &key = localCoefficients.localKey( i );
        const unsigned int block = count[ key.codim() ][ key.subEntity() ];
        assert( block > 0 );
        assert( (key.index() >= 0) && (key.index() < code[ block-1 ]) );
        code[ block + key.index() ] = i;
      }

      // clean up and return

      delete[] count[ 0 ];

      return code;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_COMPILE_HH
