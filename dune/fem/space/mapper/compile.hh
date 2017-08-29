#ifndef DUNE_FEM_DOFMAPPER_COMPILE_HH
#define DUNE_FEM_DOFMAPPER_COMPILE_HH

#include <algorithm>
#include <type_traits>
#include <utility>

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

    template< class RefElement,
              std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().size( 0 ) ) >, int >::value, int > = 0,
              std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().type( 0, 0 ) ) >, GeometryType >::value, int > = 0 >
    inline DofMapperCode generateCodimensionCode ( const RefElement &refElement, int codim, unsigned int blockSize = 1 )
    {
      unsigned int count = refElement.size( codim );
      DofMapperCodeWriter code( count, count*blockSize );
      unsigned int pos = 0;
      for( unsigned int i = 0; i < count; ++i )
      {
        code[ pos++ ] = GlobalGeometryTypeIndex::index( refElement.type( i, codim ) );
        code[ pos++ ] = i;
        code[ pos++ ] = blockSize;
        for( unsigned int j = 0; j < blockSize; ++j )
          code[ pos++ ] = i*blockSize + j;
      }
      return code;
    }



    // compile (for LocalCoefficients)
    // -------------------------------

    template< class RefElement, class LocalCoefficients,
              std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().size( 0 ) ) >, int >::value, int > = 0,
              std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().type( 0, 0 ) ) >, GeometryType >::value, int > = 0 >
    inline DofMapperCode compile ( const RefElement &refElement, const LocalCoefficients &localCoefficients )
    {
      const int dim = RefElement::dimension;

      const std::size_t numDofs = localCoefficients.size(); // total number of DoFs

      // count number of keys per subentity

      // total number of all sub-entities
      unsigned int numSubEntities = 0;
      for( int codim = 0; codim <= dim; ++codim )
        numSubEntities += refElement.size( codim );
      assert( numSubEntities > 0 );

      // form a "matrix" with variable lenght rows. This is the usual
      // approach: pre-allocate the needed storage once and then
      // insert the proper offsets into the row-pointer. After
      // completion count[codim] is an array with one entry for each
      // sub-entity for the given codim. It is initialized with zeros.
      unsigned int *count[ dim+1 ];
      count[ 0 ] = new unsigned int[ numSubEntities ];
      assert( count[ 0 ] );
      std::fill( count[ 0 ], count[ 0 ] + numSubEntities, (unsigned int)0 );
      for( int codim = 0; codim < dim; ++codim )
        count[ codim+1 ] = count[ codim ] + refElement.size( codim );

      // Now do the actual counting. After completion
      // cound[codim][subEntity] will contain the number of DoFs
      // attached to the particular sub-entity.
      //
      // numBlocks is the actual number of __USED__
      // sub-entities. E.g. for continuous Lagrange-1 on a triangle numBlocks
      // would be 3, after counting (only the vertices carry DoFs).
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
      //
      // I.e.: count[cd][subEntIdx] = local index offset for start of
      // DoFs attached to sub entity

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
      //
      // Format of the code-array is described in code.hh
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
