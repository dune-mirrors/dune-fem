#ifndef DUNE_FEM_DOFMAPPER_CODE_HH
#define DUNE_FEM_DOFMAPPER_CODE_HH

#include <algorithm>
#include <cassert>
#include <iostream>

namespace Dune
{

  namespace Fem
  {

    // DofMapperCode
    // -------------

    class DofMapperCode
    {
    protected:
      typedef const unsigned int *ConstIterator;
      typedef unsigned int *Iterator;

      DofMapperCode ( unsigned int numBlocks, unsigned int numDofs )
      {
        code_ = new unsigned int[ size( numBlocks, numDofs ) ];
        code_[ 0 ] = numBlocks;
        code_[ 1 ] = numDofs;
      }

    public:
      DofMapperCode ()
      {
        code_ = new unsigned int[ size( 0, 0 ) ];
        code_[ 1 ] = code_[ 0 ] = 0;
      }

      DofMapperCode ( const DofMapperCode &other )
      {
        code_ = new unsigned int[ other.size() ];
        std::copy( ConstIterator( other.code_ ), other.end(), code_ );
      }

      ~DofMapperCode ()
      {
        delete[] code_;
      }

      const DofMapperCode &operator= ( const DofMapperCode &other )
      {
        if( size() != other.size() )
        {
          delete[] code_;
          code_ = new unsigned int[ other.size() ];
        }
        std::copy( ConstIterator( other.code_ ), other.end(), code_ );
        return *this;
      }

      template< class Functor >
      void operator() ( Functor f ) const
      {
        for( ConstIterator it = begin(); it != end(); )
        {
          const unsigned int gtIndex = *(it++);
          const unsigned int subEntity = *(it++);
          unsigned int nDofs = *(it++);
          f( gtIndex, subEntity, ConstIterator( it ), ConstIterator( it + nDofs ) );
          it += nDofs;
        }
      }

      unsigned int numBlocks () const { return code_[ 0 ]; }
      unsigned int numDofs () const { return code_[ 1 ]; }

      friend std::ostream &operator<< ( std::ostream &out, const DofMapperCode &code )
      {
        out << "{ " << code.numBlocks() << ", " << code.numDofs();
        for( DofMapperCode::ConstIterator it = code.begin(); it != code.end(); ++it )
          out << ", " << *it;
        return out << " }";
      }

    protected:
      ConstIterator begin () const { return code_ + 2; }
      Iterator begin () { return code_ + 2; }
      ConstIterator end () const { return code_ + size(); }
      Iterator end () { return code_ + size(); }

      std::size_t size () const { return size( numBlocks(), numDofs() ); }

      static std::size_t size ( unsigned int numBlocks, unsigned int numDofs )
      {
        return 2 + 3*numBlocks + numDofs;
      }

      unsigned int *code_;
    };



    // DofMapperCodeWriter
    // -------------------

    class DofMapperCodeWriter
    : public DofMapperCode
    {
    public:
      DofMapperCodeWriter ( unsigned int numBlocks, unsigned int numDofs )
      : DofMapperCode( numBlocks, numDofs )
      {}

      const unsigned int &operator[] ( unsigned int i ) const
      {
        assert( (ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }

      unsigned int &operator[] ( unsigned int i )
      {
        assert( (ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_CODE_HH
