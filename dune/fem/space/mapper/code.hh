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

      /**
       * \brief execute DoF mapper code
       *
       * The functor has to be a copyable object satisfying the following
       * interface:
       * \code
       * struct Functor
       * {
       *   void operator() ( unsigned int gtIndex, unsigned int subEntity, ConstIterator begin, ConstIterator end )
       * };
       * \endcode
       * The type ConstIterator is defined by the DofMapperCode.
       * It is passed the following arguments:
       *   - gtIndex: global geometry type index of a subentity
       *   - subEntity: local number of the subentity (wrt. the reference element)
       *   - begin / end: iterator pair returning the local indices (wrt. to the
       *     element) of the DoFs associated to the subentity
       *   .
       */
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

      /**@internal
       *
       * @param[in] numBlocks The number of subentities which carry DoFs
       *
       * @param[in] numDofs The total number of DoFs
       */
      static std::size_t size ( unsigned int numBlocks, unsigned int numDofs )
      {
        return 2 + 3*numBlocks + numDofs;
      }

      // Format of the code_ array:
      //
      // code_[0]: number of subentities with DoFs
      // code_[1]: total number of DoFs
      //
      // For all k = 0 ... (numBlocks-1)
      // (NB: k corresponds to a subentity with DoFs)
      // It follows a variable size block (offset_0 := 2):
      //
      // code_[offset_k + 0]: global geometry type index of the subentity
      //                      (this also encodes the codim)
      // code_[offset_k + 1]: local number i_k of subentity (0 <= i_k < refElem.size( codim ))
      // code_[offset_k + 2]: #DoFs n_k attached to this subentity
      //
      // code_[offset_k + 3 + j]:
      // for 0 <= j < n_k, the local index of the given DoF, where "local" now
      // means the number of the corresponding local basis function of the bulk
      // element, i.e., not the numbering inside the entity.
      //
      // offset_(k+1) := (offset_k + 3 + n_k) is then just the start of the
      // next block.
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
        assert( (std::ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }

      unsigned int &operator[] ( unsigned int i )
      {
        assert( (std::ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_CODE_HH
