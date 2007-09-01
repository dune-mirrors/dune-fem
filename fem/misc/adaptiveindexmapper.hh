#ifndef DUNE_FEM_ADAPTIVEINDEXMAPPER_HH
#define DUNE_FEM_ADAPTIVEINDEXMAPPER_HH

#include <dune/fem/storage/array.hh>

namespace Dune
{

  class IndexMapperHole
  {
  private:
    typedef IndexMapperHole ThisType;

  protected:
    unsigned int oldIndex_;
    unsigned int newIndex_;

  public:
    inline IndexMapperHole ()
    {
    }
    
    inline IndexMapperHole ( const ThisType &other )
    : oldIndex_( other.oldIndex_ ),
      newIndex_( other.newIndex_ )
    {
    }
    
    inline ThisType &operator= ( const ThisType &other )
    {
      oldIndex_ = other.oldIndex_;
      newIndex_ = other.newIndex_;
      return *this;
    }
    
    inline unsigned int newIndex () const
    {
      return newIndex_;
    }

    inline unsigned int oldIndex () const
    {
      return oldIndex_;
    }

    inline void setNewIndex ( const unsigned int i )
    {
      newIndex_ = i;
    }

    inline void setOldIndex ( const unsigned int i )
    {
      oldIndex_ = i;
    }
  };



  class AdaptiveIndexMapper
  {
  private:
    typedef AdaptiveIndexMapper ThisType;

  protected:
    enum IndexState { Unused = -1, Used = 1, New = 2, Deleted = -2 };

    typedef DynamicArray< unsigned int, DefaultArrayOverAllocator > IndexArrayType;
    typedef DynamicArray< char, DefaultArrayOverAllocator > StateArrayType;

  public:
    typedef IndexMapperHole HoleType;
    typedef DynamicArray< HoleType > HoleArrayType;

  protected:
    IndexArrayType index_;

    StateArrayType state_;

    unsigned int size_;

  public:
    inline AdaptiveIndexMapper ()
    : size_( 0 )
    {
    }

  private:
    // prohibit copying
    inline AdaptiveIndexMapper ( const ThisType & );

    // prohibit copying
    inline ThisType &operator= ( const ThisType & );

  public:
    /** \brief find out, whether an index is actually mapped
     *  
     *  \note This method returns true for indices that are
     *        marked deleted. This is needed during adaption.
     *
     *  \param[in]  i  index which shall be mapped
     *
     *  \return true, if the index can still be mapped
     */
    inline bool contains ( unsigned int i ) const
    {
      return (state_[ i ] != Unused );
    }

    /** \brief return the size of the domain
     */
    inline unsigned int domainSize () const
    {
      return index_.size();
    }
    
    /** \brief map an index
     */
    inline unsigned int index ( unsigned int i ) const
    {
      assert( this->contains( i ) );
      return index_[ i ];
    }

    /** \brief obtain the number of used indices
     */
    inline unsigned int size () const
    {
      return size_;
    }

    // Methods for adaption
    // --------------------

    /** \brief remove all indices
     */
    inline void clear ()
    {
      state_.assign( Unused );
      size_ = 0;
    }
   
    /** \brief compress the index map
     *
     *  After adding and removing indices, compress must be called. This method
     *  makes the index map continuous again and returns the necessary move
     *  operations to adapt a storage array to the new map.
     *  
     *  \param[in]  holes  array to hold the hole mappings
     *
     *  \returns true, if at least one hole was closed
     */
    inline bool compress ( HoleArrayType &holes );

    /** \brief resize the domain of the map
     *
     *  Before adapting the indices, this method is called to suitably resize
     *  the domain of the index map.
     */
    inline void resize ( const unsigned int newSize )
    {
      index_.resize( newSize );
      state_.resize( newSize, Unused );
    }

    /** \brief insert an index into the map
     *
     *  \note The index is marked as new, so that later calls to remove do not
     *        remove this index
     *
     *  \param[in]  i  index to insert
     */
    inline void insert ( const unsigned int i )
    {
      assert( i < domainSize() );

      if( state_[ i ] == Unused )
        index_[ i ] = size_++;
      state_[ i ] = New;
    }

    /** \brief remove an index from the map
     *
     *  \note Newly added indices are not removed by this method.
     *
     *  \param[in]  i  index to remove
     */
    inline void remove ( const unsigned int i )
    {
      assert( i < domainSize() );
      if( state_[ i ] == Used )
        state_[ i ] = Deleted;
    }
  };
  
}

#include "adaptiveindexmapper_inline.hh"

#endif
