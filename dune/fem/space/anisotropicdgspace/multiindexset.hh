#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_MULTIINDEXSET_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_MULTIINDEXSET_HH

// C++ includes
#include <cstddef>

// dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Provides an iterator over all multi indices in an AnisotropicDGSpace
*/


namespace AnisotropicDG
{

  // MultiIndexSet
  // -------------

  template< int dimension, int order >
  class MultiIndexSet
  {
    class Iterator;

  public:
    //! \brief type of multi index
    typedef Dune::FieldVector< int, dimension > MultiIndexType; 
    //! \brief return
    typedef Iterator IteratorType;

    //! \brief return number of multi indices
    static std::size_t size ()
    {
      return Dune::StaticPower< order+1, dimension >::power;
    }

    //! \brief return true, if index is contained in this set
    static bool contains ( const MultiIndexType &multiIndex )
    {
      for( int i = 0; i < dimension; ++i )
      {
        const int k = multiIndex[ i ];
        if( k < 0 || k > order )
          return false;
      }
      return true;
    }

    //! \brief return begin iterator
    static IteratorType begin () { return IteratorType::begin(); }
    //! \brief return end iterator
    static IteratorType end () { return IteratorType::end(); }
  };



  // Implementation of MultiIndexSet::Iterator
  // -----------------------------------------

  template< int dimension, int order >
  class MultiIndexSet< dimension, order >::Iterator
  {
    typedef Iterator ThisType;

  protected:
    // index type
    typedef int IndexType;

    // constructor
    explicit Iterator ( IndexType n )
    : multiIndex_( n )
    {}

    static IndexType invalidIndex ()
    {
      return std::numeric_limits< IndexType >::max();
    }

  public:
    //! \brief return begin iterator
    static ThisType begin () { return ThisType( 0 ); }
    //! \brief return end iterator
    static ThisType end () { return ThisType( invalidIndex() ); }

    //! \brief iterator increment
    ThisType operator++ ()
    {
      // try to increment and leave eventually
      for( int i = 0; i < dimension; ++i )
      {
        const int j = dimension-i-1;
        if( ++multiIndex_[ j ] <= order )
          return *this;
        multiIndex_[ j ] = 0;
      }
      
      // otherwise, reset this iterator to end iterator
      *this = end(); 
      return *this;
    }

    //! \brief check for equality
    bool operator== ( const ThisType &other ) const
    {
      return ( multiIndex_ == other.multiIndex_ );
    }

    //! \brief check for inequality
    bool operator!= ( const ThisType &other ) const
    {
      return !( *this == other );
    }

    //! \brief iterator dereference
    const MultiIndexType &operator* () const
    {
      return multiIndex_;
    }

    //! \brief iterator dereference
    const MultiIndexType *operator-> () const
    {
      return &multiIndex_;
    }

  private:
    MultiIndexType multiIndex_;
  };

} // namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_MULTIINDEXSET_HH
