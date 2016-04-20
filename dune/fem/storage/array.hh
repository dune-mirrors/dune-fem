#ifndef DUNE_FEM_ARRAY_HH
#define DUNE_FEM_ARRAY_HH

#include <algorithm>
#include <cassert>
#include <type_traits>

#include <dune/common/iteratorfacades.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{
  namespace Fem
  {

    /** \class ArrayInterface
     *  \ingroup VectorClasses
     *  \brief abstract array interface
     */
    template< class AT >
    class ArrayInterface
    : public BartonNackmanInterface< ArrayInterface< AT >, typename AT::ArrayType >
    {
      typedef ArrayInterface< AT > ThisType;
      typedef BartonNackmanInterface< ThisType, typename AT::ArrayType > BaseType;

    public:
      //! type of the traits
      typedef AT Traits;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::ArrayType ArrayType;

      //! type of this interface
      typedef ThisType ArrayInterfaceType;

      //! type of the array elements
      typedef typename Traits::ElementType ElementType;

      //! make consistent with std::vector
      typedef ElementType value_type ;

      //! type of constant iterator
      typedef typename Traits::ConstIteratorType ConstIteratorType;

      //! type of (non-constant) iterator
      typedef typename Traits::IteratorType IteratorType;

      //! type of constant iterator
      typedef ConstIteratorType const_iterator;

      //! type of (non-constant) iterator
      typedef IteratorType iterator;

      //! type of unsigned integral type of indexing
      typedef unsigned int  size_type;

    protected:
      using BaseType::asImp;

    public:
      /** \brief access an array element
       *
       *  \param[in]  index  index of the array element to access
       *
       *  \returns a const reference to the array element
       */
      const ElementType &operator[] ( unsigned int index ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      /** \brief access an array element
       *
       *  \param[in]  index  index of the array element to access
       *
       *  \returns a reference to the array element
       */
      ElementType &operator[] ( unsigned int index )
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      /** \brief fill the array with copies of an element
       *
       *  \param[in]  element  element wich shall be copied into every array
       *                       entry
       */
      void assign ( const ElementType &element )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( element ) );
      }

      /** \brief copy another array to this one
       *
       *  Copies the data from another array to this one. Both arrays must be of
       *  the same size.
       *
       *  \param[in]  other  array to copy
       */
      template< class T >
      void assign( const ArrayInterface< T > &other )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( other ) );
      }

      /** \brief obtain begin iterator
       *
       *  \returns an iterator pointing to the first array element
       */
      ConstIteratorType begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      /** \brief obtain begin iterator
       *
       *  \returns an iterator pointing to the first array element
       */
      IteratorType begin ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      /** \brief obtain end iterator
       *
       *  \returns an iterator pointing behind the last array element
       */
      ConstIteratorType end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      /** \brief obtain end iterator
       *
       *  \returns an iterator pointing behind the last array element
       */
      IteratorType end ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      /** obtain the size of the array
       *
       *  \returns the size of the array
       */
      unsigned int size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
        return asImp().size();
      }
    };

    template< class Element, class Array >
    class ArrayDefaultIterator
      : public ForwardIteratorFacade< ArrayDefaultIterator< Element, Array >, Element >
    {
      typedef ArrayDefaultIterator< Element, Array > ThisType;

    public:
      typedef Element ElementType;

      typedef Array ArrayType;

    protected:
      ArrayType &array_;
      unsigned int index_;

    public:
      ArrayDefaultIterator ( ArrayType &array, unsigned int index )
      : array_( array ), index_( index )
      {
        assert( index <= array.size() );
      }

      ArrayDefaultIterator( const ThisType &other )
      : array_( other.array_ ), index_( other.index_ )
      {}

      ThisType &operator= ( const ThisType &other )
      {
        assert( &(other.array_) == &array_ );
        index_ = other.index_;
      }

      ElementType& dereference() const
      {
        assert( index_ < array_.size() );
        return array_[ index_ ];
      }

      void increment()
      {
        assert( index_ < array_.size() );
        ++index_;
      }

      bool equals( const ThisType &other ) const
      {
        assert( &(other.array_) == &array_ );
        return index_ == other.index_;
      }

      unsigned int index () const
      {
        return index_;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ARRAY_HH
