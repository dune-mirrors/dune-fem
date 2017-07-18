#ifndef DUNE_FEM_COMMON_GEOMETRYTYPEMAP_HH
#define DUNE_FEM_COMMON_GEOMETRYTYPEMAP_HH

#include <array>
#include <type_traits>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune
{

  // GeometryTypeMap
  // ---------------

  /**
   * \class GeometryTypeMap
   * \brief associative container assigning values to each GeometryType
   *
   *  The GeometryTypeMap is an associative container similar to
   *  std::map< GeometryType, T >.
   *
   *  \tparam  T             value type
   *  \tparam  TypeIndexSet  type index set (see note)
   *
   *  \note T must be default constructable, copy constructable and copy
   *        assignable.
   *
   *  \note A \c TypeIndexSet is required to implement the
   *        following interface:
   *  \code
   *  struct TypeIndexSet
   *  {
   *    static constexpr std::size_t size () noexcept;
   *
   *    static std::size_t index ( const GeometryType &type ) noexcept;
   *  };
   *  \endcode
   */
  template< class T, int dim, class TypeIndex >
  class GeometryTypeMap
  {
    typedef GeometryTypeMap< T, dim, TypeIndex > This;
    typedef std::array< T, TypeIndex::size( dim ) > Container;

  public:
    //! \brief value type
    typedef typename Container::value_type value_type;
    //! \brief reference type
    typedef typename Container::reference reference;
    //! \brief const reference type
    typedef typename Container::const_reference const_reference;
    //! \brief pointer type
    typedef typename Container::pointer pointer;
    //! \brief const pointer type
    typedef typename Container::const_pointer const_pointer;
    //! \brief iterator type
    typedef typename Container::iterator iterator;
    //! \brief const iterator type
    typedef typename Container::const_iterator const_iterator;
    //! \brief size type
    typedef typename Container::size_type size_type;
    //! \brief difference type
    typedef typename Container::difference_type difference_type;

    //! \brief value type
    typedef value_type Value;
    //! \brief iterator type
    typedef iterator Iterator;
    //! \brief iterator type
    typedef const_iterator ConstIterator;
    //! \brief size type
    typedef size_type Size;

    /**
     * \name Construction and Assignment
     * \{
     */

    /** \brief default constructor */
    GeometryTypeMap () : container_() {}

    /** \brief copy constructor */
    GeometryTypeMap ( const This & ) = default;

    /** \brief move constructor */
    GeometryTypeMap ( This && ) = default;

    /** \brief copy assignment */
    This &operator= ( const This & ) = default;

    /** \brief move assignment */
    This &operator= ( This && ) = default;

    /** \} */


    /** \name Iterators
     *  \{
     */

    //! \brief return iterator to beginning
    Iterator begin () noexcept { return container_.begin(); }
    //! \brief return iterator to end
    Iterator end () noexcept { return container_.end(); }
    //! \brief return iterator to beginning
    ConstIterator begin () const noexcept { return container_.begin(); }
    //! \brief return iterator to end
    ConstIterator end () const noexcept { return container_.end(); }

    //! \brief return const_iterator to beginning
    ConstIterator cbegin () const noexcept { return container_.cbegin(); }
    //! \brief return const_iterator to end
    ConstIterator cend () const noexcept { return container_.cend(); }

    /** \} */


    /** \name Capacity
     *  \{
     */

    //! \brief return size
    constexpr Size size () const noexcept { return TypeIndex::size( dim ); }
    //! \brief return maximum size
    constexpr Size max_size () const noexcept { return size(); }
    //! \brief test whether container is empty
    constexpr bool empty () const noexcept { return (size() > 0); }

    /** \} */


    /** \name Element access
     *  \{
     */
    //! \brief access element
    Value &operator[] ( const GeometryType &type ) { return container_[ TypeIndex::index( type ) ]; }
    //! \brief access element
    const Value &operator[] ( const GeometryType &type ) const { return container_[ TypeIndex::index( type ) ]; }

    //! \brief access element
    Value &at ( const GeometryType &type ) { return container_.at( TypeIndex::index( type ) ); }
    //! \brief access element
    const Value &at ( const GeometryType &type ) const { return container_.at( TypeIndex::index( type ) ); }

    /** \brief access first element */
    Value &front () { return container_.front(); }
    /** \brief access first element */
    const Value &front () const { return container_.front(); }

    /** \brief access last element */
    Value &back () { return container_.back(); }
    /** \brief access last element */
    const Value &back () const { return container_.back(); }

    //! \brief get pointer to data
    Value *data () noexcept { return container_.data(); }
    //! \brief get pointer to data
    const Value *data () const noexcept { return container_.data(); }

    /** \}*/


    /** \name Modifiers
     *  \{
     */
    //! \brief fill container with value
    void fill ( const Value &value ) { container_.fill( value ); }

    //! \brief swap content
    void swap ( This &other ) noexcept(noexcept(swap(std::declval<Value &>(), std::declval<Value &>())))
    {
      container_.swap( other.container_ );
    }

    /** \} */

  private:
    Container container_;
  };



  // LocalGeometryTypeMap
  // --------------------

  /**
   * \class LocalGeometryTypeMap
   *
   * \brief Please doc me.
   *
   * \tparam  T        value type
   * \tparam  dim      dimension
   */
  template< class T, int dim >
  using LocalGeometryTypeMap = GeometryTypeMap< T, dim, LocalGeometryTypeIndex >;



  // GlobalGeometryTypeMap
  // ---------------------

  /**
   * \class GlobalGeometryTypeMap
   *
   * \brief Please doc me.
   *
   * \tparam  T        value type
   * \tparam  maxdim   maximum dimension
   */
  template< class T, int maxdim >
  using GlobalGeometryTypeMap = GeometryTypeMap< T, maxdim, GlobalGeometryTypeIndex >;

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_GEOMETRYTYPEMAP_HH
