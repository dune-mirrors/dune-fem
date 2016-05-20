#ifndef DUNE_FEM_SIMPLEBLOCKVECTOR_HH
#define DUNE_FEM_SIMPLEBLOCKVECTOR_HH

#include <algorithm>
#include <cassert>
#include <memory>

#include <dune/common/densevector.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/misc/debug.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/storage/subvector.hh>
#include <dune/fem/space/common/arrays.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#endif

namespace Dune {

  namespace Fem {

  // tag for block vectors
  struct IsBlockVector {};



  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  Imp     Implementation
  *   \tparam  Field   Field type of all dofs
  */
  template< class Imp, class Field >
  class BlockVectorInterface
  : public IsBlockVector
  {
  protected:
    typedef DebugCounter<size_t> CounterType;

    BlockVectorInterface () {}

    //! Type of derived class (implementation)
    typedef Imp    ThisType;
  public:
    //! Type of the field the dofs lie in
    typedef Field  FieldType;

    /** \brief Copy assignment operator */
    const ThisType& operator= ( const ThisType &other )
    {
      if( &asImp() != &other )
      {
        assign( other );
        sequence_ = other.sequence_;
      }
      return asImp();
    }

    /** \brief Add another block vector to *this */
    const ThisType &operator+= ( const ThisType &other )
    {
      assert( asImp().size() == other.size() );
      const auto endit = asImp().end();
      auto oit = other.begin();
      for( auto it = asImp().begin(); it != endit; ++it, ++oit )
        *it += *oit;

      ++sequence_;
      return asImp();
    }

    /** \brief Subtract another block vector from *this */
    const ThisType &operator-= ( const ThisType &other )
    {
      assert( asImp().size() == other.size() );
      const auto endit = asImp().end();
      auto oit = other.begin();
      for( auto it = asImp().begin(); it != endit; ++it, ++oit )
        *it -= *oit;

      ++sequence_;
      return asImp();
    }

    /** \brief Scalar product between *this and another block vector */
    FieldType operator* ( const ThisType &other ) const
    {
      assert( asImp().size() == other.size() );
      FieldType sum( 0 );
      const auto endit = asImp().end();
      auto oit = other.asImp().begin();
      for( auto it = asImp().begin(); it != endit; ++it, ++oit )
        sum += (*it * *oit);

      return sum;
    }

    /** \brief  Scale this block vector
     *
     *  \param[in] scalar   Factor for the scaling
     *  \return   Constant reference to *this
     */
    const ThisType &operator*= ( const FieldType &scalar )
    {
      const auto endit = asImp().end();
      for( auto it = asImp().begin(); it != endit; ++it )
        *it *= scalar;

      ++sequence_;
      return asImp();
    }

    /** \brief Add a scalar multiple of another block vector to this block vector.
     *
     *    Semantic in pseudocode: " *this = *this + scalar*v "
     *
     *  \param[in] scalar  Scalar factor by which v is multiplied before it is added to *this
     *  \param[in] other   The other block vector
     */
    void axpy ( const FieldType &scalar, const ThisType &other )
    {
      assert( asImp().size() == other.size() );
      const auto endit = asImp().end();
      auto oit = other.begin();
      for( auto it = asImp().begin(); it != endit; ++it, ++oit )
        *it += scalar * (*oit);

      ++sequence_;
    }

    /** \brief  Clear this block vector, i.e. set each dof to 0 */
    void clear ()
    {
      std::fill( asImp().begin(), asImp().end(), FieldType( 0 ) );
      ++sequence_;
    }

  protected:
    // Copy block vectors.
    //    Note: No '++sequence_' here, sequence_ is only changed in public methods
    void assign ( const ThisType &other )
    {
      assert( asImp().size() == other.size() );
      std::copy( other.begin(), other.end(), asImp().begin() );
    }

    ThisType& asImp() { return static_cast< ThisType& > (*this); }
    const ThisType& asImp() const { return static_cast< ThisType& > (*this); }

    mutable CounterType sequence_; // for consistency checks...
  };



  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  Container   Container type
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Container, int BlockSize >
  class SimpleBlockVector
  : public BlockVectorInterface< SimpleBlockVector< Container, BlockSize>, typename Container::value_type >
  {
    typedef BlockVectorInterface< SimpleBlockVector< Container, BlockSize>, typename Container::value_type > BaseType;
    typedef SimpleBlockVector< Container, BlockSize > ThisType;
    typedef Container                                 ArrayType;

    typedef DebugCounter<size_t>                      CounterType;

  public:
    typedef ArrayType DofContainerType;

    //! Type of the field the dofs lie in
    typedef typename ArrayType::value_type                  FieldType;
    //! Iterator to iterate over the dofs
    typedef typename ArrayType::iterator                    IteratorType;
    //! Constant iterator to iterate over the dofs
    typedef typename ArrayType::const_iterator              ConstIteratorType;
    //! Used for indexing the blocks, for example
    typedef typename ArrayType::size_type                   SizeType;

    //! Typedef to make this class STL-compatible
    typedef FieldType           value_type;
    //! Typedef to make this class STL-compatible
    typedef SizeType            size_type;

    //! Type of one (mutable) block
    typedef SubVector< DofContainerType, StaticOffsetSubMapper< BlockSize > > DofBlockType;
    //! Type of one constant block
    typedef DofBlockType ConstDofBlockType;

    typedef Fem::Envelope< DofBlockType >       DofBlockPtrType;
    typedef Fem::Envelope< ConstDofBlockType >  ConstDofBlockPtrType;

    //! Size of each block
    enum { blockSize = BlockSize };

    /** \brief Constructor */
    explicit SimpleBlockVector ( ArrayType& array )
    : array_( array )
    {}

    /** \brief Copy assignment operator */
    const ThisType& operator= ( const ThisType& other )
    {
      BaseType::operator=( other );
      return *this;
    }

    /** \brief Constant access the i-th block */
    ConstDofBlockType operator[] ( const unsigned int i ) const
    {
      assert( i < size() );
      return ConstDofBlockType( array_, StaticOffsetSubMapper< blockSize >( i*blockSize ) );
    }

    /** \brief Access the i-th block */
    DofBlockType operator[] ( const unsigned int i )
    {
      assert( i < size() );
      return DofBlockType( array_, StaticOffsetSubMapper< blockSize >( i*blockSize ) );
    }

    /** \brief Constant access for the i-th block */
    ConstDofBlockPtrType blockPtr( const unsigned int i ) const
    {
      return ConstDofBlockPtrType( this->operator[] ( i ) );
    }

    /** \brief Access the i-th block */
    DofBlockPtrType blockPtr( const unsigned int i )
    {
      return DofBlockPtrType( this->operator[] ( i ) );
    }

    /** \brief Iterator pointing to the first dof */
    IteratorType begin() { return array().begin(); }

    /** \brief Const-iterator pointing to the first dof */
    ConstIteratorType begin() const { return array().begin(); }

    /** \brief Iterator pointing to the last dof */
    IteratorType end() { return array().end(); }

    /** \brief Const-iterator pointing to the last dof */
    ConstIteratorType end() const { return array().end(); }

    /** \brief Number of blocks */
    SizeType size () const { return array().size() / blockSize; }

    /** \brief Number of dofs in the block vector */
    SizeType numDofs() const { return array().size(); }

    FieldType* data() { return array().data(); }
    const FieldType* data() const { return array().data(); }

    const ArrayType &array () const { return array_; }
    ArrayType &array () { return array_; }

  protected:
    ArrayType& array_;
  };



  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  Container   Container type
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Container, unsigned int BlockSize >
  class MutableBlockVector
  : public SimpleBlockVector< Container, BlockSize >
  {
    typedef MutableBlockVector< Container, BlockSize > ThisType;
    typedef SimpleBlockVector < Container, BlockSize > BaseType;

    typedef Container                                  ArrayType;
    using BaseType :: array_;
    using BaseType :: sequence_;
  public:

    using BaseType :: array;
    using BaseType :: blockSize ;
    typedef typename BaseType :: SizeType SizeType;

    /** \brief Construct a block vector with 'size' blocks (not initialized) */
    explicit MutableBlockVector ( SizeType size )
    : BaseType( *(new Container( size*blockSize ) ) )
    {}

    /** \brief Copy constructor */
    MutableBlockVector ( const ThisType &other )
    : BaseType( *(new Container( other.array().size() ) ) )
    {
      assign( other );
    }

    ~MutableBlockVector()
    {
      delete &array_;
    }

    /** \brief Reserve memory.
     *
     *  This method is a no-op. It is defined here to make the block vector
     *  compatible to the managed dof storage mechanisms used by
     *  Dune::Fem::BlockVectorDiscreteFunction
     *
     *  \param[in] size  Number of blocks
     */
    void reserve ( const int size )
    {
      array().reserve( size*blockSize );
    }

    /** \brief Resize the block vector */
    void resize ( SizeType size )
    {
      array().resize( size*blockSize );
      ++sequence_;
    }
  };



  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  Field       Field type of all dofs
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Field, unsigned int BlockSize >
  class MutableBlockVector< MutableArray< Field >, BlockSize >
  : public SimpleBlockVector< StaticArray< Field >, BlockSize >
  {
    typedef StaticArray< Field >          StaticContainer ;
    typedef MutableArray< Field >         MutableContainer ;
    typedef SimpleBlockVector< StaticContainer, BlockSize >   BaseType;
    typedef MutableBlockVector< MutableContainer, BlockSize > ThisType;

  protected:
    using BaseType :: array_;
    using BaseType :: sequence_;

    std::unique_ptr< MutableContainer > container_;
  public:
    using BaseType :: blockSize ;
    typedef typename BaseType :: SizeType SizeType;

    /** \brief Construct a block vector with 'size' blocks (not initialized) */
    explicit MutableBlockVector ( SizeType size )
    : BaseType( allocateContainer( size*blockSize ) ),
      container_( static_cast< MutableContainer* > (&array_) )
    {}

    /** \brief Copy constructor */
    MutableBlockVector ( const ThisType &other )
    : BaseType( allocateContainer( other.array().size() ) ),
      container_( static_cast< MutableContainer* > (&array_) )
    {
      assign( other );
    }

    /** \brief Reserve memory.
     *
     *  This method is a no-op. It is defined here to make the block vector
     *  compatible to the managed dof storage mechanisms used by
     *  Dune::Fem::BlockVectorDiscreteFunction
     *
     *  \param[in] size  Number of blocks
     */
    void reserve ( const int size )
    {
      assert( container_ );
      container_->reserve( size*blockSize );
    }

    /** \brief Resize the block vector */
    void resize ( SizeType size )
    {
      assert( container_ );
      container_->resize( size*blockSize );
      ++sequence_;
    }

  protected:
    StaticContainer& allocateContainer( const SizeType size )
    {
      MutableContainer* container = new MutableContainer( size );
      return *container;
    }
  };



  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  */
  template< class DofBlock >
  class ISTLBlockVector
  : public BlockVectorInterface< ISTLBlockVector< DofBlock >, typename DofBlock :: value_type >
  {
    ISTLBlockVector ( const ISTLBlockVector& );
    typedef ISTLBlockVector< DofBlock>                            ThisType;
#if HAVE_DUNE_ISTL
    typedef BlockVector< DofBlock >                               ArrayType;
#else
    // fallback in case dune-istl is not present
    typedef Dune::DynamicVector< DofBlock >                       ArrayType;
#endif
    typedef BlockVectorInterface< ISTLBlockVector< DofBlock >, typename DofBlock :: value_type  >  BaseType;


    using BaseType :: sequence_;
  public:
    typedef ArrayType DofContainerType;

    enum { blockSize = DofBlock :: dimension };

    typedef typename DofBlock :: value_type FieldType;
  protected:
    template <class EmbeddedIterator, class V>
    class Iterator
      : public ForwardIteratorFacade< Iterator< EmbeddedIterator,V >, V >
    {
    public:
      typedef V FieldType;
    protected:
      mutable EmbeddedIterator it_;
#ifndef NDEBUG
      const   EmbeddedIterator end_;
#endif
      int index_;
    public:
      //! Default constructor
      Iterator( const EmbeddedIterator& it
#ifndef NDEBUG
              , const EmbeddedIterator& end = EmbeddedIterator()
#endif
              )
        : it_( it ),
#ifndef NDEBUG
          end_( end ),
#endif
          index_(0)
      {}

      //! return dof
      FieldType& dereference () const
      {
        assert( it_ != end_ );
        assert( index_ < blockSize );
        return (*it_)[ index_ ];
      }

      //! go to next dof
      void increment ()
      {
        ++index_;
        if( index_ >= blockSize )
        {
          index_ = 0;
          ++it_;
        }
      }

      //! compare
      bool equals ( const Iterator &other ) const
      {
        return (it_ == other.it_) && (index_ == other.index_);
      }

    }; // end DofIteratorBlockVectorDiscreteFunction

  public:
    typedef Iterator< typename ArrayType::Iterator, FieldType >            IteratorType;
    typedef Iterator< typename ArrayType::ConstIterator, const FieldType > ConstIteratorType;

    typedef DofBlock DofBlockType;
    typedef const DofBlock ConstDofBlockType;

    typedef DofBlockType*  DofBlockPtrType;
    typedef ConstDofBlockType*  ConstDofBlockPtrType;

    typedef typename ArrayType::size_type   SizeType;
    //! Typedef to make this class STL-compatible
    typedef typename ArrayType::value_type  value_type;

    /** \brief Constructor */
    explicit ISTLBlockVector ( ArrayType& array )
    : array_( &array )
    {}

    ISTLBlockVector () : array_( 0 )
    {}

    /** \brief Copy assignment operator */
    const ThisType& operator= ( const ThisType& other )
    {
      if( this != &other )
      {
        array() = other.array();
      }
      return *this;
    }

    DofBlockPtrType blockPtr(const unsigned int i ) { return &array()[ i ]; }
    ConstDofBlockPtrType blockPtr(const unsigned int i ) const { return &array()[ i ]; }

    DofBlockType& operator[] (const unsigned int i ) { return array()[ i ]; }
    ConstDofBlockType& operator[] (const unsigned int i ) const { return array()[ i ]; }

    IteratorType begin() { return IteratorType( array().begin()
#ifndef NDEBUG
                                              , array().end()
#endif
                                              ); }
    ConstIteratorType begin() const
    {
      return ConstIteratorType( array().begin()
#ifndef NDEBUG
                              , array().end()
#endif
                              ); }

    IteratorType end() { return IteratorType( array().end() ); }
    ConstIteratorType end() const  { return ConstIteratorType( array().end() ); }

    SizeType size() const { return array().size(); }

    /** \brief Reserve memory.
     *
     *  This method is a no-op. It is defined here to make the block vector
     *  compatible to the managed dof storage mechanisms used by
     *  Dune::Fem::BlockVectorDiscreteFunction
     *
     *  \param[in] size  Number of blocks
     */
    void reserve ( const int size )
    {
      array().reserve( size );
    }

    /** \brief Resize the block vector */
    void resize ( SizeType size )
    {
      array().resize( size );
      ++sequence_;
    }

    ArrayType& array() { assert( array_ ); return *array_; }
    const ArrayType& array() const { assert( array_ ); return *array_; }

  protected:
    // ISTL BlockVector
    ArrayType* array_;
  };

} // namespace Fem
} // namespace Dune

#endif // DUNE_FEM_REFERENCEBLOCKVECTOR_HH
