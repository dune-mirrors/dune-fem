#ifndef DUNE_FEM_SIMPLEBLOCKVECTOR_HH
#define DUNE_FEM_SIMPLEBLOCKVECTOR_HH

#include <algorithm>
#include <cassert>

#include <dune/common/densevector.hh>

#include <dune/fem/misc/debug.hh> // for DebugCounter
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/space/common/arrays.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#endif

namespace Dune {
  namespace Fem {
    // Forward declaration
    template< class BlockVector >
    class SimpleBlockVectorBlock;
  }

  template< class BlockVector >
  struct DenseMatVecTraits< Fem::SimpleBlockVectorBlock< BlockVector > >
  {
    typedef Fem::SimpleBlockVectorBlock< BlockVector > derived_type;
    typedef Fem::SimpleBlockVectorBlock< BlockVector > container_type;
    typedef typename BlockVector :: FieldType     value_type;
    typedef typename BlockVector :: SizeType      size_type;
  };

namespace Fem {

  // tag for block vectors
  struct IsBlockVector {};

  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Imp, class Field >
  class BlockVectorInterface
  : public IsBlockVector
  {
  protected:
    typedef DebugCounter<size_t>                      CounterType;

    BlockVectorInterface () {}

    //! Type of derived class (implementation)
    typedef Imp    ThisType;
  public:
    //! Type of the field the dofs lie in
    typedef Field  FieldType;

    /*
     * ########## operators ##############################
     */
    /** \brief Copy assignment operator
     */
    const ThisType& operator= ( const ThisType &other )
    {
      if( &asImp() != &other )
      {
        assign( other );
        sequence_ = other.sequence_;
      }
      return asImp();
    }

    /** \brief Add another block vector to *this
     *
     *  \param[in] other    Other block vector to add
     *  \return Constant reference to *this
     */
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

    /** \brief Subtract another block vector from *this
     *
     *  \param[in] other    Other block vector to subtract
     *  \return Constant reference to *this
     */
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

    /** \brief Scalar product *this with another block vector
     *
     *  \param[in] other  Other block vector
     *  \return Returns the scalar product " (*this)*other"
     */
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

    /*
     * ########## methods ##############################
     */

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
      // TODO: revise
      for( auto it = asImp().begin(); it != endit; ++it, ++oit )
      {
        *it += scalar * (*oit);
      }

      ++sequence_;
    }

    /** \brief Clear this block vector, i.e. set each dof to 0
     */
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
  *   \tparam  F           The ground fields. All dofs are elements of this field.
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

    friend class SimpleBlockVectorBlock< ThisType >;
    friend class SimpleBlockVectorBlock< const ThisType >;

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
    typedef SimpleBlockVectorBlock< ThisType >        DofBlockType;
    //! Type of one constant block
    typedef SimpleBlockVectorBlock< ThisType >  ConstDofBlockType;

    typedef Fem::Envelope< DofBlockType >       DofBlockPtrType;
    typedef Fem::Envelope< ConstDofBlockType >  ConstDofBlockPtrType;

    //! Size of each block
    enum { blockSize = BlockSize };


    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit SimpleBlockVector ( ArrayType& array )
    : array_( array )
    {}

    /*
     * ########## operators ##############################
     */
    /** \brief Copy assignment operator
     */
    const ThisType& operator= ( const ThisType& other )
    {
      BaseType::operator=( other );
      return *this;
    }

    /** \brief Accessor for a constant block
     *
     *  \param[in]  i         Index of the block
     *  \return   The i-th block, constant
     */
    ConstDofBlockType operator[] ( const unsigned int i ) const
    {
      assert( i < size() );
      return ConstDofBlockType( *this, i*blockSize );
    }

    /** \brief Accessor for a block
     *
     *  \param[in]  i         Index of the block
     *  \return   The i-th block
     */
    DofBlockType operator[] ( const unsigned int i )
    {
      assert( i < size() );
      return DofBlockType( *this, i*blockSize );
    }

    /** \brief Accessor for a constant block
     *
     *  \param[in]  i         Index of the block
     *  \return   The i-th block, constant
     */
    ConstDofBlockPtrType blockPtr( const unsigned int i ) const
    {
      return ConstDofBlockPtrType( this->operator[] ( i ) );
    }

    /** \brief Accessor for a block
     *
     *  \param[in]  i         Index of the block
     *  \return   The i-th block
     */
    DofBlockPtrType blockPtr( const unsigned int i )
    {
      return DofBlockPtrType( this->operator[] ( i ) );
    }

    /** \brief Iterator pointing to the first dof
     *
     *  \return Iterator pointing to the first dof
     */
    IteratorType begin() { return array().begin(); }

    /** \brief Const-iterator pointing to the first dof
     *
     *  \return Const-iterator pointing to the first dof
     */
    ConstIteratorType begin() const { return array().begin(); }

    /** \brief Iterator pointing to the last dof
     *
     *  \return Iterator pointing to the last dof
     */
    IteratorType end() { return array().end(); }

    /** \brief Const-iterator pointing to the last dof
     *
     *  \return Const-iterator pointing to the last dof
     */
    ConstIteratorType end() const { return array().end(); }

    /** \brief Returns the number of blocks
     *
     *  \return Number of blocks
     */
    SizeType size () const { return array().size() / blockSize; }

    /** \brief Returns the number of dofs in the block vector
     *
     *  \return Number of dofs
     */
    SizeType numDofs() const { return array().size(); }

    FieldType* data() { return array().data(); }
    const FieldType* data() const { return array().data(); }

    const ArrayType &array () const { return array_; }
    ArrayType &array () { return array_; }

  protected:
    /*
     * data fields
     */
    ArrayType& array_;
  };

  /** \class SimpleBlockVectorBlock
  *   \brief This is the implementation of a block of SimpleBlockVector
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class BlockVector >
  class SimpleBlockVectorBlock : public Dune::DenseVector< SimpleBlockVectorBlock< BlockVector > >
  {
    typedef BlockVector  BlockVectorType;
    typedef typename BlockVectorType :: FieldType     FieldType;
    typedef typename BlockVectorType::CounterType     CounterType;

    typedef SimpleBlockVectorBlock< BlockVectorType > ConstBlockType;
    typedef SimpleBlockVectorBlock< BlockVectorType > ThisType;

  public:
    typedef typename BlockVectorType :: SizeType size_type;

    //! The block size
    static const unsigned int blockSize = BlockVector :: blockSize ;

    /** \brief Standard constructor for SimpleBlockVectorBlocks
     *
     *  \param[in]  blockVector   The block vector in which this block lives
     *  \param[in]  blockBegin    Beginning index of this block in the block vector's array (implementation detail)
     */
    SimpleBlockVectorBlock ( const BlockVectorType &blockVector, unsigned int blockBegin )
    : blockVector_( const_cast< BlockVectorType& > (blockVector) ),
      blockBegin_( blockBegin )
    {}

    /** \brief Copy constructor
     */
    SimpleBlockVectorBlock ( const SimpleBlockVectorBlock< BlockVectorType > &other )
    : blockVector_( const_cast< BlockVectorType& > (other.blockVector_)),
      blockBegin_( other.blockBegin_ )
    {}

    size_type size () const { return blockSize; }
    size_type vec_size () const { return blockSize; }

    /** \brief Copy assignment operator for constant blocks
     *
     *  \param[in] other  Other block (constant) which should be assigned to *this
     *  \return  Constant reference to *this
     */
    ThisType &operator= ( const ConstBlockType &other )
    {
      copy( other );
      return *this;
    }
#if 0

    /** \brief Add another block to *this
     *
     *  \param[in] other  Other block to add
     *  \return Constant reference to *this
     */
    ThisType& operator+= ( const ConstBlockType& other )
    {
      assert( sequence_ == blockVector_.sequence_ );
      for ( unsigned int i = 0; i < blockSize; ++i )
      {
        (*this)[i] += other[i];
      }
      ++sequence_;
      return *this;
    }

    /** \brief Subtract another block from *this
     *
     *  \param[in] other  Other block to subtract
     *  \return Constant reference to *this
     */
    ThisType& operator-= ( const ConstBlockType& other )
    {
      assert( sequence_ == blockVector_.sequence_ );
      for ( unsigned int i = 0; i < blockSize; ++i )
      {
        (*this)[i] -= other[i];
      }
      ++sequence_;
      return *this;
    }

    /** \brief Calculate the scalar product of this block with another block
     *
     *  \param[in] other  Other block to scalar-multiply this block by
     *  \return The value of the scalar product
     */
    FieldType operator* ( const ConstBlockType& other ) const
    {
      assert( sequence_ == blockVector_.sequence_ );
      FieldType sum( 0 );
      for ( unsigned int i = 0; i < blockSize; ++i )
      {
        sum += (*this)[i] * other[i];
      }
      return sum;
    }

    /** \brief Scale this block
     *
     *  \param[in] scalar   Scalar to use for the scaling
     *  \return Constant reference to *this
     */
    ThisType& operator*= ( const FieldType& scalar )
    {
      assert( sequence_ == blockVector_.sequence_ );
      for ( unsigned int i = 0; i < blockSize; ++i )
      {
        (*this)[i] *= scalar;
      }
      ++sequence_;
      return *this;
    }
#endif

    /** \brief Obtain a dof inside this block
     *
     *  \param[in] index   Index of the dof
     *  \return Reference to the dof
     */
    FieldType& vec_access(unsigned int index)
    {
      assert(index < blockSize);
      return blockVector_.array()[blockBegin_ + index];
    }

    /** \brief Obtain a dof inside this block
     *
     *  \param[in] index   Index of the dof
     *  \return Constant reference to the dof
     */
    const FieldType& vec_access(unsigned int index) const
    {
      assert(index < blockSize);
      return blockVector_.array()[blockBegin_ + index];
    }

    /** \brief Returns the size of the block
     */
    int dim() const { return blockSize; }

  private:

    // An empty constructor does not make sense in this case
    SimpleBlockVectorBlock();

    template< class Block >
    void copy ( const Block &other )
    {
      assert( &blockVector_ == &other.blockVector_ );
      for( unsigned int i=0; i < blockSize; ++i )
        (*this)[ i ] = other[ i ];
    }

    // data fields
    BlockVectorType &blockVector_;
    const unsigned int blockBegin_;
  };

  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Container, unsigned int BlockSize >
  class MutableBlockVector
  : public SimpleBlockVector< Container, BlockSize >
  {
    typedef MutableBlockVector< Container, BlockSize > ThisType;
    typedef SimpleBlockVector< Container, BlockSize >  BaseType;
    typedef Container                                  ArrayType;
    using BaseType :: array_;
    using BaseType :: sequence_;
  public:

    using BaseType :: array;
    using BaseType :: blockSize ;
    typedef typename BaseType :: SizeType SizeType;

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit MutableBlockVector ( SizeType size )
    : BaseType( *(new Container( size*blockSize ) ) )
    {}

    /** \brief Copy constructor
     */
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

    /** \brief Resize the block vector
     *
     *  \param[in] size  New number of blocks
     */
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
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class F, unsigned int BlockSize >
  class MutableBlockVector< MutableArray< F >, BlockSize >
  : public SimpleBlockVector< StaticArray< F >, BlockSize >
  {
    typedef StaticArray< F >          StaticContainer ;
    typedef MutableArray< F >         MutableContainer ;
    typedef SimpleBlockVector< StaticContainer, BlockSize >   BaseType;
    typedef MutableBlockVector< MutableContainer, BlockSize > ThisType;

  protected:
    using BaseType :: sequence_;

    MutableContainer* container_;
  public:
    using BaseType :: blockSize ;
    typedef typename BaseType :: SizeType SizeType;

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit MutableBlockVector ( SizeType size )
    : BaseType( allocateContainer( size*blockSize ) )
    {}

    /** \brief Copy constructor
     */
    MutableBlockVector ( const ThisType &other )
    : BaseType( allocateContainer( other.array().size() ) )
    {
      assign( other );
    }

    ~MutableBlockVector()
    {
      delete container_;
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
      container_->reserve( size*blockSize );
    }

    /** \brief Resize the block vector
     *
     *  \param[in] size  New number of blocks
     */
    void resize ( SizeType size )
    {
      container_->resize( size*blockSize );
      ++sequence_;
    }

  protected:
    StaticContainer& allocateContainer( const SizeType size )
    {
      container_ = new MutableContainer( size );
      return *container_;
    }
  };

#if HAVE_DUNE_ISTL
  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class DofBlock >
  class ISTLBlockVector
  : public BlockVectorInterface< ISTLBlockVector< DofBlock >, typename DofBlock :: value_type >
  {
    ISTLBlockVector ( const ISTLBlockVector& );
    typedef ISTLBlockVector< DofBlock>                            ThisType;
    typedef BlockVector< DofBlock >                               ArrayType;
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

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit ISTLBlockVector ( ArrayType& array )
    : array_( &array )
    {}

    ISTLBlockVector () : array_( 0 )
    {}

    /*
     * ########## operators ##############################
     */
    /** \brief Copy assignment operator
     */
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

    /** \brief Resize the block vector
     *
     *  \param[in] size  New number of blocks
     */
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
#endif // #if HAVE_DUNE_ISTL

} // namespace Fem
} // namespace Dune

#endif // DUNE_FEM_REFERENCEBLOCKVECTOR_HH
