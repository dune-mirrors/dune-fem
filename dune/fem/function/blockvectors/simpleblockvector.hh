#ifndef DUNE_FEM_SIMPLEBLOCKVECTOR_HH
#define DUNE_FEM_SIMPLEBLOCKVECTOR_HH

#include <algorithm>
#include <cassert>

#include <dune/fem/misc/debug.hh> // for DebugCounter

namespace Dune {
namespace Fem {

  // Forward declaration
  template< class BlockVector >
  class SimpleBlockVectorBlock;

  // tag for block vectors
  struct IsBlockVector {};

  /** \class SimpleBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class Container, unsigned int BlockSize >
  class SimpleBlockVector
  : public IsBlockVector
  {
    typedef SimpleBlockVector< Container, BlockSize > ThisType;
    typedef Container                                 ArrayType;

    typedef DebugCounter<size_t>                      CounterType;

    friend class SimpleBlockVectorBlock< ThisType >;
    friend class SimpleBlockVectorBlock< const ThisType >;

  public:

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
    typedef SimpleBlockVectorBlock< const ThisType >  ConstDofBlockType;

    //! Size of each block
    enum { blockSize = BlockSize };


    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit SimpleBlockVector ( SizeType size )
    : array_( size*blockSize ),
      size_( size )
    {}

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  All the dofs are set to 'initialValue'.
     *
     *  \param[in]  size          Number of blocks
     *  \param[in]  initialValue  This is the value to which each dof is set
     */
    SimpleBlockVector ( SizeType size, const FieldType& initialValue )
    : array_( size*blockSize, initialValue ),
      size_( size )
    {}

    /** \brief Copy constructor
     */
    SimpleBlockVector ( const ThisType &other )
    : array_( other.size()*blockSize ),
      size_( other.size() ),
      sequence_( other.sequence_ )
    {
      assign( other );
    }

    /*
     * ########## operators ##############################
     */
    /** \brief Copy assignment operator
     */
    const ThisType& operator= ( const ThisType &other )
    {
      if( this != &other )
      {
        assign( other );
        sequence_ = other.sequence_;
      }
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

    /** \brief Add another block vector to *this
     *
     *  \param[in] other    Other block vector to add
     *  \return Constant reference to *this
     */
    const ThisType &operator+= ( const ThisType &other )
    {
      assert( size() == other.size() );
      for( SizeType i=0; i < size(); ++i )
        (*this)[ i ] += other[ i ];

      ++sequence_;
      return *this;
    }

    /** \brief Subtract another block vector from *this
     *
     *  \param[in] other    Other block vector to subtract
     *  \return Constant reference to *this
     */
    const ThisType &operator-= ( const ThisType &other )
    {
      assert( size() == other.size() );
      for( SizeType i=0; i < size(); ++i )
        (*this)[ i ] -= other[ i ];

      ++sequence_;
      return *this;
    }

    /** \brief Scalar product *this with another block vector
     *
     *  \param[in] other  Other block vector
     *  \return Returns the scalar product " (*this)*other"
     */
    FieldType operator* ( const ThisType &other ) const
    {
      assert( size() == other.size() );
      FieldType sum( 0 );
      for( SizeType i=0; i < size(); ++i )
        sum += (*this)[ i ] * other[ i ];

      return sum;
    }

    /** \brief  Scale this block vector
     *
     *  \param[in] scalar   Factor for the scaling
     *  \return   Constant reference to *this
     */
    const ThisType &operator*= ( const FieldType &scalar )
    {
      for( SizeType i=0; i < size(); ++i )
        (*this)[ i ] *= scalar;

      ++sequence_;
      return *this;
    }

    /*
     * ########## methods ##############################
     */

    /** \brief Add a scalar multiple of another block vector to this block vector.
     *
     *    Semantic in pseudocode: " *this = *this + scalar*v "
     *
     *  \param[in] v       The other block vector
     *  \param[in] scalar  Scalar factor by which v is multiplied before it is added to *this
     */
    void addScaled ( const ThisType &v, const FieldType &scalar )
    {
      assert( size() == v.size() );
      for( SizeType i=0; i < size(); ++i )
      {
        DofBlockType thisBlock = (*this)[ i ];
        ConstDofBlockType otherBlock = v[ i ];
        for ( unsigned int bi=0; bi < blockSize; ++bi )
        {
          thisBlock[bi] += scalar*otherBlock[bi];
        }
      }
      ++sequence_;
    }

    /** \brief Clear this block vector, i.e. set each dof to 0
     */
    void clear ()
    {
      std::fill( begin(), end(), FieldType( 0 ) );
      ++sequence_;
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

    /** \brief Reserve memory.
     *
     *  This method is a no-op. It is defined here to make the block vector
     *  compatible to the managed dof storage mechanisms used by
     *  Dune::Fem::BlockVectorDiscreteFunction
     *
     *  \param[in] size  Number of blocks
     */
    void reserve ( const int size )
    {}

    /** \brief Resize the block vector
     *
     *  \param[in] size  New number of blocks
     */
    void resize ( SizeType size )
    {
      size_ = size;
      array().resize( size_*blockSize );
      ++sequence_;
    }

    /** \brief Returns the number of blocks
     *
     *  \return Number of blocks
     */
    SizeType size () const { return size_; }

  protected:

    /** \brief Returns the number of dofs in the block vector
     *
     *  \return Number of dofs
     */
    unsigned int numDofs() const { return size() * blockSize; }


  private:

    // Copy block vectors.
    //    Note: No '++sequence_' here, sequence_ is only changed in public methods
    void assign ( const ThisType &other )
    {
      assert( size() == other.size() );
      std::copy( other.begin(), other.end(), begin() );
    }

    const ArrayType &array () const { return array_; }
    ArrayType &array () { return array_; }

    /*
     * data fields
     */
    ArrayType array_;
    SizeType size_;
    mutable CounterType sequence_; // for consistency checks...
  };


  /** \class SimpleBlockVectorBlock
  *   \brief This is the implementation of a block of SimpleBlockVector
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< class BlockVector >
  class SimpleBlockVectorBlock
  {
    typedef BlockVector  BlockVectorType;
    typedef typename BlockVectorType :: FieldType      FieldType;
    typedef typename BlockVectorType::CounterType  CounterType;

    typedef typename std::remove_const< BlockVectorType >::type NonConstBlockVectorType;
    typedef SimpleBlockVectorBlock< const NonConstBlockVectorType > ConstBlockType;
    typedef SimpleBlockVectorBlock<       NonConstBlockVectorType > NonConstBlockType;

    typedef SimpleBlockVectorBlock< BlockVectorType > ThisType;

    /*
     * friends
     */
    friend class SimpleBlockVectorBlock< const NonConstBlockVectorType >;
    friend class SimpleBlockVectorBlock<       NonConstBlockVectorType >;

  public:
    //! The block size
    static const unsigned int blockSize = BlockVector :: blockSize ;

    /** \brief Standard constructor for SimpleBlockVectorBlocks
     *
     *  \param[in]  blockVector   The block vector in which this block lives
     *  \param[in]  blockBegin    Beginning index of this block in the block vector's array (implementation detail)
     */
    SimpleBlockVectorBlock ( BlockVectorType &blockVector, unsigned int blockBegin )
    : blockVector_( blockVector ),
      blockBegin_( blockBegin ),
      sequence_( blockVector_.sequence_ )
    {}

    /** \brief Copy constructor
     */
    SimpleBlockVectorBlock ( const SimpleBlockVectorBlock< BlockVectorType > &other )
    : blockVector_( other.blockVector_),
      blockBegin_( other.blockBegin_ ),
      sequence_( other.sequence_ )
    {}

    /** \brief Copy assignment operator for constant blocks
     *
     *  \param[in] other  Other block (constant) which should be assigned to *this
     *  \return  Constant reference to *this
     */
    const ThisType &operator= ( const ConstBlockType &other )
    {
      assert( sequence_ == blockVector_.sequence_ );
      copy( other );
      sequence_ = other.sequence_;
      return *this;
    }

    /** \brief Copy assignment operator non-constant blocks
     *
     *  \param[in] other  Other block (non-constant) which should be assigned to *this
     *  \return  Constant reference to *this
     */
    const ThisType &operator= ( const NonConstBlockType &other )
    {
      assert( sequence_ == blockVector_.sequence_ );
      copy( other );
      sequence_ = other.sequence_;
      return *this;
    }

    /** \brief Add another block to *this
     *
     *  \param[in] other  Other block to add
     *  \return Constant reference to *this
     */
    const ThisType& operator+= ( const ConstBlockType& other )
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
    const ThisType& operator-= ( const ConstBlockType& other )
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
    const ThisType& operator*= ( const FieldType& scalar )
    {
      assert( sequence_ == blockVector_.sequence_ );
      for ( unsigned int i = 0; i < blockSize; ++i )
      {
        (*this)[i] *= scalar;
      }
      ++sequence_;
      return *this;
    }

    /** \brief Obtain a dof inside this block
     *
     *  \param[in] index   Index of the dof
     *  \return Reference to the dof
     */
    FieldType& operator[] (unsigned int index)
    {
      assert( sequence_ == blockVector_.sequence_ );
      assert(index < blockSize);
      return blockVector_.array()[blockBegin_ + index];
    }

    /** \brief Obtain a dof inside this block
     *
     *  \param[in] index   Index of the dof
     *  \return Constant reference to the dof
     */
    const FieldType& operator[] (unsigned int index) const
    {
      assert( sequence_ == blockVector_.sequence_ );
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
    mutable CounterType sequence_;
  };

} // namespace Fem
} // namespace Dune

#endif // DUNE_FEM_REFERENCEBLOCKVECTOR_HH
