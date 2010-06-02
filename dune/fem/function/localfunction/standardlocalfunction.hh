#ifndef DUNE_FEM_STANDARDLOCALFUNCTION_HH
#define DUNE_FEM_STANDARDLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{
  
  // StandardLocalFunctionImpl
  // -------------------------

  /** \class StandardLocalFunctionImpl
   *  \brief standard implementation of a local function
   */
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  class StandardLocalFunctionImpl
  : public LocalFunctionDefault
    < DiscreteFunctionSpace, StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > >
  {
    typedef StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpace, ThisType > BaseType;

  public:
    //! type of discrete function the local function belongs to
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of  discrete function space the local function belongs to
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    //! type of grid
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    
    //! type of underlying function space
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
   
    //! field type for domain vectors
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    //! field type for range vectors
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    //! type of domain vectors
    typedef typename FunctionSpaceType::DomainType DomainType;
    //! type of range vectors
    typedef typename FunctionSpaceType::RangeType RangeType;
    //! type of the Jacobian
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
      BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType::template Codim< 0 >::Entity EntityType;

    //! constructor
    explicit StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction );
    
    //! copy constructor
    StandardLocalFunctionImpl ( const ThisType &other );

    /** \copydoc Dune::LocalFunction::operator[](const int num) const */
    const RangeFieldType &operator[] ( const int num ) const;

    /** \copydoc Dune::LocalFunction::operator[](const int num) */
    RangeFieldType &operator[] ( const int num );

    /** \copydoc Dune::LocalFunction::order() const */
    int order () const;

    /** \copydoc Dune::LocalFunction::baseFunctionSet() const */
    const BaseFunctionSetType &baseFunctionSet() const;

    /** \copydoc Dune::LocalFunction::entity() const */
    const EntityType &entity () const;

    //! initialize local function 
    void init ( const EntityType &entity );
    
    /** \copydoc Dune::LocalFunction::numDofs() const */
    int numDofs () const;

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  protected:
    DiscreteFunctionType &discreteFunction_;
    
    // array holding pointer to local dofs 
    DynamicArray< RangeFieldType* > values_;

     // base function set 
    BaseFunctionSetType baseFunctionSet_;

    // actual entity
    const EntityType *entity_;

    // number of local dofs
    unsigned int numDofs_;

    bool needCheckGeometry_;
  };



  // StandardLocalFunctionFactory
  // ----------------------------

  template< class DiscreteFunctionTraits >
  class StandardLocalFunctionFactory
  {
    typedef StandardLocalFunctionFactory< DiscreteFunctionTraits > ThisType;

  public:
    typedef typename DiscreteFunctionTraits::DiscreteFunctionType DiscreteFunctionType;
    typedef typename DiscreteFunctionTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef StandardLocalFunctionImpl< DiscreteFunctionType, DiscreteFunctionSpaceType > ObjectType;

    explicit StandardLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {}

    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_ );
    }

  protected:
    DiscreteFunctionType &discreteFunction_;
  };



  // Implementation of StandardLocalFunctionImpl
  // -------------------------------------------

  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
  : discreteFunction_( discreteFunction ),
    values_( discreteFunction_.space().maxNumLocalDofs() ),
    baseFunctionSet_(),
    entity_( 0 ),
    numDofs_( 0 ),
    needCheckGeometry_( true )
  {}


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::StandardLocalFunctionImpl ( const ThisType &other )
  : discreteFunction_( other.discreteFunction_ ),
    values_( other.values_ ),
    baseFunctionSet_( other.baseFunctionSet_ ),
    entity_( other.entity_ ),
    numDofs_( other.numDofs_ ),
    needCheckGeometry_( other.needCheckGeometry_ )
  {}


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::RangeFieldType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::operator[] ( const int num ) const
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }
  
  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::RangeFieldType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::operator[] ( const int num )
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::order () const
  {
    return discreteFunction_.space().order();
  }


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::BaseFunctionSetType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: baseFunctionSet () const
  {
    assert( entity_ != 0 );
    return baseFunctionSet_;
  }

  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::EntityType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::entity () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline void
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::init ( const EntityType &entity )
  {
    typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
    typedef typename BlockMapperType :: DofMapIteratorType DofMapIteratorType;
    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

    typedef typename DiscreteFunctionType :: DofBlockPtrType DofBlockPtrType;
    
    const DiscreteFunctionSpaceType &space = discreteFunction_.space();
    const bool multipleBaseSets = space.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      // if multiple base sets skip geometry call
      bool updateBaseSet = true;
      if( !multipleBaseSets && (entity_ != 0) )
        updateBaseSet = (baseFunctionSet_.geometryType() != entity.type());
      
      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = space.baseFunctionSet( entity );

        // note, do not use baseFunctionSet() here, entity might no have been set
        numDofs_ = baseFunctionSet_.numBaseFunctions();

        needCheckGeometry_ = space.multipleGeometryTypes();
      }
    }

    // cache entity
    entity_ = &entity;
    assert( baseFunctionSet_.geometryType() == entity.type() );

    assert( numDofs_ <= values_.size() );
    const BlockMapperType &mapper = space.blockMapper();
    const DofMapIteratorType end = mapper.end( entity );
    for( DofMapIteratorType it = mapper.begin( entity ); it != end; ++it )
    {
      assert( it.global() == mapper.mapToGlobal( entity, it.local() ) );
      
      DofBlockPtrType blockPtr = discreteFunction_.block( it.global() );
      
      const unsigned int localBlock = it.local() * blockSize;
      for( unsigned int i = 0; i < blockSize; ++i )
        values_[ localBlock + i ] = &((*blockPtr)[ i ]);
    }
  }


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::numDofs () const
  {
    return numDofs_;
  }

}

#endif // #ifndef DUNE_FEM_STANDARDLOCALFUNCTION_HH
