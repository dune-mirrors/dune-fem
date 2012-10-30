#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_STANDARD_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_STANDARD_HH

// dune-fem includes
#include <dune/fem/function/localfunction/default.hh>
#include <dune/fem/storage/array.hh>


namespace Dune
{

  namespace Fem
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
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType
        BasisFunctionSetType;

      //! entity type is specified by space 
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

    protected:
      typedef  Fem :: DynamicArray< RangeFieldType * > ValuesArrayType ;

      struct AssignDofs;

    public:  
      //! constructor
      explicit StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction );
      
      //! copy constructor
      StandardLocalFunctionImpl ( const ThisType &other );

      /** \copydoc Dune::Fem::LocalFunction::operator[](const int num) const */
      const RangeFieldType &operator[] ( const int num ) const;

      /** \copydoc Dune::Fem::LocalFunction::operator[](const int num) */
      RangeFieldType &operator[] ( const int num );

      /** \copydoc Dune::Fem::LocalFunction::order() const */
      int order () const;

      /** \copydoc Dune::Fem::LocalFunction::basisFunctionSet() const */
      const BasisFunctionSetType &basisFunctionSet() const;

      /** \copydoc Dune::Fem::LocalFunction::entity() const */
      const EntityType &entity () const;

      //! initialize local function 
      void init ( const EntityType &entity );
      
      /** \copydoc Dune::Fem::LocalFunction::numDofs() const */
      int numDofs () const;

    private:
      // prohibit assignment
      ThisType &operator= ( const ThisType & );

    protected:
      DiscreteFunctionType &discreteFunction_;
      
      // array holding pointer to local dofs 
      ValuesArrayType values_;

      // base function set 
      BasisFunctionSetType basisFunctionSet_;

      // actual entity
      const EntityType *entity_;

      // number of local dofs
      unsigned int numDofs_;

      bool needCheckGeometry_;
    };



    // StandardLocalFunctionImpl::AssignDofs
    // -------------------------------------

    template< class DiscreteFunction, class DiscreteFunctionSpace >
    struct StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::AssignDofs
    {
      AssignDofs ( DiscreteFunctionType &discreteFunction, ValuesArrayType &values )
      : discreteFunction_( discreteFunction ), values_( values )
      {}

      template < class GlobalKey >
      void operator () ( const std::size_t local, const GlobalKey& globalKey )
      {
        typedef typename DiscreteFunctionType::DofBlockPtrType DofBlockPtrType;
        static const unsigned int blockSize = DiscreteFunctionSpaceType::localBlockSize;
      
        DofBlockPtrType blockPtr = discreteFunction_.block( globalKey );
        for( unsigned int i = 0; i < blockSize; ++i )
          values_[ (local*blockSize) + i ] = &((*blockPtr)[ i ]);
      }

    private:
      DiscreteFunctionType &discreteFunction_;
      ValuesArrayType &values_;
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
      values_( DiscreteFunctionSpace::localBlockSize * discreteFunction_.space().blockMapper().maxNumDofs() ),
      basisFunctionSet_(),
      entity_( 0 ),
      numDofs_( 0 ),
      needCheckGeometry_( true )
    {}


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::StandardLocalFunctionImpl ( const ThisType &other )
    : discreteFunction_( other.discreteFunction_ ),
      values_( other.values_ ),
      basisFunctionSet_( other.basisFunctionSet_ ),
      entity_( other.entity_ ),
      numDofs_( other.numDofs_ ),
      needCheckGeometry_( other.needCheckGeometry_ )
    {}


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::RangeFieldType &
    StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::operator[] ( const int num ) const
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
      return discreteFunction_.space().order( entity() );
    }


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::BasisFunctionSetType &
    StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::basisFunctionSet () const
    {
      assert( entity_ != 0 );
      return basisFunctionSet_;
    }

    
    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::EntityType &
    StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::entity () const
    {
      assert( entity_ != 0 );
      return *entity_;
    }


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline void
    StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::init ( const EntityType &entity )
    {
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;

      const DiscreteFunctionSpaceType &space = discreteFunction_.space();
      const bool multipleBaseSets = space.multipleBasisFunctionSets();

      if( multipleBaseSets || needCheckGeometry_ )
      {
        // if multiple base sets skip geometry call
        bool updateBaseSet = true;
        if( !multipleBaseSets && (entity_ != 0) )
          updateBaseSet = (basisFunctionSet_.type() != entity.type());
        
        if( multipleBaseSets || updateBaseSet )
        {
          basisFunctionSet_ = space.basisFunctionSet( entity );

          // note, do not use basisFunctionSet() here, entity might no have been set
          numDofs_ = basisFunctionSet_.size();

          needCheckGeometry_ = space.multipleGeometryTypes();
        }
      }

      // cache entity
      entity_ = &entity;
      assert( basisFunctionSet_.entity().type() == entity.type() );

      assert( numDofs_ <= values_.size() );
      space.blockMapper().mapEach( entity, AssignDofs( discreteFunction_, values_ ) );
    }


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline int
    StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::numDofs () const
    {
      return numDofs_;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_STANDARD_HH
