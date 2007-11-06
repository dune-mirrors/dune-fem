#ifndef DUNE_STANDARDLOCALFUNCTION_HH
#define DUNE_STANDARDLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/function/common/localfunction.hh>

namespace Dune
{
  
  //- Forward declarations of Combined Space 
  template< class, int, DofStoragePolicy >
  class CombinedSpace;



  /** \class StandardLocalFunction
   *  \brief standard implementation of a local function
   */
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  class StandardLocalFunction
  : public LocalFunctionDefault
    < DiscreteFunctionSpaceImp,
      StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    >
  {
  public:
    //! type of discrete function this local function belongs to
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! type of  discrete function space this local function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef StandardLocalFunction< DiscreteFunctionType, DiscreteFunctionSpaceType >
      ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;

  public:
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    
    //! type of underlying function space
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
   
    //! field type for domain vectors
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type for range vectors
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors
    typedef typename FunctionSpaceType :: DomainType DomainType;
    //! type of range vectors
    typedef typename FunctionSpaceType :: RangeType RangeType;
    //! type of the Jacobian
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType :: DimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;

  public:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  protected:
    DiscreteFunctionType &discreteFunction_;
    
    // array holding pointer to local dofs 
    DynamicArray< RangeFieldType* > values_;

     // base function set 
    BaseFunctionSetType baseFunctionSet_;

    // actual entity
    const EntityType *entity_;

    // number of local dofs
    int numDofs_;

    bool needCheckGeometry_;

  public:
    //! constructor
    inline StandardLocalFunction ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      values_(),
      baseFunctionSet_(),
      entity_( 0 ),
      numDofs_( 0 ),
      needCheckGeometry_( true )
    {
    }
    
    //! copy constructor
    inline StandardLocalFunction ( const ThisType &other )
    : discreteFunction_( other.discreteFunction_ ),
      values_( other.values_ ),
      baseFunctionSet_( other.baseFunctionSet_ ),
      entity_( other.entity_ ),
      numDofs_( other.numDofs_ ),
      needCheckGeometry_( other.needCheckGeometry_ )
    {
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) const */
    inline const RangeFieldType &operator[] ( const int num ) const
    {
      assert( entity_ != 0 );
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }

    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) */
    inline RangeFieldType &operator[] ( const int num )
    {
      assert( entity_ != 0 );
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }

    /** \copydoc LocalFunctionInterface::numDofs */
    inline int numDofs() const
    {
      return numDofs_; 
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate( const DomainType &x,
                          RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian ( const DomainType &x,
                           JacobianRangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian( const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian( const QuadratureType &quadrature,
                          const int quadPoint,
                          JacobianRangeType &ret ) const;
   
    /** \copydoc Dune::LocalFunctionDefault::baseFunctionSet */
    const BaseFunctionSetType &baseFunctionSet() const
    {
      assert( entity_ != 0 );
      return baseFunctionSet_;
    }
    
    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const RangeType &factor) */
    inline void axpy( const DomainType &x,
                      const RangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const RangeType &factor) */
    template< class QuadratureType >
    inline void axpy( const QuadratureType &quadrature,
                      const int quadPoint,
                      const RangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const JacobianRangeType &factor) */
    inline void axpy( const DomainType &x,
                      const JacobianRangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &factor) */
    template< class QuadratureType >
    inline void axpy( const QuadratureType &quadrature,
                      const int quadPoint,
                      const JacobianRangeType &factor );
    
    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    inline void axpy( const DomainType &x,
                      const RangeType &factor1,
                      const JacobianRangeType &factor2 );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 );

    //! initialize local function 
    void init ( const EntityType &entity );

  protected:
    inline const EntityType &entity () const
    {
      assert( entity_ != 0 );
      return *entity_;
    }

    inline void rightMultiply ( const JacobianRangeType &factor,
                                const DomainType &x,
                                JacobianRangeType &result );
  };



  /** \copydoc Dune::StandardLocalFunction
   *
   *  Specialised version for CombinedSpaces
   */
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  class StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
  : public LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpaceImp, N, policy >,
      StandardLocalFunction< DiscreteFunctionImp, 
                             CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    >
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef CombinedSpace< ContainedFunctionSpaceImp, N, policy >
      DiscreteFunctionSpaceType;

  private:
    typedef StandardLocalFunction< DiscreteFunctionType, DiscreteFunctionSpaceType >
      ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;
     
  public:
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    
    //! type of underlying function space
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
   
    //! field type for domain vectors
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type for range vectors
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors
    typedef typename FunctionSpaceType :: DomainType DomainType;
    //! type of range vectors
    typedef typename FunctionSpaceType :: RangeType RangeType;
    //! type of the Jacobian
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType :: DimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;

  protected:
    typedef typename DiscreteFunctionSpaceType :: ContainedRangeType ScalarRangeType;
    typedef typename DiscreteFunctionSpaceType :: ContainedJacobianRangeType
      ScalarJacobianRangeType;
      
    template< DofStoragePolicy >
    struct DofStoragePolicyType
    {
    };

  protected:
    DiscreteFunctionType &discreteFunction_;
    
    // array holding pointer to local dofs 
    DynamicArray< RangeFieldType* > values_;

     // base function set 
    BaseFunctionSetType baseFunctionSet_;

    // actual entity
    const EntityType *entity_;

    // number of local dofs (in the scalar case)
    int numScalarDofs_;

    bool needCheckGeometry_;

  public:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! constructor
    inline StandardLocalFunction ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      values_(),
      baseFunctionSet_(),
      entity_( 0 ),
      numScalarDofs_( 0 ),
      needCheckGeometry_( true )
    {
    }
    
    //! copy constructor
    inline StandardLocalFunction ( const ThisType &other )
    : discreteFunction_( other.discreteFunction_ ),
      values_( other.values_ ),
      baseFunctionSet_( other.baseFunctionSet_ ),
      entity_( other.entity_ ),
      numScalarDofs_( other.numScalarDofs_ ),
      needCheckGeometry_( other.needCheckGeometry_ )
    {
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) const */
    inline const RangeFieldType &operator[] ( const int num ) const
    {
      assert( entity_ != 0 );
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }

    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) */
    inline RangeFieldType &operator[] ( const int num )
    {
      assert( entity_ != 0 );
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }

    /** \copydoc LocalFunctionInterface::numDofs */
    inline int numDofs() const
    {
      return N * numScalarDofs_; 
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate( const DomainType &x,
                          RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian ( const DomainType &x,
                           JacobianRangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian( const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian( const QuadratureType &quadrature,
                          const int quadPoint,
                          JacobianRangeType &ret ) const;
   
    /** \copydoc Dune::LocalFunctionDefault::baseFunctionSet */
    const BaseFunctionSetType &baseFunctionSet() const
    {
      assert( entity_ != 0 );
      return baseFunctionSet_;
    }
    
    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const RangeType &factor) */
    inline void axpy( const DomainType &x,
                      const RangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const RangeType &factor) */
    template< class QuadratureType >
    inline void axpy( const QuadratureType &quadrature,
                      const int quadPoint,
                      const RangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const JacobianRangeType &factor) */
    inline void axpy( const DomainType &x,
                      const JacobianRangeType &factor );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &factor) */
    template< class QuadratureType >
    inline void axpy( const QuadratureType &quadrature,
                      const int quadPoint,
                      const JacobianRangeType &factor );
    
    /** \copydoc Dune::LocalFunctionInterface::axpy(const DomainType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    inline void axpy( const DomainType &x,
                      const RangeType &factor1,
                      const JacobianRangeType &factor2 );

    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 );

    //! initialize local function 
    void init ( const EntityType &entity );

  protected:
    inline const EntityType &entity () const
    {
      assert( entity_ != 0 );
      return *entity_;
    }

    inline void mapLocalDofs ( const DofStoragePolicyType< PointBased > p,
                               const EntityType &entity,
                               const DiscreteFunctionSpaceType &space )
    {
      const int numDofs = this->numDofs();
      for( int i = 0; i < numDofs; i += N )
      {
        const int baseindex = space.mapToGlobal( entity, i );
        for( int j = 0; j < N; ++j )
        {
          const int dof = i+j;
          const int index = baseindex + j;
          assert( index == space.mapToGlobal( entity, dof ) );
          values_[ dof ] = &(discreteFunction_.dof( index ));
        }
      }
    }

    inline void mapLocalDofs ( const DofStoragePolicyType< VariableBased > p,
                               const EntityType &entity,
                               const DiscreteFunctionSpaceType &space )
    {
      const int numDofs = this->numDofs();
      const int scalarSize = space.containedSpace().size();
      for( int i = 0; i < numDofs; i += N )
      {
        const int baseindex = space.mapToGlobal( entity, i );
        for( int j = 0; j < N; ++j )
        {
          const int dof = i + j;
          const int index = baseindex + j * scalarSize;
          assert( index == space.mapToGlobal( entity, dof ) );
          values_[ dof ] = &(discreteFunction_.dof( index ));
        }
      }
    }

    inline void rightMultiply ( const JacobianRangeType &factor,
                                const DomainType &x,
                                JacobianRangeType &result );
  };



  template< class DiscreteFunctionTraits >
  class StandardLocalFunctionFactory
  {
  private:
    typedef StandardLocalFunctionFactory< DiscreteFunctionTraits > ThisType;

  public:
    typedef typename DiscreteFunctionTraits :: DiscreteFunctionType
      DiscreteFunctionType;

    typedef typename DiscreteFunctionTraits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef StandardLocalFunction< DiscreteFunctionType, DiscreteFunctionSpaceType >
      ObjectType;

  protected:
    DiscreteFunctionType &discreteFunction_;

  public:
    inline explicit StandardLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {
    }

    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_ );
    }
  };
  
}

#include "standardlocalfunction_inline.hh"

#endif
