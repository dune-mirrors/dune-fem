#ifndef DUNE_FEM_STANDARDLOCALFUNCTION_HH
#define DUNE_FEM_STANDARDLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{
  
  //- Forward declarations of Combined Space 
  template< class, int, DofStoragePolicy >
  class CombinedSpace;



  /** \class StandardLocalFunctionImpl
   *  \brief standard implementation of a local function
   */
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  class StandardLocalFunctionImpl
  : public LocalFunctionDefault
    < DiscreteFunctionSpace,
      StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > >
  {
  public:
    //! type of discrete function the local function belongs to
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of  discrete function space the local function belongs to
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  private:
    typedef StandardLocalFunctionImpl< DiscreteFunctionType, DiscreteFunctionSpaceType >
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
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
#if 0
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;
#endif

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

  public:
    //! constructor
    inline explicit StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction );
    
    //! copy constructor
    inline StandardLocalFunctionImpl ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \copydoc Dune::LocalFunction::operator[](const int num) const */
    inline const RangeFieldType &operator[] ( const int num ) const;

    /** \copydoc Dune::LocalFunction::operator[](const int num) */
    inline RangeFieldType &operator[] ( const int num );

    /** \copydoc Dune::LocalFunction::baseFunctionSet() const */
    inline const BaseFunctionSetType &baseFunctionSet() const;

    /** \copydoc Dune::LocalFunction::entity() const */
    inline const EntityType &entity () const;

    //! initialize local function 
    inline void init ( const EntityType &entity );
    
    /** \copydoc Dune::LocalFunction::numDofs() const */
    inline int numDofs() const;
  };



  /** \copydoc Dune::StandardLocalFunctionImpl
   *
   *  Specialised version for CombinedSpaces
   */
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  class StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
  : public LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >,
      StandardLocalFunctionImpl< DiscreteFunction,
                                 CombinedSpace< ContainedFunctionSpace, N, policy > > >
  {
  public:
    //! type of discrete function the local function belongs to
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of discrete function space the local function belongs to
    typedef CombinedSpace< ContainedFunctionSpace, N, policy >
      DiscreteFunctionSpaceType;

  private:
    typedef StandardLocalFunctionImpl< DiscreteFunctionType, DiscreteFunctionSpaceType >
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
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
#if 0
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;
#endif

  protected:
    typedef typename DiscreteFunctionSpaceType :: ContainedRangeType ScalarRangeType;
    typedef typename DiscreteFunctionSpaceType :: ContainedJacobianRangeType
      ScalarJacobianRangeType;
      
  protected:
    DiscreteFunctionType &discreteFunction_;
    
    // array holding pointer to local dofs 
    DynamicArray< RangeFieldType* > values_;

     // base function set 
    BaseFunctionSetType baseFunctionSet_;

    // actual entity
    const EntityType *entity_;

    // number of local dofs (in the scalar case)
    unsigned int numScalarDofs_;

    bool needCheckGeometry_;

  public:
    //! constructor
    inline StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction );

  private:
    // prohibit copying
    inline StandardLocalFunctionImpl ( const ThisType &other );

    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \copydoc Dune::LocalFunction::operator[](const int num) const */
    inline const RangeFieldType &operator[] ( const int num ) const;

    /** \copydoc Dune::LocalFunction::operator[](const int num) */
    inline RangeFieldType &operator[] ( const int num );
   
    /** \copydoc Dune::LocalFunction::baseFunctionSet() const */
    inline const BaseFunctionSetType &baseFunctionSet () const;

    /** \copydoc Dune::LocalFunction::entity() const */
    inline const EntityType &entity () const;
      
    //! initialize local function 
    inline void init ( const EntityType &entity );

    /** \copydoc Dune::LocalFunction::numDofs() const */
    inline int numDofs () const;
    
    inline int numScalarDofs () const;
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

    typedef StandardLocalFunctionImpl
      < DiscreteFunctionType, DiscreteFunctionSpaceType >
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
