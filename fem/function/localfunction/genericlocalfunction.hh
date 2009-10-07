#ifndef DUNE_FEM_GENERICLOCALFUNCTION_HH
#define DUNE_FEM_GENERICLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{
  
  /** \class GenericLocalFunctionImpl
   *  \brief implementation of a local function based on a generic
   *  localfunction space
   */
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  class GenericLocalFunctionImpl
  : public LocalFunctionDefault
    < DiscreteFunctionSpace,
      GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > >
  {
    typedef GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpace, ThisType > BaseType;

  public:
    //! type of discrete function the local function belongs to
    typedef DiscreteFunction                                         DiscreteFunctionType;

    //! type of  discrete function space the local function belongs to
    typedef DiscreteFunctionSpace                                    DiscreteFunctionSpaceType;

  public:
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType           GridType;
    
    //! type of underlying function space
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
   
    //! field type for domain vectors
    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    //! field type for range vectors
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;
    //! type of domain vectors
    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    //! type of range vectors
    typedef typename FunctionSpaceType :: RangeType                  RangeType;
    //! type of the Jacobian
    typedef typename FunctionSpaceType :: JacobianRangeType          JacobianRangeType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };
    
    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType  BaseFunctionSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity       EntityType;

  public:
    //! constructor
    explicit GenericLocalFunctionImpl ( DiscreteFunctionType &discreteFunction );
    
    //! copy constructor
    GenericLocalFunctionImpl ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
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

  protected:
    DiscreteFunctionType &discreteFunction_;
    
    // std::vector to hold the indices
    std::vector< unsigned int > indices_;

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



  template< class DiscreteFunctionTraits >
  class GenericLocalFunctionFactory
  {
    typedef GenericLocalFunctionFactory< DiscreteFunctionTraits > ThisType;

  public:
    typedef typename DiscreteFunctionTraits :: DiscreteFunctionType  DiscreteFunctionType;

    typedef typename DiscreteFunctionTraits
              :: DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;

    typedef GenericLocalFunctionImpl< DiscreteFunctionType,
                                      DiscreteFunctionSpaceType >    ObjectType;

  public:
    explicit GenericLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {}

    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_ );
    }

  protected:
    DiscreteFunctionType &discreteFunction_;
  };
  
}

#include "genericlocalfunction_inline.hh"

#endif // DUNE_FEM_GENERICLOCALFUNCTION_HH
