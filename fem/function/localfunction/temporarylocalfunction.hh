#ifndef DUNE_FEM_TEMPORARYLOCALFUNCTION_HH
#define DUNE_FEM_TEMPORARYLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{
  
  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator = DefaultArrayAllocator >
  class TemporaryLocalFunctionImpl;

  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator = DefaultArrayAllocator >
  class TemporaryLocalFunction;



  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator >
  struct TemporaryLocalFunctionTraits
  {
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
    typedef TemporaryLocalFunctionImpl
      < DiscreteFunctionSpaceType, ArrayAllocator >
      LocalFunctionImpType;
    typedef TemporaryLocalFunction
      < DiscreteFunctionSpaceType, ArrayAllocator >
      LocalFunctionUserType;
  };


  
  /** \ingroup LocalFunction
   *  \class TemporaryLocalFunction
   *  \brief A temporary function carrying values for one entity
   *
   *  A TemporaryLocalFunction is a LocalFunction which is not associated with
   *  any DiscreteFunction. It can be used when generating discrete functions
   *  to temporarily store values for one entity.
   *
   *  \note Local DoF numbers correspond directly to array indices. Hence it
   *  may be more cache efficient to generate a TemporaryLocalFunction and then
   *  do only one update step on the discrete function's LocalFunction.
   *
   *  \param DiscreteFunctionSpaceImp type of the discrete function space, the
   *                                  local function shall belong to
   */
  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator >
  class TemporaryLocalFunction
  : public LocalFunction
    < TemporaryLocalFunctionTraits< DiscreteFunctionSpace, ArrayAllocator > >
  {
  public:
    typedef TemporaryLocalFunctionTraits
      < DiscreteFunctionSpace, ArrayAllocator >
      Traits;

    //! type of the discrete function space
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
    //! type of the local function implementation (engine concept)
    typedef typename Traits :: LocalFunctionImpType LocalFunctionImpType;

  private:
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType, ArrayAllocator >
      ThisType;
    typedef LocalFunction< Traits > BaseType;

    friend class EngineWrapper< LocalFunctionImpType, ThisType >;

  public:
    typedef typename BaseType :: EntityType EntityType;

  protected:
    LocalFunctionImpType impl_;

  public:
    /** \brief constructor creating a local function without binding it to an 
     *         entity
     *
     *  Creates the local function without initializing the fields depending on
     *  the current entity.
     *
     *  \note Before using the local function it must be initilized by
     *  \code
     *  localFunction.init( entity );
     *  \endcode
     *
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     */
    inline explicit
    TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace )
    : impl_( dfSpace )
    {
    }
    
    /** \brief constructor creating a local function and binding it to an
     *         entity
     *
     *  Creates the local function and initilizes the fields depending on the
     *  current entity. It is not necessary, though allowed, to call init
     *  before using the discrete function.
     *
     *  \note The degrees of freedom are not initialized by this function.
     *  
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     *  \param[in] entity  entity for initialize the local function to
     */
    inline TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace,
                                    const EntityType &entity )
    : impl_( dfSpace, entity )
    {
    }

    /** \brief copy constructor
     *
     *  Creates the local function as a copy of the specified one. It is bound
     *  to an entity if and only if the copied local function is bound to an
     *  entity.
     *
     *  \note The degrees of freedom are always copied, even if the copied
     *        local function is not bound to an entity.
     * 
     *  \param[in]  other  TemporaryLocalFunction to copy
     */
    inline TemporaryLocalFunction ( const ThisType &other )
    : impl_( other.impl_ )
    {
    }

  protected:
    const LocalFunctionImpType &asImp () const
    { 
      return impl_;
    } 

    LocalFunctionImpType &asImp () 
    {
      return impl_;
    } 
  };



  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator >
  class TemporaryLocalFunctionImpl
  : public LocalFunctionDefault
    < DiscreteFunctionSpace,
      TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator > >
  {
  public:
    //! type of the discrete function space
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  private:
    typedef TemporaryLocalFunctionImpl< DiscreteFunctionSpaceType, ArrayAllocator >
      ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;
    
  protected:
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of the entity, this local function is associated with
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

    //! type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    //! type of the jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    
    //! field type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

  protected:
    typedef DynamicArray< RangeFieldType, ArrayAllocator > DofArrayType;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const EntityType *entity_;

    BaseFunctionSetType baseFunctionSet_;

    DofArrayType dofs_;

    bool needCheckGeometry_;

  public:
    /** \brief constructor creating a local function without binding it to an 
     *         entity
     *
     *  Creates the local function without initializing the fields depending on
     *  the current entity.
     *
     *  \note Before using the local function it must be initilized by
     *  \code
     *  localFunction.init( entity );
     *  \endcode
     *
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     */
    inline explicit
    TemporaryLocalFunctionImpl ( const DiscreteFunctionSpaceType &dfSpace );
    
    /** \brief constructor creating a local function and binding it to an
     *         entity
     *
     *  Creates the local function and initilizes the fields depending on the
     *  current entity. It is not necessary, though allowed, to call init
     *  before using the discrete function.
     *
     *  \note The degrees of freedom are not initialized by this function.
     *  
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     *  \param[in] entity  entity for initialize the local function to
     */
    inline TemporaryLocalFunctionImpl ( const DiscreteFunctionSpaceType &dfSpace,
                                        const EntityType &entity );

    /** \brief copy constructor
     *
     *  Creates the local function as a copy of the specified one. It is bound
     *  to an entity if and only if the copied local function is bound to an
     *  entity.
     *
     *  \note The degrees of freedom are always copied, even if the copied
     *        local function is not bound to an entity.
     * 
     *  \param[in]  other  TemporaryLocalFunction to copy
     */
    inline TemporaryLocalFunctionImpl ( const ThisType &other );

    /** \copydoc Dune::LocalFunction::operator[]( const int num ) const */
    inline const RangeFieldType &operator[] ( const int num ) const;

    /** \copydoc Dune::LocalFunction::operator[]( const int num ) */
    inline RangeFieldType &operator[] ( const int num );

    /** \copydoc Dune::LocalFunction::baseFunctionSet() const */
    inline const BaseFunctionSetType &baseFunctionSet () const;

    /** \copydoc Dune::LocalFunction::entity() const */
    inline const EntityType &entity () const;

    /** \brief initialize the local function for an entity
     *
     *  Binds the local function to an entity.
     *
     *  \note A local function must be initialized to an entity before it can
     *        be used.
     *        
     *  \note This function can be called multiple times to use the local
     *        function for more than one entity.
     *
     *  \param[in] entity entity to bind the local function to
     */
    inline void init ( const EntityType &entity );

    /** \copydoc Dune::LocalFunction::numDofs() const */
    inline int numDofs () const;
  };


  
  template< class DiscreteFunctionSpace,
            template< class > class ArrayAllocator = DefaultArrayAllocator >
  class TemporaryLocalFunctionFactory
  {
  public:
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  private:
    typedef TemporaryLocalFunctionFactory
      < DiscreteFunctionSpaceType, ArrayAllocator >
      ThisType;

  public:
    typedef TemporaryLocalFunctionImpl< DiscreteFunctionSpaceType, ArrayAllocator >
      ObjectType;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

  public:
    inline explicit
    TemporaryLocalFunctionFactory ( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace )
    {
    }

    inline ObjectType *newObject () const
    {
      return new ObjectType( discreteFunctionSpace_ );
    }
  };

}

#include "temporarylocalfunction_inline.hh"

#endif
