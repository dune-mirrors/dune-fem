#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_LOCALSUBFUNCTION_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_LOCALSUBFUNCTION_HH

#include <dune/common/typetraits.hh>
#include <dune/common/tuples.hh>
// dune-fem includes
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/function/localfunction/default.hh>


namespace Dune
{

  namespace Fem 
  {

    template< class LFTraits, int nr >
    class LocalSubFunctionImpl;

    template< class LFTraits, int nr >
    class LocalSubFunction;


    // LocalSubFunctionTraits 
    // ----------------------
    template< class LFTraits, int nr >
    struct LocalSubFunctionTraits
    {
      //! extract type of discrete function space
      typedef typename Dune :: tuple_element< nr, typename LFTraits :: DiscreteFunctionSpaceType :: SpaceTupleType >  :: type 
        DiscreteFunctionSpaceType;

      //! type of DoF use in local function 
      typedef typename LFTraits :: DofType DofType;

      typedef LocalSubFunctionImpl < LFTraits, nr > LocalFunctionImpType;
      typedef LocalSubFunction < LFTraits, nr > LocalFunctionUserType;
    };


    // LocalSubFunction
    // ----------------
    
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
    template< class LFTraits, int nr >
    class LocalSubFunction
    : public LocalFunction
      < LocalSubFunctionTraits< LFTraits, nr > >     
    {
    public:
      typedef LocalSubFunctionTraits< LFTraits, nr > Traits;

      //! type of the discrete function space
      typedef typename Traits :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      
      //! type of the local function implementation (engine concept)
      typedef typename Traits :: LocalFunctionImpType LocalFunctionImpType;

    private:
      typedef LocalSubFunction< LFTraits, nr > ThisType;
      typedef LocalFunction< Traits > BaseType;

      friend class EngineWrapper< LocalFunctionImpType, ThisType >;
    public:
      typedef typename BaseType :: EntityType EntityType;

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
      LocalSubFunction ( LocalFunction< LFTraits > &hostLocalFunction )  
      : impl_( hostLocalFunction )
      {}
      
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
      LocalSubFunction ( const ThisType &other )
      : impl_( other.impl_ )
      {}

    protected:
      const LocalFunctionImpType &asImp () const
      { 
        return impl_;
      } 

      LocalFunctionImpType &asImp () 
      {
        return impl_;
      } 

      LocalFunctionImpType impl_;
    };


    template<class LFTraits, int nr >
    class ConstLocalSubFunction 
    : public LocalSubFunction< LFTraits, nr >
    {
      typedef LocalSubFunction< LFTraits, nr > BaseType;
      typedef ConstLocalSubFunction< LFTraits, nr > ThisType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename BaseType :: EntityType EntityType;
      typedef typename BaseType :: RangeFieldType RangeFieldType;
      typedef typename BaseType :: DofType DofType;

      ConstLocalSubFunction ( const LocalFunction< LFTraits > &hostLocalFunction ) 
      : BaseType( const_cast<LocalFunction<LFTraits> &>( hostLocalFunction ) )
      {}

      ConstLocalSubFunction ( const ThisType &other )
      : BaseType( other )
      {}

      using BaseType::operator [];

    protected:
      DofType& operator[] ( const int num ) 
      {
        return asImp()[ num ];
      }

      using BaseType::clear;
      using BaseType::asImp;
      using BaseType::assign;
      using BaseType::operator +=;
      using BaseType::operator -=;
      using BaseType::axpy;
    };



    template<class LFTraits, int nr >
    struct LocalSubFunctionImplTraits
    {
      // here we need the switch between the two spaces
      //! extract type of discrete function space
      typedef typename Dune :: tuple_element< nr, typename LFTraits :: DiscreteFunctionSpaceType :: SpaceTupleType >  :: type 
        DiscreteFunctionSpaceType;
    };


    // LocalSubFunctionImpl
    // --------------------

    template< class LFTraits, int nr >
    class LocalSubFunctionImpl
    : public LocalFunctionDefault
      < typename LocalSubFunctionImplTraits< LFTraits, nr > :: DiscreteFunctionSpaceType, 
        LocalSubFunctionImpl< LFTraits, nr > >
    {
    public:
      //! type of the discrete function space
      typedef typename LocalSubFunctionImplTraits< LFTraits, nr > :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    private:
      typedef LocalSubFunctionImpl< LFTraits, nr > ThisType;
      typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;
      
    protected:
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    public:
      //! type of the base function set
      typedef typename DiscreteFunctionSpaceType :: BasisFunctionSetType
        BasisFunctionSetType;

      //! type of DoF used 
      typedef typename LFTraits :: DofType DofType;

      //! type of the entity, this local function is associated with
      typedef typename DiscreteFunctionSpaceType :: EntityType EntityType;

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
      LocalSubFunctionImpl ( LocalFunction< LFTraits > &hostLocalfunction ) 
      : hostLocalfunction_( hostLocalfunction )
      {}
      
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
      LocalSubFunctionImpl ( const ThisType &other )
      : hostLocalfunction_( other.hostLocalfunction_ )
      {}

      /** \copydoc Dune::Fem::LocalFunction::operator[]( const int num ) const */
      const DofType &operator[] ( const int num ) const 
      { 
        return hostLocalfunction_[ index( num ) ];
      }

      /** \copydoc Dune::Fem::LocalFunction::operator[]( const int num ) */
      DofType &operator[] ( const int num ) 
      {
        return hostLocalfunction_[ index( num ) ];
      } 

      /** \copydoc Dune::Fem::LocalFunction::basisFunctionSet() const */
      const BasisFunctionSetType &basisFunctionSet () const 
      {
        return hostLocalfunction_.basisFunctionSet().template subBasisFunctionSet<nr>();
      }

      /** \copydoc Dune::Fem::LocalFunction::entity() const */
      const EntityType &entity () const { return hostLocalfunction_.entity(); }

      /** \copydoc Dune::Fem::LocalFunction::numDofs() const */
      int numDofs () const { return basisFunctionSet().size(); } 

    protected:
      const int index ( const int num ) const
      {
        return num + hostLocalfunction_.basisFunctionSet().template offset<nr>();
      }

      LocalFunction< LFTraits > &hostLocalfunction_;
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_LOCALSUBFUNCTION_HH
