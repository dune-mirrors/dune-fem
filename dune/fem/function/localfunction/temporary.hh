#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_TEMPORARY_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_TEMPORARY_HH

#include <dune/common/ftraits.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/common/intersectionside.hh>


namespace Dune
{

  namespace Fem
  {

    // internal forward declarations
    // -----------------------------

    template< class DiscreteFunctionSpace, class DoFVector >
    class BasicTemporaryLocalFunction;

    template< class DiscreteFunctionSpace, class Dof >
    class TemporaryLocalFunction;


    /** \ingroup LocalFunction
        \class BasicTemporaryLocalFunction
        \brief A temporary function carrying values for one entity

        A BasicTemporaryLocalFunction is a LocalFunction which is not associated with
        any DiscreteFunction. It can be used when generating discrete functions
        to temporarily store values for one entity.

        \param DiscreteFunctionSpace type of the discrete function space, the
                                     local function shall belong to
        \param DofVector type of vector the degrees of freedom for an entity are stored in
     */
    template< class DiscreteFunctionSpace, class DofVector >
    class BasicTemporaryLocalFunction
    : public LocalFunction < typename DiscreteFunctionSpace :: BasisFunctionSetType, DofVector >
    {
      typedef BasicTemporaryLocalFunction< DiscreteFunctionSpace, DofVector > ThisType;
      typedef LocalFunction < typename DiscreteFunctionSpace :: BasisFunctionSetType, DofVector > BaseType;

    public:
      //! type of the discrete function space
      typedef DiscreteFunctionSpace  DiscreteFunctionSpaceType;

      //! type of Entity
      typedef typename BaseType :: EntityType EntityType;

      //! type of BasisFunctionSet
      typedef typename BaseType :: BasisFunctionSetType BasisFunctionSetType;

      //! type of LocalDofVector
      typedef typename BaseType :: LocalDofVectorType LocalDofVectorType;

      /* \copydoc Dune::Fem::LocalFunction :: localDofVector */
      using BaseType::localDofVector;

      /** \brief constructor creating a local function without binding it to an
                 entity

          Creates the local function without initializing the fields depending on
          the current entity.

          \note Before using the local function it must be initilized by
          \code
          localFunction.init( entity );
          \endcode

          \param[in] dfSpace discrete function space the local function shall
                             belong to
       */
      explicit BasicTemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace,
                                             const LocalDofVectorType &dofVector = LocalDofVectorType() )
      : BaseType( dofVector ),
        dfSpace_( dfSpace )
      {
        localDofVector().reserve( DiscreteFunctionSpaceType::localBlockSize * dfSpace_.blockMapper().maxNumDofs() );
      }

      /** \brief constructor creating a local function and binding it to an
                 entity

          Creates the local function and initilizes the fields depending on the
          current entity. It is not necessary, though allowed, to call init
          before using the discrete function.

          \note The degrees of freedom are not initialized by this function.

          \param[in] dfSpace discrete function space the local function shall
                             belong to
          \param[in] entity  entity for initialize the local function to
       */
      BasicTemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace, const EntityType &entity,
                                    const LocalDofVectorType &dofVector = LocalDofVectorType() )
      : BaseType( dofVector ),
        dfSpace_( dfSpace )
      {
        localDofVector().reserve( DiscreteFunctionSpaceType::localBlockSize * dfSpace_.blockMapper().maxNumDofs() );
        init( entity );
      }

      /** \brief initialize the local function for an entity
       *
          Binds the local function to an entity.

          \note A local function must be initialized to an entity before it can
                be used.

          \note This function can be called multiple times to use the local
                function for more than one entity.

          \param[in] entity entity to bind the local function to
       */
      void init ( const EntityType &entity )
      {
        BaseType::init( dfSpace_.basisFunctionSet( entity ) );
      }

      /** \brief initialize the local function for an entity
       *
          Binds the local function to an entity.

          \note A local function must be initialized to an entity before it can
                be used.

          \note This function can be called multiple times to use the local
                function for more than one entity.

          \param[in] entity entity to bind the local function to
       */
      void bind ( const EntityType &entity ) { init( entity ); }

      /** \brief Unbinds a local function from an entity.
       */
      void unbind ()
      {
        BaseType::unbind();
      }

      /** \brief initialize the local function for an entity adjacent to the
       * intersection
       *
       *  Binds the local function to an entity.
       *
       *  \note A local function must be initialized to an entity before it can
       *        be used.
       *
       *  \note This function can be called multiple times to use the local
       *        function for more than one entity.
       *
       *   \param[in] intersection to bind the local function to
       *              either inside or outside entity
       *   \param[in] side  side of intersection, i.e. in or out
       */
      template <class IntersectionType>
      void bind(const IntersectionType &intersection, IntersectionSide side)
      {
        // store local copy to avoid problems with casting to temporary types
        const EntityType entity = side==IntersectionSide::in? intersection.inside(): intersection.outside();
        bind( entity );
      }

      /** \brief return discrete function space this local function belongs to
       */
      const DiscreteFunctionSpaceType &space() const
      {
        return dfSpace_;
      }

    protected:
      const DiscreteFunctionSpaceType &dfSpace_;
    };


    /** \ingroup LocalFunction
        \class TemporaryLocalFunction
        \brief A temporary function carrying values for one entity

        A TemporaryLocalFunction is a LocalFunction which is not associated with
        any DiscreteFunction. It can be used when generating discrete functions
        to temporarily store values for one entity.

        \note Local DoF numbers correspond directly to array indices. Hence it
        may be more cache efficient to generate a TemporaryLocalFunction and then
        do only one update step on the discrete function's LocalFunction.

        \param DiscreteFunctionSpaceImp type of the discrete function space, the
                                        local function shall belong to
     */

    template< class DiscreteFunctionSpace, class Dof >
    class TemporaryLocalFunction;
  }

  template< class DiscreteFunctionSpace, class Dof  >
  struct FieldTraits< Fem::TemporaryLocalFunction<DiscreteFunctionSpace,Dof> >
  : public FieldTraits< Dof >
  {};

  namespace Fem
  {
    template< class DiscreteFunctionSpace, class Dof = typename DiscreteFunctionSpace::RangeFieldType >
    class TemporaryLocalFunction
    : public BasicTemporaryLocalFunction< DiscreteFunctionSpace, Dune::DynamicVector< Dof > >
    {
      typedef TemporaryLocalFunction< DiscreteFunctionSpace, Dof > ThisType;
      typedef BasicTemporaryLocalFunction< DiscreteFunctionSpace, Dune::DynamicVector< Dof > > BaseType;

    public:
      //! type of Entity
      typedef typename BaseType :: EntityType EntityType;

      //! type of the discrete function space
      typedef typename BaseType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      /** \brief constructor creating a local function without binding it to an
                 entity

          Creates the local function without initializing the fields depending on
          the current entity.

          \note Before using the local function it must be bound by
          \code
          localFunction.bind( entity );
          \endcode

          \param[in] dfSpace discrete function space the local function shall
                             belong to
       */
      explicit TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( dfSpace ) {}

      /** \brief constructor creating a local function and binding it to an
                 entity

          Creates the local function and initializes the fields depending on the
          current entity. It is not necessary, though allowed, to call bind
          before using the discrete function.

          \note The degrees of freedom are not initialized by this function.

          \param[in] dfSpace discrete function space the local function shall
                             belong to
          \param[in] entity  entity for initialize the local function to
       */
      TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace, const EntityType &entity )
      : BaseType( dfSpace, entity ){}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_TEMPORARY_HH
