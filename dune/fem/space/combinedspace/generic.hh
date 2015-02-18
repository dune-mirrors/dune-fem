#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH

#include <algorithm>
#include <type_traits>

#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // GenericCombinedDiscreteFunctionSpace
    // ------------------------------------

    template< class Traits >
    class GenericCombinedDiscreteFunctionSpace
      : public DiscreteFunctionSpaceDefault< Traits >
    {
      typedef GenericCombinedDiscreteFunctionSpace< Traits > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

      typedef typename Traits::DiscreteFunctionSpaceTupleType DiscreteFunctionSpaceTupleType;

    public:
      // type of i-th contained sub space
      template< int i >
      using SubDiscreteFunctionSpace = typename Traits::template SubDiscreteFunctionSpace< i >;

      //! extract grid informations, it is assumed the both spaces are living on the
      //! same gridPart
      typedef typename Traits::GridPartType GridPartType;
      typedef typename Traits::GridType GridType;
      typedef typename GridPartType::IntersectionType IntersectionType;
      //! extract informations about IndexSet and Iterators
      typedef typename Traits::IndexSetType IndexSetType;
      typedef typename Traits::IteratorType IteratorType;
      //! dimension of the grid (not the world)
      enum { dimension = GridType::dimension };

      //! the underlaying Analytical function space
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits::polynomialOrder };

      //! type of the base function set(s)
      typedef typename Traits::BasisFunctionSetType BasisFunctionSetType;

      //! mapper used to for block vector function
      typedef typename Traits::BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits::localBlockSize };

      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 669;

      //! type of DofManager
      typedef DofManager< GridType > DofManagerType;

      //! default communication interface
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

      //! default communication direction
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      GenericCombinedDiscreteFunctionSpace ( GridPartType &gridPart,
                                             const InterfaceType commInterface = defaultInterface,
                                             const CommunicationDirection commDirection = defaultDirection )
        : BaseType( gridPart, commInterface, commDirection ),
          spaceTuple_( Traits::createSpaces( gridPart, commInterface, commDirection ) ),
          blockMapper_( Traits::getBlockMapper( spaceTuple_ ) )
      {
        DofManagerType::instance( gridPart.grid() ).addIndexSet( *blockMapper_ );
      }

      // prohibit copy constructor and copy assignment
      GenericCombinedDiscreteFunctionSpace ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

      /** \brief Destructor (freeing base functions and mapper)
          \return
       **/
      ~GenericCombinedDiscreteFunctionSpace ()
      {
        DofManagerType::instance( gridPart().grid() ).removeIndexSet( *blockMapper_ );
        Traits::deleteBlockMapper( blockMapper_ );
        Traits::deleteSpaces( spaceTuple_ );
      }

      using BaseType::gridPart;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains */
      bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there
        return blockMapper().contains( codim );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return Traits::continuous( spaceTuple_ );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const
      {
        return Traits::continuous( intersection, spaceTuple_ );
      }

      /** \brief get the type of this discrete function space
          \return DFSpaceIdentifier
       **/
      DFSpaceIdentifier type () const
      {
        return DFSpaceIdentifier( -1 );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template< class Entity >
      int order ( const Entity &entity ) const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return Traits::getBasisFunctionSet( entity, spaceTuple_ );
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
       **/
      BlockMapperType &blockMapper () const
      {
        return *blockMapper_;
      }

      //! obtain the i-th subspace
      template< int i >
      const typename SubDiscreteFunctionSpace< i >::Type& subDiscreteFunctionSpace() const
      {
        return SubDiscreteFunctionSpace< i >::subDiscreteFunctionSpace( spaceTuple_ );
      }

    private:
      //! tuple of spaces
      DiscreteFunctionSpaceTupleType spaceTuple_;

      //! BlockMapper
      BlockMapperType *blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH
