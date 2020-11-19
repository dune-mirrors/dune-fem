#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH

#include <algorithm>
#include <memory>

#include <dune/common/hybridutilities.hh>

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

      //! the underlaying analytical function space
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;

      //! type of the base function set(s)
      typedef typename Traits::BasisFunctionSetType BasisFunctionSetType;

      //! mapper used to for block vector function
      typedef typename Traits::BlockMapperType BlockMapperType;

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
       *  \param[in]  gridPart       reference to the grid part
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       *
       */
      GenericCombinedDiscreteFunctionSpace ( GridPartType &gridPart,
                                             const InterfaceType commInterface = defaultInterface,
                                             const CommunicationDirection commDirection = defaultDirection )
        : BaseType( gridPart, commInterface, commDirection ),
          spaceTuple_( Traits::createSpaces( gridPart, commInterface, commDirection ) ),
          blockMapper_( Traits::getBlockMapper( spaceTuple_ ) )
      {}

    protected:

      GenericCombinedDiscreteFunctionSpace ( DiscreteFunctionSpaceTupleType &&spaceTuple )
        : BaseType(
            Traits::template SubDiscreteFunctionSpace< 0 >::subDiscreteFunctionSpace( spaceTuple ).gridPart(),
            Traits::template SubDiscreteFunctionSpace< 0 >::subDiscreteFunctionSpace( spaceTuple ).communicationInterface(),
            Traits::template SubDiscreteFunctionSpace< 0 >::subDiscreteFunctionSpace( spaceTuple ).communicationDirection() ),
          spaceTuple_( std::move( spaceTuple ) ),
          blockMapper_( Traits::getBlockMapper( spaceTuple_ ) )
      {}

    public:

      GenericCombinedDiscreteFunctionSpace ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

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
        return Traits::accumulate( spaceTuple_, true, [] ( bool c, const auto &s ) { return c && s.continuous(); } );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const
      {
        return Traits::accumulate( spaceTuple_, true, [ &intersection ] ( bool c, const auto &s ) { return c && s.continuous( intersection ); } );
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
        return Traits::accumulate( spaceTuple_, int( 0 ), [] ( int o, const auto &s ) { return std::max( o, s.order() ); } );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template< class Entity >
      int order ( const Entity &entity ) const
      {
        return Traits::accumulate( spaceTuple_, int( 0 ), [ &entity ] ( int o, const auto &s ) { return std::max( o, s.order( entity ) ); } );
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
      std::unique_ptr< BlockMapperType > blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_GENERIC_HH
