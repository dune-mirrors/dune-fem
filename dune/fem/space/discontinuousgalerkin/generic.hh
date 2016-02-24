#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH

#include <utility>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/storage/singletonlist.hh>

#include "localinterpolation.hh"
#include "localrestrictprolong.hh"

namespace Dune
{

  namespace Fem
  {

    // GenericDiscontinuousGalerkinSpace
    // ---------------------------------

    /** \class GenericDiscontinuousGalerkinSpace
     *
     *  \brief generic implementation of a Discontinuous Galerkin space based
     *         on a fixed family of basis function sets
     *
     *  \tparam  Traits  traits class
     */
    template< class Traits >
    class GenericDiscontinuousGalerkinSpace
    : public DiscreteFunctionSpaceDefault< Traits >
    {
      typedef GenericDiscontinuousGalerkinSpace< Traits > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

    public:
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::codimension */
      static const int codimension = Traits::codimension;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::GridPartType */
      typedef typename BaseType::GridPartType GridPartType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::EntityType */
      typedef typename BaseType::EntityType EntityType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::IntersectionType */
      typedef typename BaseType::IntersectionType IntersectionType;

      /** \brief basis function sets */
      typedef typename Traits::BasisFunctionSetsType BasisFunctionSetsType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BlockMapperType */
      typedef typename BaseType::BlockMapperType BlockMapperType;

    protected:
      using BaseType::asImp;

    public:
      /** \name Construction
       *  \{
       */

      explicit GenericDiscontinuousGalerkinSpace ( GridPartType &gridPart, BasisFunctionSetsType &&basisFunctionSets,
                                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                   const CommunicationDirection commDirection = ForwardCommunication )
      : BaseType( gridPart, commInterface, commDirection ),
        basisFunctionSets_( std::forward< BasisFunctionSetsType >( basisFunctionSets ) ),
        blockMapper_( gridPart )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      GenericDiscontinuousGalerkinSpace ( const ThisType & ) = delete;

      /** \brief move constructor */
      GenericDiscontinuousGalerkinSpace ( ThisType &&other )
        : BaseType( other ),
          basisFunctionSets_( std::move( other.basisFunctionSets_ ) ),
          blockMapper_( std::move( blockMapper_ ) )
      {}

      GenericDiscontinuousGalerkinSpace &operator= ( const ThisType & ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      static DFSpaceIdentifier type () { return DGSpace_id; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return basisFunctionSets_.basisFunctionSet( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      static constexpr bool continuous () { return false; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      static constexpr bool continuous ( const IntersectionType &intersection ) { return false; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return basisFunctionSets_.order(); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order ( const EntityType &entity ) const { return basisFunctionSets_.order( entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { return blockMapper_; }

      /** \} */

    private:
      BasisFunctionSetsType basisFunctionSets_;
      mutable BlockMapperType blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH
