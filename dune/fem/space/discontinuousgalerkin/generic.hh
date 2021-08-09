#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH

#include <utility>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/storage/singletonlist.hh>

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
      typedef SingletonList< const typename GridPartType::IndexSetType*, BlockMapperType > BlockMapperProdiverType;

      typedef CachingQuadrature<GridPartType, EntityType::codimension> VolumeQuadratureType;

    public:
      /** type of discrete function space implementation */
      typedef typename Traits :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

      typedef LocalMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;
      typedef std::vector< typename BaseType::RangeType > VectorType;
      typedef std::pair< LocalMassMatrixType, VectorType >  LocalMassMatrixStorageType;

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
        basisFunctionSets_( std::move( basisFunctionSets ) ),
        // block mapper is a singleton so that the communication can be cached efficiently
        blockMapper_( &BlockMapperProdiverType::getObject( &(gridPart.indexSet() )))
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      GenericDiscontinuousGalerkinSpace ( const ThisType & ) = delete;

      /** \brief move constructor */
      GenericDiscontinuousGalerkinSpace ( ThisType &&other ) = default;

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
        return basisFunctionSets().basisFunctionSet( entity );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      static constexpr bool continuous () { return false; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      static constexpr bool continuous ( const IntersectionType &intersection ) { return false; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return basisFunctionSets().order(); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order ( const EntityType &entity ) const { return basisFunctionSets().order( entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { assert( blockMapper_ ); return *blockMapper_; }

      /** \} */

      const BasisFunctionSetsType &basisFunctionSets () const { return basisFunctionSets_; }
      BasisFunctionSetsType &basisFunctionSets () { return basisFunctionSets_; }

      LocalMassMatrixStorageType& localMassMatrixStorage() const
      {
        auto& localMassPtr = *localMassMatrixStorage_;
        if( ! localMassPtr )
        {
          localMassPtr.reset( new LocalMassMatrixStorageType( LocalMassMatrixType( asImp(), 2*order() ), VectorType() ) );
        }

        return *localMassPtr;
      }

    private:
      BasisFunctionSetsType basisFunctionSets_;
      std::unique_ptr< BlockMapperType, typename BlockMapperProdiverType::Deleter > blockMapper_;

      mutable ThreadSafeValue< std::shared_ptr< LocalMassMatrixStorageType > > localMassMatrixStorage_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_GENERIC_HH
