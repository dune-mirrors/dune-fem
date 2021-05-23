#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH

#include <cstddef>

#include <functional>
#include <memory>
#include <type_traits>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/localinterpolation.hh>

#include <dune/fem/space/common/dataprojection.hh>

#include <dune/fem/space/discontinuousgalerkin/interpolation.hh>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // DiscontinuousGalerkinSpace
      // --------------------------

      /** \brief Generic implementation of a \f$hp\f$-adaptive discontinuous finite element space
       *
       *  \tparam Traits  traits class
       *
       *  \ingroup DiscreteFunctionSpace_API
       */
      template< class Traits >
      class DiscontinuousGalerkinSpace
      : public Dune::Fem::DiscreteFunctionSpaceDefault< Traits >
      {
        using BaseType = Dune::Fem::DiscreteFunctionSpaceDefault< Traits >;

      public:
        /** \brief grid part type */
        using GridPartType =  typename BaseType::GridPartType;
        /** \brief entity type */
        using EntityType = typename BaseType::EntityType;

        /** \brief basis function sets type */
        using BasisFunctionSetsType = typename Traits::BasisFunctionSetsType;
        /** \brief key type identifying a basis function set */
        using KeyType = typename BasisFunctionSetsType::KeyType;
        /** \brief basis function set type */
        using BasisFunctionSetType =  typename BaseType::BasisFunctionSetType;

        /** \brief block mapper type */
        using BlockMapperType = typename BaseType::BlockMapperType;

        typedef typename BaseType::AuxiliaryDofsType AuxiliaryDofsType;

    protected:
      struct AuxiliaryDofsFactory
      {
        typedef std::pair< AuxiliaryDofsType, int > ObjectType;

        static ObjectType *createObject ( std::pair< GridPartType *, BlockMapperType * > key )
        {
          return new ObjectType( std::piecewise_construct, std::tie( *key.first, *key.second ), std::make_tuple( -1 ) );
        }

        static void deleteObject ( ObjectType *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, BlockMapperType * >, std::pair< AuxiliaryDofsType, int >, AuxiliaryDofsFactory > AuxiliaryDofsProviderType;



      protected:
        using BaseType::asImp;
        using BaseType::gridPart_;

      private:
        template< class DataProjection >
        struct DataProjectionWrapper;

        typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      public:
        /** \brief local interpolation type  */
        using InterpolationImplType = DiscontinuousGalerkinLocalL2Projection< GridPartType, BasisFunctionSetType >;

        /** \brief local interpolation type  */
        using InterpolationType = LocalInterpolationWrapper< DiscreteFunctionSpaceType >;

        /** \name Construction
         *  \{
         */

        template< class Function >
        DiscontinuousGalerkinSpace ( GridPartType &gridPart, const BasisFunctionSetsType &basisFunctionSets,
                                     const KeyType &value, Function function,
                                     const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                     const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, interface, direction ),
            basisFunctionSets_( basisFunctionSets ),
            blockMapper_( gridPart_, basisFunctionSets_, value, function )
        {}

        DiscontinuousGalerkinSpace ( GridPartType &gridPart, const BasisFunctionSetsType &basisFunctionSets, const KeyType &value,
                                     const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                     const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : DiscontinuousGalerkinSpace( gridPart, basisFunctionSets, value, [&value]( const EntityType &){ return value; }, interface, direction )
        {}

        /** \} */

        /** \name Deleted methods
         *  \{
         */

        /** \brief copy constructor */
        DiscontinuousGalerkinSpace ( const DiscontinuousGalerkinSpace & ) = delete;

        /** \brief assignment operator */
        const DiscontinuousGalerkinSpace &operator= ( const DiscontinuousGalerkinSpace & ) = delete;

        /** \} */

        /** \name Query methods
         *  \{
         */

        /** \brief please doc me */
        bool continuous () const { return false; }

        /** \brief please doc me */
        bool continuous ( const typename BaseType::IntersectionType &intersection ) const
        {
          return false;
        }

        /** \brief please doc me */
        bool multipleBasisFunctionSets () const { return true; }

        /** \} */

        /** \name Basis function set methods
         *  \{
         */

        /** \brief return polynomial order */
        int order () const { return basisFunctionSets().order(); }

        /** \brief return polynomial order */
        int order ( const EntityType &entity ) const
        {
          return basisFunctionSet( entity ).order();
        }

        /** \brief return basis function set */
        BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
        {
          return basisFunctionSets().basisFunctionSet( entity, key( entity ) );
        }

        /** \} */

        /** \name Interpolation
         *  \{
         */

        /** \brief return interpolation
         *
         *  \param[in]  entity  a grid part entity
         *
         *  \returns local interpolation
         */
        InterpolationType interpolation () const
        {
          return InterpolationType( static_cast< const DiscreteFunctionSpaceType& > (*this) );
        }

        /** \brief return interpolation
         *
         *  \param[in]  entity  a grid part entity
         *
         *  \returns local interpolation
         */
        InterpolationImplType interpolation ( const EntityType &entity ) const
        {
          return InterpolationImplType( basisFunctionSet( entity ) );
        }

        /** \brief return interpolation
         *
         *  \param[in]  entity  a grid part entity
         *
         *  \returns local interpolation
         */
        InterpolationImplType localInterpolation ( const EntityType &entity ) const
        {
          return InterpolationImplType( basisFunctionSet( entity ) );
        }

        /** \} */

        /** \name Block mapper
         *  \{
         */

        /** \brief return block mapper */
        BlockMapperType &blockMapper () const { return blockMapper_; }

        /** \} */

        /** \name Adaptation
         *  \{
         */

        /** \brief get identifiying basis function set key assigned to given entity
         *
         *  \param[in]  entity  grid part entity
         *
         *  \returns key
         */
        const KeyType &key ( const EntityType &entity ) const
        {
          return blockMapper().key( entity );
        }

        /** \brief assign new key to given entity
         *
         *  \param[in]  key  key identifying basis function set
         *  \param[in]  entity  grid part entity
         */
        void mark ( const KeyType &key, const EntityType &entity )
        {
          return blockMapper_.mark( key, entity );
        }

        /** \brief get key to be assigned to an entity after next call to adapt()
         *
         *  \param[in]  entity  grid part entity
         *
         *  \returns key
         */
        KeyType getMark ( const EntityType &entity ) const
        {
          return blockMapper().getMark( entity );
        }

        /** \brief please doc me */
        bool adapt () { return blockMapper_.adapt(); }

        /** \brief please doc me */
        template< class DiscreteFunctionSpace, class Implementation >
        bool adapt ( DataProjection< DiscreteFunctionSpace, Implementation > &projection )
        {
          DataProjectionWrapper< DataProjection< DiscreteFunctionSpace, Implementation > > wrapper( basisFunctionSets(), projection );
          return blockMapper_.adapt( wrapper );
        }

        /** \} */

        /** \name Deprecated methods
         *  \{
         */

        /* \brief return space identitifier */
        DFSpaceIdentifier type () const
        {
          DUNE_THROW( NotImplemented, "Method type() not implemented" );
        }

        /** \} */

        /** \name Non-interface methods
         *  \{
         */

        /** \brief return basis function sets */
        const BasisFunctionSetsType &basisFunctionSets () const
        {
          return basisFunctionSets_;
        }

        /** \} */

      protected:
        BasisFunctionSetsType basisFunctionSets_;
        mutable BlockMapperType blockMapper_;
      };



      // DiscontinuousGalerkinSpace::DataProjectionWrapper
      // -------------------------------------------------

      template< class Traits >
      template< class DataProjection >
      struct DiscontinuousGalerkinSpace< Traits >::DataProjectionWrapper
      {
        explicit DataProjectionWrapper ( const BasisFunctionSetsType &basisFunctionSets,
                                         DataProjection &dataProjection )
          : basisFunctionSets_( basisFunctionSets ),
            dataProjection_( dataProjection )
        {}

        DataProjectionWrapper ( const DataProjectionWrapper & ) = default;

        DataProjectionWrapper &operator= ( const DataProjectionWrapper & ) = default;

        void operator() ( const EntityType &entity,
                          const KeyType &prior,
                          const KeyType &present,
                          const std::vector< std::size_t > &origin,
                          const std::vector< std::size_t > &destination )
        {
          dataProjection_.get()( entity, basisFunctionSet( entity, prior ), basisFunctionSet( entity, present ), origin, destination );
        }

      protected:
        BasisFunctionSetType basisFunctionSet ( const EntityType &entity, const KeyType &key ) const
        {
          return basisFunctionSets_.get().basisFunctionSet( entity, key );
        }

        std::reference_wrapper< const BasisFunctionSetsType > basisFunctionSets_;
        std::reference_wrapper< DataProjection > dataProjection_;
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_SPACE_HH
