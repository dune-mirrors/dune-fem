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
      : public Dune::Fem::DiscreteFunctionSpaceInterface< Traits >
      {
        using BaseType = Dune::Fem::DiscreteFunctionSpaceInterface< Traits >;

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

        /** \brief communicaton manager type */
        using CommunicationManagerType = typename BaseType::CommunicationManagerType;

      protected:
        using BaseType::asImp;

      private:
        template< class DataProjection >
        struct DataProjectionWrapper;

      public:
        /** \brief local interpolation type  */
        using InterpolationType = DiscontinuousGalerkinLocalInterpolation< BasisFunctionSetsType >;

        /** \name Construction
         *  \{
         */

        template< class Function >
        DiscontinuousGalerkinSpace ( GridPartType &gridPart, const BasisFunctionSetsType &basisFunctionSets,
                                     const KeyType &value, Function function,
                                     const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                     const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : gridPart_( gridPart ),
            types_( gridPart_.get().indexSet() ),
            basisFunctionSets_( basisFunctionSets ),
            blockMapper_( gridPart_.get(), basisFunctionSets_, value, function ),
            interface_( interface ),
            direction_( direction )
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

        /** \brief please doc me */
        bool multipleGeometryTypes () const { return types_.multipleGeomTypes(); }

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
        InterpolationType interpolation ( const EntityType &entity ) const
        {
          return InterpolationType( basisFunctionSet( entity ) );
        }

        /** \brief interpolat given local function
         *
         *  \param[in]  localFunction  local function to interpolate
         *  \param[out]  localDofVector  local dof vector
         */
        template< class LocalFunction, class LocalDofVector >
        DUNE_DEPRECATED
        void interpolate ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
        {
          const EntityType &entity = localFunction.entity();
          const auto interpolation = asImp().interpolation( entity );
          interpolation( localFunction, localDofVector );
        }

        /** \} */

        /** \name Block mapper
         *  \{
         */

        /** \brief return number of dofs */
        int size () const { return BaseType::localBlockSize*blockMapper().size(); }

        /** \brief return block mapper */
        BlockMapperType &blockMapper () const { return blockMapper_; }

        int maxNumDofs () const { return blockMapper().maxNumDofs() * BaseType::localBlockSize; }

        /** \} */

        /** \name Grid part
         *  \{
         */

        /** \brief return grid part */
        GridPartType &gridPart () { return gridPart_.get(); }

        /** \brief return grid part */
        const GridPartType &gridPart () const { return gridPart_.get(); }

        /** \brief return index set */
        const typename BaseType::IndexSetType &indexSet () const
        {
          return gridPart().indexSet();
        }

        /** \brief return iterator */
        typename BaseType::IteratorType begin () const
        {
          return gridPart().template begin< EntityType::codimension >();
        }

        /** \brief return end iterator */
        typename BaseType::IteratorType end () const
        {
          return gridPart().template end< EntityType::codimension >();
        }

        /** \brief return grid */
        typename BaseType::GridType &grid () // DUNE_DEPRECATED
        {
          return gridPart().grid();
        }

        /** \brief return grid */
        const typename BaseType::GridType &grid () const // DUNE_DEPRECATED
        {
          return gridPart().grid();
        }

        /** \brief call function for each grid part element */
        template< class Function >
        DUNE_DEPRECATED
        void forEach ( Function &function ) const
        {
          auto first = gridPart().template begin< EntityType::codimension, InteriorBorder_Partition >();
          auto last = gridPart().template end< EntityType::codimension, InteriorBorder_Partition >();
          for( ; first != last; ++first )
            function( *first );
        }

        /** \} */

        /** \name Parallel interface
         *  \{
         */

        /** \brief please doc me */
        Dune::InterfaceType communicationInterface () const { return interface_; }

        /** \brief please doc me */
        Dune::CommunicationDirection communicationDirection () const { return direction_; }

        /** \brief please doc me */
        const CommunicationManagerType &communicator () const
        {
          if( !communicator_ )
            communicator_.reset( new CommunicationManagerType( asImp(), communicationInterface(), communicationDirection() ) );
          return *communicator_;
        }

        /** \brief please doc me */
        template< class DiscreteFunction >
        void communicate ( DiscreteFunction &discreteFunction ) const
        {
          //using OperationType = typename DiscreteFunction::DiscreteFunctionSpaceType::template CommDataHandle< DiscreteFunction >::OperationType;
          //communicate( discreteFunction, static_cast< OperationType * >( nullptr ) );
        }

        /** \brief please doc me */
        template< class DiscreteFunction, class Operation >
        void communicate ( DiscreteFunction &discreteFunction, const Operation operation ) const
        {
          //communicator().exchange( discreteFunction, static_cast< Operation * >( nullptr ) );
        }

        /** \brief please doc me */
        template< class DiscreteFunction, class Operation >
        typename BaseType::template CommDataHandle< DiscreteFunction, Operation >::Type
        createDataHandle ( DiscreteFunction &discreteFunction, const Operation *operation ) const
        {
          return typename BaseType::template CommDataHandle< DiscreteFunction, Operation >::Type( discreteFunction );
        }

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

        /** \brief get number in grid sequence */
        int sequence () const { return blockMapper().dofManager().sequence(); }

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
        std::reference_wrapper< GridPartType > gridPart_;
        Dune::Fem::AllGeomTypes< typename BaseType::IndexSetType, typename BaseType::GridType > types_;
        BasisFunctionSetsType basisFunctionSets_;
        mutable BlockMapperType blockMapper_;
        Dune::InterfaceType interface_;
        Dune::CommunicationDirection direction_;
        mutable std::unique_ptr< CommunicationManagerType > communicator_;
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
