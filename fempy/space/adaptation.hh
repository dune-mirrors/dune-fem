#ifndef DUNE_FEMPY_SPACE_ADAPTATION_HH
#define DUNE_FEMPY_SPACE_ADAPTATION_HH

#include <utility>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/dataprojection.hh>

#include <dune/python/common/common.hh>

#include <dune/fempy/grid/discretefunctionmanager.hh>
#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/parameter.hh>

namespace Dune
{

  namespace FemPy
  {

    // DataProjectionVector
    // -------------------

    /** \brief A DataProjection wrapping an arbitrary number of projection operators
     *
     *  \tparam DataProjections  a number of DataProjection objects
     *
     *  \ingroup DiscreteFunctionSpace_RestrictProlong
     */
    template< class DataProjection >
    class DataProjectionVector
    : public Dune::Fem::DataProjection< typename DataProjection::DiscreteFunctionSpaceType, DataProjectionVector< DataProjection > >
    {
      using ThisType = DataProjectionVector< DataProjection >;
      using BaseType = Dune::Fem::DataProjection< typename DataProjection::DiscreteFunctionSpaceType, DataProjectionVector< DataProjection > >;

    public:
      /** \copydoc Dune::Fem::hpDG::DataProjection::DiscreteFunctionSpaceType */
      using DiscreteFunctionSpaceType = typename BaseType::DiscreteFunctionSpaceType;
      /** \copydoc Dune::Fem::hpDG::DataProjection::BasisFunctionSetType */
      using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;
      /** \copydoc Dune::Fem::hpDG::DataProjection::EntityType */
      using EntityType = typename BaseType::EntityType;

    public:
      /** \name Construction
       *  \{
       */

      DataProjectionVector ( std::vector< DataProjection >&& vec )
        : vec_( std::forward( vec ) )
      {}

      DataProjectionVector () {}

      /** \} */

#ifndef DOXYGEN

      DataProjectionVector ( const ThisType & ) = delete;

      DataProjectionVector ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = delete;

      ThisType &operator= ( ThisType && ) = default;

#endif // #ifndef DOXYGEN

      /** \copydoc Dune::Fem::hpDG::DataProjection::operator() */
      void operator() ( const EntityType &entity,
                        const BasisFunctionSetType &prior,
                        const BasisFunctionSetType &present,
                        const std::vector< std::size_t > &origin,
                        const std::vector< std::size_t > &destination )
      {
        for( auto& projection : vec_ )
        {
          projection( entity, prior, present, origin, destination );
        }
      }

      template <class TemporaryStorage>
      void operator () ( TemporaryStorage& tmp )
      {
        for( auto& projection : vec_ )
        {
          projection( tmp );
        }
      }

      /** \copydoc Dune::Fem::Adaptive::DataProjection::addToList () */
      template< class Communicator >
      void addToList ( Communicator &comm )
      {
        for( auto& projection : vec_ )
        {
          projection.addToList( comm );
        }
      }

      void add( DataProjection&& dp )
      {
        vec_.emplace_back( std::forward< DataProjection > ( dp ) );
      }

      void clear() { vec_.clear(); }

    protected:
      std::vector< DataProjection > vec_; // ???
    };

    template< class DiscreteFunctionSpace, class DataProjection >
    class SpaceAdaptationManager
      : public Dune::Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjection >
    {
      using BaseType = Dune::Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjection >;
      using BaseType :: dataProjection_;
    public:
      explicit SpaceAdaptationManager ( DiscreteFunctionSpace &space, DataProjection &&dataProjection )
        : BaseType( space, std::forward< DataProjection > ( dataProjection ) )
      {}

      DataProjection& dataProjection() { return dataProjection_; }
    };



    // SpaceAdaptation
    // --------------

    template< class DF >
    struct SpaceAdaptation
    {
      typedef DF DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
      typedef Dune::Fem::DefaultDataProjection< DiscreteFunctionType >     DataProjectionType;
      typedef DataProjectionVector< DataProjectionType >                   DataProjectionVectorType;
      typedef SpaceAdaptationManager< DiscreteFunctionSpaceType, DataProjectionVectorType > AdaptationManager;
      typedef typename DiscreteFunctionSpaceType::EntityType Element;

      explicit SpaceAdaptation ( DiscreteFunctionSpaceType& space )
        : space_( space ),
          adaptationManager_( space_, DataProjectionVectorType() )
      {}

      template< class Marking, class Iterator >
      void adapt (const Marking &marking, Iterator begin, Iterator end )
      {
        // add discrete functions to data projection list
        for( Iterator it = begin; it != end; ++it )
        {
          adaptationManager_.dataProjection().add( DataProjection( *it ) );
        }

        for( element : space_)
          space_( marking(element), element );

        // ??? typedef Fem::hpDG::AdaptationManager< DiscreteFunctionSpaceType, DataProjectionVectorType > SpaceAdaptationManager;

        adaptationManager_.adapt();

        // clear list of data projections
        adaptationManager_.dataProjection().clear();
      }
      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        // add discrete functions to data projection list
        for( Iterator it = begin; it != end; ++it )
        {
          adaptationManager_.dataProjection().add( DataProjection( *it ) );
        }

        // ???? typedef Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjectionVector > SpaceAdaptationManager;

        adaptationManager_.adapt();

        // clear list of data projections
        adaptationManager_.dataProjection().clear();
      }

    private:
      DiscreteFunctionSpaceType& space_;
      AdaptationManager          adaptationManager_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_SPACE_ADAPTATION_HH
