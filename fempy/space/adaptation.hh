#ifndef DUNE_FEMPY_SPACE_ADAPTATION_HH
#define DUNE_FEMPY_SPACE_ADAPTATION_HH

#include <utility>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

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
    : public Dune::Fem::hpDG::DataProjection< typename DataProjection::DiscreteFunctionSpaceType, DataProjectionVector< DataProjection > >
    {
      using ThisType = DataProjectionVector< DataProjection >;
      using BaseType = Dune::Fem::hpDG::DataProjection< typename DataProjection::DiscreteFunctionSpaceType, DataProjectionVector< DataProjection > >;

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

      void add( DiscreteFunction& df )
      {
        vec_.emplace_back( DataProjectionImpl( df ) );
      }

      void clear() { vec_.clear(); }

    protected:
      std::vector< DataProjectionImpl > vec_;
    };

    template< class DiscreteFunctionSpace, class DataProjection >
    class SpaceAdaptationManager
      : public Dune::Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjection >
    {
      using BaseType = Dune::Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjection >;
      using BaseType :: dataProjection_;
    public:
      explicit AdaptationManager ( DiscreteFunctionSpaceType &space, DataProjection &&dataProjection )
        : BaseType( space, std::forward< DataProjection > ( dataProjection ) )
      {}

      DataProjection& dataProjection() { return dataProjection_; }
    };



    // SpaceAdaptation
    // --------------

    template< class DF >
    struct SpaceAdaptation
    {
      typedef DF DiscreteFunction;
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType   DiscreteFunctionSpace;
      typedef DataProjectionVector< DiscreteFunction > DataProjectionVector;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit SpaceAdaptation ( DiscreteFunctionSpace& space )
        : space_( space ),
          adaptationManager_( space_, DataProjectionVector() )
      {}

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        // add discrete functions to data projection list
        for( Iterator it = begin; it != end; ++it )
        {
          adaptationManager_.dataProjection().add( *it );
        }

        typedef Fem::hpDG::AdaptationManager< DiscreteFunctionSpace, DataProjectionVector > SpaceAdaptationManager;

        // todo mark space according to markings

        adaptationManager_.adapt();

        // clear list of data projections
        adaptationManager_.dataProjection().clear();
      }

    private:
      DiscreteFunctionSpace& space_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_SPACE_ADAPTATION_HH
