#ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_TUPLE_HH
#define DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_TUPLE_HH

#include <cstddef>

#include <tuple>
#include <utility>
#include <vector>

#include <dune/fem/common/forloop.hh>

#include "dataprojection.hh"
#include "default.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // DataProjectionTuple
      // -------------------

      /** \brief A DataProjection wrapping an arbitrary number of projection operators
       *
       *  \tparam DataProjections  a number of DataProjection objects
       *
       *  \ingroup DiscreteFunctionSpace_RestrictProlong
       */
      template< class... DataProjections >
      class DataProjectionTuple;

      template< class Head, class... Tail >
      class DataProjectionTuple< Head, Tail... >
      : public DataProjection< typename Head::DiscreteFunctionSpaceType, DataProjectionTuple< Head, Tail... > >
      {
        using ThisType = DataProjectionTuple< Head, Tail... >;
        using BaseType = DataProjection< typename Head::DiscreteFunctionSpaceType, DataProjectionTuple< Head, Tail... > >;

      public:
        /** \copydoc Dune::Fem::hpDG::DataProjection::DiscreteFunctionSpaceType */
        using DiscreteFunctionSpaceType = typename BaseType::DiscreteFunctionSpaceType;
        /** \copydoc Dune::Fem::hpDG::DataProjection::BasisFunctionSetType */
        using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;
        /** \copydoc Dune::Fem::hpDG::DataProjection::EntityType */
        using EntityType = typename BaseType::EntityType;

      private:
        template< int i >
        struct Operation
        {
          static void apply ( const EntityType &entity,
                              const BasisFunctionSetType &prior,
                              const BasisFunctionSetType &present,
                              const std::vector< std::size_t > &origin,
                              const std::vector< std::size_t > &destination,
                              std::tuple< Head, Tail... > &tuple )
          {
            std::get< i >( tuple )( entity, prior, present, origin, destination );
          }
        };

        template< int i >
        struct Project
        {
          template <class TemporaryStorage>
          static void apply ( TemporaryStorage& tmp,
                              std::tuple< Head, Tail... > &tuple )
          {
            std::get< i >( tuple )( tmp );
          }
        };

        template< int i >
        struct AddToList
        {
          template< class Communicator >
          static void apply ( Communicator &comm,
                              std::tuple< Head, Tail... > &tuple )
          {
            std::get< i >( tuple ).addToList( comm );
          }
        };

      public:
        /** \name Construction
         *  \{
         */

        DataProjectionTuple ( Head &&head, Tail &&...tail )
          : tuple_( std::forward< Head >( head ), std::forward< Tail >( tail )... )
        {}

        /** \} */

#ifndef DOXYGEN

        DataProjectionTuple ( const ThisType & ) = delete;

        DataProjectionTuple ( ThisType && ) = default;

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
          Dune::Fem::ForLoop< Operation, 0, sizeof...( Tail ) >::apply( entity, prior, present, origin, destination, tuple_ );
        }

        template <class TemporaryStorage>
        void operator () ( TemporaryStorage& tmp )
        {
          Dune::Fem::ForLoop< Project, 0, sizeof...( Tail ) >::apply( tmp, tuple_ );
        }

        /** \copydoc Dune::Fem::Adaptive::DataProjection::addToList () */
        template< class Communicator >
        void addToList ( Communicator &comm )
        {
          Dune::Fem::ForLoop< AddToList, 0, sizeof...( Tail ) >::apply( comm, tuple_ );
        }

      protected:
        std::tuple< Head, Tail... > tuple_;
      };



      // DefaultDataProjectionTuple
      // --------------------------

      /** \brief A DataProjection for managing an arbitrary number of discrete functions
       *
       *  \tparam DiscreteFunctions  a number of discrete functions
       *
       *  \ingroup DiscreteFunctionSpace_RestrictProlong
       */
      template< class... DiscreteFunctions >
      class DefaultDataProjectionTuple
      : public DataProjectionTuple< DefaultDataProjection< DiscreteFunctions >... >
      {
        using BaseType = DataProjectionTuple< DefaultDataProjection< DiscreteFunctions >... >;

      public:
        explicit DefaultDataProjectionTuple ( DiscreteFunctions &... discreteFunctions )
          : BaseType( DefaultDataProjection< DiscreteFunctions >( discreteFunctions )... )
        {}
      };

    } // namespace hpDG

    // forward types to Fem namespace for convenience
    using hpDG::DefaultDataProjectionTuple ;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_TUPLE_HH
