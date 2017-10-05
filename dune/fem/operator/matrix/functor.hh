#ifndef DUNE_FEM_OPERATOR_MATRIX_FUNCTOR_HH
#define DUNE_FEM_OPERATOR_MATRIX_FUNCTOR_HH

#include <utility>


namespace Dune
{

  namespace Fem
  {

    // IndexFunctor
    // ------------

    template< class Functor, class LocalIndex, class GlobalIndex >
    struct IndexFunctor
    {
      constexpr IndexFunctor ( Functor functor, const LocalIndex &localIndex, const GlobalIndex &globalIndex )
        : functor_( functor ),
          localIndex_( localIndex ),
          globalIndex_( globalIndex )
      {}

      template< class LocalKey, class GlobalKey >
      void operator() ( const LocalKey localKey, const GlobalKey &globalKey ) const
      {
        functor_( std::make_pair( localIndex_, localKey ), std::make_pair( globalIndex_, globalKey ) );
      }

    private:
      Functor functor_;
      LocalIndex localIndex_;
      GlobalIndex globalIndex_;
    };


    // PairFunctor
    // -----------

    template< class Mapper, class Entity, class Functor >
    struct PairFunctor
    {
      PairFunctor ( const Mapper &mapper, const Entity &entity, Functor functor )
        : mapper_( mapper ), entity_( entity ), functor_( std::move( functor ) )
      {}

      template< class LocalKey, class GlobalKey >
      void operator() ( const LocalKey &localKey, const GlobalKey &globalKey ) const
      {
        mapper_.mapEach( entity_, IndexFunctor< Functor, LocalKey, GlobalKey >( functor_, localKey, globalKey ) );
      }

    private:
      const Mapper &mapper_;
      const Entity &entity_;
      Functor functor_;
    };


    // makeRowFunctor
    // --------------

    template< class Mapper, class Entity, class Functor >
    PairFunctor< Mapper, Entity, Functor > makePairFunctor ( const Mapper &mapper, const Entity &entity, Functor functor )
    {
      return PairFunctor< Mapper, Entity, Functor >( mapper, entity, std::move( functor ) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_MATRIX_FUNCTOR_HH
