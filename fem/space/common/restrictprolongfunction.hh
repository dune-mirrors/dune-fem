#ifndef DUNE_RESTRICTPROLONGFUNCTION_HH
#define DUNE_RESTRICTPROLONGFUNCTION_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/grid.hh>

namespace Dune
{

  template< class LocalRestrictProlong >
  struct ProlongFunction
  {

    typedef typename LocalRestrictProlong::Entity Entity;
    typedef typename Entity::HierarchicIterator HierarchicIterator;

    template< class CoarseFunction, class FineFunction >
    void operator() ( const CoarseFunction &coarseFunction,
                      FineFunction &fineFunction ) const
    {
      typedef typename CoarseFunction::LocalFunctionType CoarseLocalFunction;
      typedef typename CoarseFunction::DiscreteFunctionSpaceType CoarseSpace;
      typedef typename CoarseSpace::IteratorType CoarseIterator;

      typedef typename FineFunction::LocalFunctionType FineLocalFunction;

      const CoarseSpace &coarseSpace = coarseFunction.space();
      const CoarseIterator end = coarseSpace.end();
      for( CoarseIterator it = coarseSpace.begin(); it != end; ++it )
      {
        const Entity &entity = *it;
        CoarseLocalFunction coarseLocalFunction = coarseFunction.localFunction( entity );

        if( isDefinedOn( fineFunction, entity ) )
        {
          FineLocalFunction fineLocalFunction = fineFunction.localFunction( entity );
          fineLocalFunction.assign( coarseLocalFunction );
        }
        else
          hierarchicProlong( coarseLocalFunction, fineFunction );
      }
    }

  private:
    template< class CoarseLocalFunction, class FineFunction >
    void hierarchicProlong ( const CoarseLocalFunction &coarseLocalFunction,
                             FineFunction &fineFunction ) const
    {
      typedef typename FineFunction::LocalFunctionType FineLocalFunction;

      const Entity &father = coarseLocalFunction.entity();
      const int childLevel = father.level()+1;

      const HierarchicIterator hend = father.hend( childLevel );
      for( HierarchicIterator hit = father.hbegin( childLevel ); hit != hend; ++hit )
      {
        const Entity &son = *hit;
        if( isDefinedOn( fineFunction, son ) )
        {
          FineLocalFunction fineLocalFunction = fineFunction.localFunction( son );
          localRestrictProlong_.prolongLocal( coarseLocalFunction, fineLocalFunction );
        }
        else
          DUNE_THROW( GridError, "Cannot prolong over more than one level." );
      }
    }

    template< class Function >
    static bool isDefinedOn ( const Function &function, const Entity &entity )
    {
      typedef typename Function::DiscreteFunctionSpaceType::GridPartType::IndexSetType IndexSet;
      const IndexSet &indexSet = function.space().gridPart().indexSet();
      return indexSet.contains( entity );
    }

  private:
    LocalRestrictProlong localRestrictProlong_;
  };

}

#endif // #ifndef DUNE_RESTRICTPROLONGFUNCTION_HH
