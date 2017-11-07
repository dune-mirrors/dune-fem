#ifndef DUNE_FEM_RESTRICTPROLONGFUNCTION_HH
#define DUNE_FEM_RESTRICTPROLONGFUNCTION_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/grid.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class   ProlongFunction
     *  \ingroup Adaptation
     *  \brief   prolong discrete functions between grid levels
     *
     *  \tparam LRP local restriction and prolongation operator
     *              (e.g., LocalLagrangeRestrictProlong)
     */
    template< class LRP >
    struct ProlongFunction
    {
      //! type of the local restriction and prolongation operator
      typedef LRP LocalRestrictProlong;

      /** \brief prolong a discrete function to finer grid level
       *
       *  \note The grid parts modelling the levels need not be of same type.
       *
       *  \param[in]  coarseFunction  discrete function on the coarser level
       *  \param[out] fineFunction    discrete function on the finer level
       */
      template< class CoarseFunction, class FineFunction >
      void operator() ( const CoarseFunction &coarseFunction,
                        FineFunction &fineFunction ) const
      {
        typedef ConstLocalFunction<CoarseFunction> CoarseLocalFunction;
        typedef typename CoarseFunction::DiscreteFunctionSpaceType CoarseSpace;

        typedef LocalContribution< FineFunction, Assembly::Set > FineLocalFunction;

        const CoarseSpace &coarseSpace = coarseFunction.space();
        CoarseLocalFunction coarseLocalFunction(coarseFunction);
        FineLocalFunction fineLocalFunction(fineFunction);
        for( const auto& entity : coarseSpace )
        {
          auto coarseBind = bindGuard(coarseLocalFunction, entity);

          if( isDefinedOn( fineFunction, entity ) )
          {
            auto fineBind = bindGuard(fineLocalFunction,entity);
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
        typedef typename CoarseLocalFunction::EntityType Entity;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        typedef LocalContribution< FineFunction, Assembly::Set > FineLocalFunction;

        const Entity &father = coarseLocalFunction.entity();
        const int childLevel = father.level()+1;

        FineLocalFunction fineLocalFunction(fineFunction);
        const HierarchicIterator hend = father.hend( childLevel );
        for( HierarchicIterator hit = father.hbegin( childLevel ); hit != hend; ++hit )
        {
          const Entity &son = *hit;
          if( isDefinedOn( fineFunction, son ) )
          {
            auto fineBind = bindGuard(fineLocalFunction,son);
            localRestrictProlong_.prolongLocal( coarseLocalFunction, fineLocalFunction, son.geometryInFather(), true );
          }
          else
            DUNE_THROW( GridError, "Cannot prolong over more than one level." );
        }
      }

      template< class Function >
      static bool isDefinedOn ( const Function &function, const typename Function::GridPartType::template Codim< 0 >::EntityType &entity )
      {
        typedef typename Function::GridPartType::IndexSetType IndexSet;
        const IndexSet &indexSet = function.gridPart().indexSet();
        return indexSet.contains( entity );
      }

    private:
      LocalRestrictProlong localRestrictProlong_;
    };



    /** \class   RestrictFunction
     *  \ingroup Adaptation
     *  \brief   restrict discrete functions between grid levels
     *
     *  \tparam LRP local restriction and prolongation operator
     *              (e.g., LocalLagrangeRestrictProlong)
     */
    template< class LRP >
    struct RestrictFunction
    {
      //! type of the local restriction and prolongation operator
      typedef LRP LocalRestrictProlong;

    public:
      /** \brief restrict a discrete function to coarser grid level
       *
       *  \note The grid parts modelling the levels need not be of same type.
       *
       *  \param[in]  fineFunction   discrete function on the finer level
       *  \param[out] coarseFunction discrete function on the coarser level
       */
      template< class FineFunction, class CoarseFunction >
      void operator() ( const FineFunction &fineFunction,
                        CoarseFunction &coarseFunction ) const
      {
        typedef LocalContribution< CoarseFunction, Assembly::Set > CoarseLocalFunction;
        typedef typename CoarseFunction::DiscreteFunctionSpaceType CoarseSpace;

        typedef ConstLocalFunction<FineFunction> FineLocalFunction;

        const CoarseSpace &coarseSpace = coarseFunction.space();
        CoarseLocalFunction coarseLocalFunction(coarseFunction);
        FineLocalFunction fineLocalFunction(fineFunction);

        for( const auto& entity : coarseSpace )
        {
          auto coarseBind = bindGuard(coarseLocalFunction,entity);

          if( isDefinedOn( fineFunction, entity ) )
          {
            auto fineBind = bindGuard(fineLocalFunction,entity);
            coarseLocalFunction.assign( fineLocalFunction );
          }
          else
            hierarchicRestrict( fineFunction, coarseLocalFunction );
        }
      }

    private:
      template< class FineFunction, class CoarseLocalFunction >
      void hierarchicRestrict ( const FineFunction &fineFunction,
                                CoarseLocalFunction &coarseLocalFunction ) const
      {
        typedef typename CoarseLocalFunction::EntityType Entity;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        ConstLocalFunction<FineFunction> fineLocalFunction(fineFunction);

        const Entity &father = coarseLocalFunction.entity();
        const int childLevel = father.level()+1;

        bool initialize = true;
        const HierarchicIterator hend = father.hend( childLevel );
        for( HierarchicIterator hit = father.hbegin( childLevel ); hit != hend; ++hit )
        {
          const Entity &son = *hit;
          if( isDefinedOn( fineFunction, son ) )
          {
            auto fineBind = bindGuard(fineLocalFunction,son);
            localRestrictProlong_.restrictLocal( coarseLocalFunction, fineLocalFunction, son.geometryInFather(), initialize );
          }
          else
            DUNE_THROW( GridError, "Cannot restrict over more than one level." );
          initialize = false;
        }
      }

      template< class Function >
      static bool isDefinedOn ( const Function &function, const typename Function::GridPartType::template Codim< 0 >::EntityType &entity )
      {
        typedef typename Function::GridPartType::IndexSetType IndexSet;
        const IndexSet &indexSet = function.gridPart().indexSet();
        return indexSet.contains( entity );
      }

    private:
      LocalRestrictProlong localRestrictProlong_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_RESTRICTPROLONGFUNCTION_HH
