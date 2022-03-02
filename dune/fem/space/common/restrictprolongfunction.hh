#ifndef DUNE_FEM_RESTRICTPROLONGFUNCTION_HH
#define DUNE_FEM_RESTRICTPROLONGFUNCTION_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/grid.hh>
#include <dune/fem/function/localfunction/const.hh>

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
        typedef typename CoarseFunction::DiscreteFunctionSpaceType CoarseSpace;

        ConstLocalFunction< CoarseFunction > coarseLocalFunction( coarseFunction );
        MutableLocalFunction< FineFunction > fineLocalFunction( fineFunction );

        const CoarseSpace &coarseSpace = coarseFunction.space();
        for( const auto& entity : coarseSpace )
        {
          auto cg = bindGuard( coarseLocalFunction, entity );

          if( isDefinedOn( fineFunction, entity ) )
          {
            auto fg = bindGuard( fineLocalFunction, entity );
            fineLocalFunction.assign( coarseLocalFunction );
          }
          else
            hierarchicProlong( coarseLocalFunction, fineLocalFunction );
        }
      }

    private:
      template< class CoarseLocalFunction, class FineLocalFunction >
      void hierarchicProlong ( const CoarseLocalFunction &coarseLocalFunction,
                               FineLocalFunction &fineLocalFunction ) const
      {
        typedef typename CoarseLocalFunction::EntityType Entity;
        typedef typename Entity::HierarchicIterator HierarchicIterator;

        const Entity &parent = coarseLocalFunction.entity();
        const int childLevel = parent.level()+1;

        const HierarchicIterator hend = parent.hend( childLevel );
        for( HierarchicIterator hit = parent.hbegin( childLevel ); hit != hend; ++hit )
        {
          const Entity &child = *hit;
          if( isDefinedOn( fineLocalFunction.discreteFunction(), child ) )
          {
            auto guard = bindGuard( fineLocalFunction, child );
            localRestrictProlong_.prolongLocal( coarseLocalFunction, fineLocalFunction, child.geometryInFather(), true );
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
        typedef typename CoarseFunction::DiscreteFunctionSpaceType CoarseSpace;

        ConstLocalFunction< FineFunction > fineLocalFunction( fineFunction );
        MutableLocalFunction< CoarseFunction > coarseLocalFunction( coarseFunction );

        const CoarseSpace &coarseSpace = coarseFunction.space();
        for( const auto& entity : coarseSpace )
        {
          auto cg = bindGuard( coarseLocalFunction, entity );

          if( isDefinedOn( fineFunction, entity ) )
          {
            auto fg = bindGuard( fineLocalFunction, entity );
            coarseLocalFunction.assign( fineLocalFunction );
          }
          else
            hierarchicRestrict( fineLocalFunction, coarseLocalFunction );
        }
      }

    private:
      template< class FineLocalFunction, class CoarseLocalFunction >
      void hierarchicRestrict ( FineLocalFunction &fineLocalFunction,
                                CoarseLocalFunction &coarseLocalFunction ) const
      {
        typedef typename CoarseLocalFunction::EntityType Entity;
        typedef typename Entity::HierarchicIterator HierarchicIterator;

        const Entity &parent = coarseLocalFunction.entity();
        const int childLevel = parent.level()+1;

        bool initialize = true;
        const HierarchicIterator hend = parent.hend( childLevel );
        for( HierarchicIterator hit = parent.hbegin( childLevel ); hit != hend; ++hit )
        {
          const Entity &child = *hit;
          if( isDefinedOn( fineLocalFunction.discreteFunction(), child ) )
          {
            auto guard = bindGuard( fineLocalFunction, child );
            localRestrictProlong_.restrictLocal( coarseLocalFunction, fineLocalFunction, child.geometryInFather(), initialize );
          }
          else
            DUNE_THROW( GridError, "Cannot restrict over more than one level." );
          initialize = false;
        }
        localRestrictProlong_.restrictFinalize(parent);
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
