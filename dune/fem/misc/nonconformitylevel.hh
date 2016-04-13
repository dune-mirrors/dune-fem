#ifndef DUNE_FEM_NONCONFORMITYLEVEL_HH
#define DUNE_FEM_NONCONFORMITYLEVEL_HH

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <dune/common/timer.hh>

namespace Dune
{

  namespace Fem
  {

    /** \brief mark entities such that non-conformity is a given number
     *
     *  \param      gridPart         given grid part, i.e. grid to mark for
     *  \param[in]  levelDifference  maximum allowed level difference
     *  \param[in]  verbose          if true some output is given (default is false)
     */
    template <class GridPartType>
    static inline void
    makeNonConformity(GridPartType& gridPart,
                      const int levelDifference,
                      const bool verbose = false)
    {
      // measure time
      Dune::Timer timer;

      // type of our standard grid iterator
      typedef typename GridPartType :: GridType  GridType;
      typedef typename GridPartType :: template Codim<0>::IteratorType IteratorType;
      typedef typename GridPartType :: GridType :: template Codim<0>:: Entity EntityType;
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

      // end iterator
      const IteratorType endit = gridPart.template end<0> ();

      // grid reference
      GridType& grid = gridPart.grid();

      // get number of elements
      const int gridsize = gridPart.indexSet().size(0);

      // allowed level difference
      const int levelAllowed = levelDifference - 1;

      // make non-conformity only of levelDifferenc
      bool finished = false;
      int count = 0;
      while ( ! finished )
      {
        finished = true;
        for(IteratorType it = gridPart.template begin<0>();
            it != endit; ++it)
        {
          // get entity
          const EntityType & en = *it;
          // get marker
          const int enMarker = grid.getMark(en);

#ifndef NDEBUG
          {
            // make sure we have only one level difference
            IntersectionIteratorType endnit = gridPart.iend(en);
            for(IntersectionIteratorType nit = gridPart.ibegin(en);
                nit != endnit; ++nit)
            {
              const typename IntersectionIteratorType::Intersection &intersec = *nit;
              // check level difference
              if(intersec.neighbor())
              {
                int diff = std::abs(intersec.outside().level() - en.level());
                assert( diff <= levelDifference );
                if( diff > levelDifference )
                {
                  std::cerr << "makeNonConformity: " << diff << " level difference to large! \n";
                  abort();
                }
              }
            }
          }
#endif

          // if entity will be refined, check nothing
          if( enMarker > 0 ) continue;

          // make sure we have only one level difference
          IntersectionIteratorType endnit = gridPart.iend(en);
          for(IntersectionIteratorType nit = gridPart.ibegin(en);
              nit != endnit; ++nit)
          {
            const typename IntersectionIteratorType::Intersection &intersec = *nit;
            if(intersec.neighbor())
            {
              assert( enMarker <= 0 );
              EntityType nb = intersec.outside();
              const int nbMarker = grid.getMark(nb);
              const int levelDiff = nb.level() - en.level();

              // if level difference and refine on neighbor also refine here
              if(levelDiff > levelAllowed)
              {
                // get new marker
                const int newMarker = std::max(enMarker,std::max(nbMarker,0));
                // check whether we have to iterate once more
                finished = (enMarker == newMarker) ? finished : false;
                // mark entity with new marker
                grid.mark(newMarker, en);

                // in case of refinement break
                if( newMarker > 0 ) break;
              }
              else if( (levelDiff == 0) && (nbMarker > 0) )
              {
                // get new marker
                const int newMarker = std::max(enMarker,0);
                // check whether we have to iterate once more
                finished = (enMarker == newMarker) ? finished : false;
                // mark entity with new marker
                grid.mark(newMarker, en);

                // in case of refinement break
                if( newMarker > 0 ) break;
              }
            }

          } // end intersections
        } // end element loop

        ++count;
        if(count > gridsize)
        {
          std::cerr << "makeNonConformity: Break Adaptation loop because not terminating! \n";
          break;
        }
      } // end while

      // output time if verbosity mode
      if(verbose)
      {
        std::cout << "Making non-conformity level took ";
        std::cout << timer.elapsed() << " seconds. \n";
      }
    } //  makeNonConformity

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NONCONFORMITYLEVEL_HH
