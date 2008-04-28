#ifndef DUNE_NONCONFORMITYLEVEL_HH
#define DUNE_NONCONFORMITYLEVEL_HH

#include <cassert>
#include <iostream>

//- Dune includes 
#include <dune/common/timer.hh>

namespace Dune
{

/** \brief For a given grid entities are marked such that the level of
 non-conformity is adjusted to a given number, for example 1. 
 \param[inout] gridPart given grid part, i.e. grid to mark for 
 \param[in] levelDifference maximum allowed level difference 
 \param[in] verbose if true some output is given (default is false)
*/
template <class GridPartType>
static inline void
makeNonConformity(GridPartType& gridPart, 
                  const int levelDifference,
                  const bool verbose = false) 
{
  // measure time 
  Timer timer; 

  // type of our standard grid iterator 
  typedef typename GridPartType :: GridType  GridType; 
  typedef typename GridPartType :: template Codim<0>::IteratorType IteratorType;
  typedef typename GridPartType :: GridType :: template Codim<0>:: Entity EntityType;
  typedef typename GridPartType :: GridType :: template Codim<0>:: EntityPointer EntityPointerType;
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
      const int enMarker = grid.getMark(it); 

#ifndef NDEBUG
      {
        // make sure we have only one level difference 
        IntersectionIteratorType endnit = gridPart.iend(en);
        for(IntersectionIteratorType nit = gridPart.ibegin(en);
            nit != endnit; ++nit)
        {
          // check level difference 
          if(nit.neighbor())
          {
            int diff = std::abs(nit.outside().level() - en.level());
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
        if(nit.neighbor())
        {
          assert( enMarker <= 0 );
          EntityPointerType ep = nit.outside();
          const int nbMarker = grid.getMark(ep); 
          const int levelDiff = ep->level() - en.level();

          // if level difference and refine on neighbor also refine here 
          if(levelDiff > levelAllowed) 
          {
            // get new marker 
            const int newMarker = std::max(enMarker,std::max(nbMarker,0));
            // check whether we have to iterate once more 
            finished = (enMarker == newMarker) ? finished : false;
            // mark entity with new marker 
            grid.mark(newMarker, it);

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
            grid.mark(newMarker, it);
        
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
} // end makeNonConformity 

} // end namespace Dune 
#endif
