#ifndef DUNE_FEM_THEADFILTER_HH
#define DUNE_FEM_THEADFILTER_HH

//- system includes
#include <algorithm>

//- dune-fem includes
#include <dune/fem/gridpart/filter/filter.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune
{
  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class > class FilterDefaultImplementation;
    template< class > class ThreadFilter;


    // ThreadFilterTraits
    // ------------------------

    template< class GridPartImp >
    struct ThreadFilterTraits
    {
      //! \brief grid part type
      typedef GridPartImp GridPartType;

      //! \brief filter type
      typedef ThreadFilter< GridPartType > FilterType;

      //! \brief entity types
      template < int cd >
      struct Codim 
      {
        //! \brief entity type for given codimension
        typedef typename GridPartType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief entity type for codimension 0 
      typedef typename Codim< 0 >::EntityType EntityType;   
    };


    // ThreadFilter
    // ------------------

    template< class GridPartImp >
    class ThreadFilter
    : public FilterDefaultImplementation< ThreadFilterTraits< GridPartImp > >
    {
      // type of grid part 
      typedef GridPartImp GridPartType;
 
      // type of traits
      typedef ThreadFilterTraits< GridPartType > Traits;

      // this type 
      typedef ThreadFilter< GridPartType > ThisType;
      
      // base type
      typedef FilterDefaultImplementation< Traits > BaseType;

      static const int dimension = GridPartType::GridType::dimension;

      static const int nCodim = dimension+1;

    public:
      //! \brief type of the filter implementation
      typedef typename Traits::FilterType FilterType;

      //! \brief type of index set 
      typedef typename GridPartType :: IndexSetType IndexSetType;
    
      template< int cd >
      struct Codim
      {
        typedef typename Traits::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of codim 0 entity 
      typedef typename Traits::EntityType EntityType;

      typedef MutableArray< int >  ThreadArrayType ;
       
      //! \brief constructor
      ThreadFilter ( const GridPartType & gridPart, 
                     const ThreadArrayType& threadNum,
                     const int thead )
      : indexSet_( gridPart.indexSet() ),
        threadNum_( threadNum ),
        thread_( thead )
      { }

      //! \brief copy constructor
      ThreadFilter ( const ThisType & other )
      : indexSet_( other.indexSet_ ),
        threadNum_( other.threadNum_ ),
        thread_( other.thread_ )
      { 
      }

      //! \brief return false since all interior intersections should be skipped 
      template< class Intersection >
      bool interiorIntersection( const Intersection &intersection ) const
      {
        return false;
      }
    
      //! \brief returns true if the given entity has the correct thread number  
      //! for higher codims false is returned 
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        if( cd == 0 ) 
        {
          return (thread_ == threadNum_[ indexSet_.index( entity ) ]);
        }
        else 
          return false ;
      }

      //! \brief returns true if the given entity has the correct thread number  
      //! for higher codims false is returned 
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        enum { cc = Entity::codimension };
        return contains< cc >( entity );
      }
 
      //! \brief returns true if an intersection is a boundary intersection 
      template< class Intersection >
      bool intersectionBoundary( const Intersection & intersection ) const
      {
        return intersection.boundary();
      }
     
      //! \brief returns the boundary id for an intersection 
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection & intersection ) const
      {
        return intersection.boundaryId();
      }

      //! \brief returns true if for an intersection a neighbor exsits 
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & intersection ) const
      {
        return intersection.neighbor();
      }

    protected:
      const IndexSetType& indexSet_;
      const ThreadArrayType& threadNum_ ;
      const int thread_;
    };

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_THREADFILTER_HH
