#ifndef DUNE_FEM_GRIDPART_FILTER_BASICFILTERWRAPPER_HH
#define DUNE_FEM_GRIDPART_FILTER_BASICFILTERWRAPPER_HH

#include <algorithm>
#include <vector>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/filter/filter.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class > class FilterDefaultImplementation;
    template< class, class > class BasicFilterWrapper;


    // BasicFilterWrapperTraits
    // ------------------------

    template< class GridPartImp, class BasicFilterImp >
    struct BasicFilterWrapperTraits
    {
      //! \brief grid part type
      typedef GridPartImp GridPartType;

      //! \brief export basic filter type
      typedef BasicFilterImp BasicFilterType;

      //! \brief filter type
      typedef BasicFilterWrapper< GridPartType, BasicFilterType > FilterType;

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


    // BasicFilterWrapper
    // ------------------

    template< class GridPartImp, class BasicFilterImp >
    class BasicFilterWrapper
    : public FilterDefaultImplementation< BasicFilterWrapperTraits< GridPartImp, BasicFilterImp > >
    {
      // basic filter type
      typedef BasicFilterImp BasicFilterType;

      // type of grid part
      typedef GridPartImp GridPartType;

      // type of traits
      typedef BasicFilterWrapperTraits< GridPartType, BasicFilterType > Traits;

      // this type
      typedef BasicFilterWrapper< GridPartType, BasicFilterType > ThisType;

      // base type
      typedef FilterDefaultImplementation< Traits > BaseType;

      static const int dimension = GridPartType::GridType::dimension;

      static const int nCodim = dimension+1;

      template< int codim, class GridPart, class BasicFilter >
      struct Contains
      {
        // entity type
        typedef typename ThisType::template Codim< codim >::EntityType EntityType;

        // return true, if entity is contained in filtered gridpart
        static inline bool value ( const EntityType & entity,
                                   const GridPart & gridPart,
                                   const BasicFilter & filter,
                                   std::vector< bool > & contains )
        {
          if( contains.size() != size_t(gridPart.indexSet().size(codim)) )
            update< All_Partition >( gridPart, filter, contains );

          // get index of entity
          typedef typename GridPartType::IndexSetType IndexSetType;
          const IndexSetType & indexSet = gridPart.indexSet();
          size_t index = size_t( indexSet.index( entity ) );

          return contains[ index ];
        }

        // update vector
        template< PartitionIteratorType pitype >
        static inline void update ( const GridPart & gridPart,
                                    const BasicFilter & filter,
                                    std::vector< bool > & contains )
        {
          // type of index set
          typedef typename GridPartType::IndexSetType IndexSetType;

          // get index set
          const IndexSetType & indexSet = gridPart.indexSet();

          // resize vector
          contains.resize( indexSet.size( codim ) );

          // fill vector
          std::fill( contains.begin(), contains.end(), false );

         // codim 0 iterator type
          typedef typename GridPart::template Codim< 0 >::template Partition< pitype >::IteratorType IteratorType;

         // traverse grid
          IteratorType it = gridPart.template begin< 0, pitype >();
          const IteratorType end = gridPart.template end< 0, pitype >();
          for( ; it != end; ++it )
          {
            const typename IteratorType::Entity & entity = *it;

            // continue, if codim 0 entity is not contained in filtered grid part
            if( !filter.contains( entity ) )
              continue;

            const int count = entity.subEntities( codim );
            for( int i = 0; i < count; ++i )
            {
              size_t subIndex = size_t( indexSet.subIndex( entity, i , codim ) );
              contains[ subIndex ] = true;
            }
          }
        }
      };

      template< class GridPart, class BasicFilter >
      struct Contains< 0, GridPart, BasicFilter >
      {
        // entity type
        typedef typename ThisType::template Codim< 0 >::EntityType EntityType;

        // call BasicFilter::contains()
        static inline bool value ( const EntityType & entity,
                                   const GridPart &,
                                   const BasicFilter & filter,
                                   std::vector< bool > & )
        {
          return filter.contains( entity );
        }
      };

    public:
      //! \brief type of the filter implementation
      typedef typename Traits::FilterType FilterType;

      template< int cd >
      struct Codim
      {
        typedef typename Traits::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of codim 0 entity
      typedef typename Traits::EntityType EntityType;

      //! \brief constructor
      BasicFilterWrapper ( const GridPartType & gridPart, const BasicFilterType & filter )
      : gridPart_( gridPart ),
        filter_( filter )
      { }

      //! \brief constructor
      template< typename ... Args >
      BasicFilterWrapper ( const GridPartType & gridPart, Args && ... args )
      : gridPart_( gridPart ),
        filter_( args ... )
      { }

      //! \brief copy constructor
      BasicFilterWrapper ( const ThisType & other )
      : gridPart_( other.gridPart_ ),
        filter_( other.filter_ )
      {
        reset();
      }

      //! \brief assignment operator
      ThisType & operator= ( const ThisType & other )
      {
        gridPart_ = other.gridPart_;
        filter_ = other.filter_;
        reset();
        return *this;
      }

      //! \brief default implementation returns contains from neighbor
      template< class Intersection >
      bool interiorIntersection( const Intersection &intersection ) const
      {
        return filter().interiorIntersection( intersection );
      }

      //! \brief returns true if the given entity of the pointer in the domain
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        return Contains< cd, GridPartType, BasicFilterType >::value( entity, gridPart_, filter_, contains_[ cd ] );
      }

      //! \brief returns true if the given entity of the pointer in the domain
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
        return filter().intersectionBoundary( intersection );
      }

      //! \brief returns the boundary id for an intersection
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection & intersection ) const
      {
        return filter().intersectionBoundaryId( intersection );
      }

      //! \brief returns true if for an intersection a neighbor exsits
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & intersection ) const
      {
        return filter().intersectionNeighbor( intersection );
      }

      //! \brief reset cached values
      void reset ()
      {
        for( int codim = 0; codim < nCodim; ++codim )
          contains_[ codim ].clear();
      }

    private:
      // reference to basic filter
      const BasicFilterType & filter () const
      {
        return filter_;
      }

      const GridPartType & gridPart_;
      BasicFilterType filter_;
      mutable std::vector< bool > contains_[ nCodim ];
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_BASICFILTERWRAPPER_HH
