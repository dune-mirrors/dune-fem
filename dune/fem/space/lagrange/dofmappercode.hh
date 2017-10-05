#ifndef DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH
#define DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

// dune-fem includes
#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/compile.hh>


namespace Dune
{

  namespace Fem
  {

    // LagrangeDofMapperCodeFactory
    // ----------------------------

    template< class LagrangePointSetContainer >
    struct LagrangeDofMapperCodeFactory
    {
      explicit LagrangeDofMapperCodeFactory ( const LagrangePointSetContainer &lagrangePointSets )
      : lagrangePointSets_( lagrangePointSets )
      {}

      template< class RefElement,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().size( 0 ) ) >, int >::value, int > = 0,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().type( 0, 0 ) ) >, GeometryType >::value, int > = 0 >
      DofMapperCode operator() ( const RefElement &refElement ) const
      {
        const GeometryType type = refElement.type();
        if( lagrangePointSets_.exists( type ) )
          return compile( refElement, lagrangePointSets_[ type ] );
        else
          return DofMapperCode();
      }

    private:
      const LagrangePointSetContainer &lagrangePointSets_;
    };




    // LagrangeLocalDofMapping
    // -----------------------

    // for non-Cartesian grids check edge twist and reverse DoF order
    // in this case. This is the case for Lagrange Spaces with polOrder > 2.
    //
    // Notice: This is a purely 2-dimensional hack!
    template< class GridPart >
    class LagrangeLocalDofMapping
    {
      static const int dimension = GridPart::dimension;

      struct Mapping
      {
        explicit Mapping ( bool reverse = false ) : reverse_( reverse ) {}

        template< class Iterator, class Functor >
        void operator() ( std::size_t index, unsigned int numDofs, Iterator begin, Iterator end, Functor functor ) const
        {
          if( reverse_ )
          {
            index += numDofs;
            while( begin != end )
              functor( *(begin++), --index );
          }
          else
          {
            while( begin != end )
              functor( *(begin++), index++ );
          }
        }

        bool reverse_;
      };

    public:
      LagrangeLocalDofMapping ( const GridPart &gridPart )
        : localIdSet_( gridPart.grid().localIdSet() )
      {}

      Mapping operator() ( const typename GridPart::template Codim< 0 >::EntityType &element, unsigned int subEntity, unsigned int codim ) const
      {
        if( codim == dimension-1 )
        {
          const auto &refElement = ReferenceElements< typename GridPart::ctype, GridPart::dimension >::general( element.type() );
          assert( refElement.size( subEntity, codim, dimension ) == 2 );
          const int vx[ 2 ] = { refElement.subEntity( subEntity, codim, 0, dimension ), refElement.subEntity( subEntity, codim, 1, dimension ) };
          return Mapping( localIdSet_.subId( gridEntity( element ), vx[ 1 ], dimension ) < localIdSet_.subId( gridEntity( element ), vx[ 0 ], dimension ) );
        }
        else
          return Mapping();
      }

    private:
      const typename GridPart::GridType::LocalIdSet &localIdSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH
