#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ENTITY_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/extendedentity.hh>
#include <dune/fem/gridpart/idgridpart/entity.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartEntity for codimension 0
    // ----------------------------------------

    template< int cd, int dim, class GridFamily >
    class FilteredGridPartEntityCodim0 : public IdEntityBasic< cd, dim, GridFamily >
    {
      static_assert( cd == 0, "This should only be used for codim 0");

      typedef IdEntityBasic< cd, dim, GridFamily > BaseType ;
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;
    protected:
      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

      using BaseType :: data;
    public:
      using BaseType :: hostEntity ;
      using BaseType :: codimension ;

      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      /** \} */

      /** \brief construct a null entity */
      FilteredGridPartEntityCodim0 () = default;

      FilteredGridPartEntityCodim0 ( ExtraData data, HostEntityType hostEntity )
      : BaseType( data, hostEntity )
      {}

      FilteredGridPartEntityCodim0 ( HostEntityType hostEntity )
      : BaseType( ExtraData(), hostEntity )
      {
        assert( false );
        std::abort();
      }

      unsigned int subEntities( const unsigned int codim ) const
      {
        return hostEntity().subEntities( codim );
      }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      template< int codim >
      typename Traits::template Codim< codim >::Entity
      subEntity ( int i ) const
      {
        if constexpr ( codim == 0 )
        {
          return *this;
        }
        else
        {
          return hostEntity().template subEntity< codim >( i );
        }
      }

      bool hasBoundaryIntersections () const
      {
        return data()->hasBoundaryIntersections( hostEntity() );
      }
    };

    // FilteredGridPartEntity (reuse host entity for codim > 0)
    // --------------------------------------------------------

    template <int codim, int dim, class GridFamily >
    struct FilteredGridPartEntitySelector
    {
      typedef typename GridFamily::HostGridPart::template Codim< codim >::Entity type;
    };

    template <int dim, class GridFamily >
    struct FilteredGridPartEntitySelector< 0, dim, GridFamily >
    {
      typedef Dune::ExtendedEntity< 0, dim, GridFamily, FilteredGridPartEntityCodim0 > type;
    };

    template < int codim, int dim, class GridFamily >
    using FilteredGridPartEntity = typename FilteredGridPartEntitySelector< codim, dim, GridFamily >::type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_ENTITY_HH
