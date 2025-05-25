#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ENTITY_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/idgridpart/entity.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartEntity (reuse host entity for codim > 0)
    // --------------------------------------------------------

    template <int codim, int dim, class GridFamily >
    using FilteredGridPartEntity = typename GridFamily::Traits::HostGridPartType::template Codim< codim >::Entity;

    // FilteredGridPartEntity for codimension 0
    // ----------------------------------------

    template< int dim, class GridFamily >
    class FilteredGridPartEntity< 0, dim, GridFamily> : public IdEntityBasic< 0, dim, GridFamily >
    {
      typedef IdEntityBasic< 0, dim, GridFamily > BaseType ;
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;
    protected:
      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

    public:
      using BaseType :: codimension ;

      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      /** \} */

      /** \brief construct a null entity */
      FilteredGridPartEntity () = default;

      FilteredGridPartEntity ( ExtraData data, HostEntityType hostEntity )
      : BaseType( data, hostEntity )
      {}

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
        // TODO, implement correctly.
        return hostEntity().hasBoundaryIntersections();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_ENTITY_HH
