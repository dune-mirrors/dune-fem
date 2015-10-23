#ifndef DUNE_FEM_GRIDPART_COMMON_DEFAULTGRIDPARTENTITY_HH
#define DUNE_FEM_GRIDPART_COMMON_DEFAULTGRIDPARTENTITY_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  namespace Fem
  {

    // DefaultGridPartEntity
    // ---------------------

    template< int codim, int dim, class GridFamily >
    class DefaultGridPartEntity
    {
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = dim;
      static const int mydimension = dimension - codimension;

      int level () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
    };



    // DefaultGridPartEntity for codimension 0
    // ---------------------------------------

    template< int dim, class GridFamily >
    class DefaultGridPartEntity< 0, dim, GridFamily >
    {
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = dim;
      static const int mydimension = dimension - codimension;

      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
      typedef typename Traits::template Codim< codimension >::Entity Entity;

      typedef typename Traits::HierarchicIterator HierarchicIterator;

      int level () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool isLeaf () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      Entity father () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool hasFather () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      const LocalGeometry &geometryInFather () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      HierarchicIterator hbegin ( int maxLevel ) const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      HierarchicIterator hend ( int maxLevel ) const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool isRegular () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool isNew () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool mightVanish () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_DEFAULTGRIDPARTENTITY_HH
