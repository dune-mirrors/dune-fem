#ifndef DUNE_FEM_MISC_HASBOUNDARYINTERSECTION_HH
#define DUNE_FEM_MISC_HASBOUNDARYINTERSECTION_HH

namespace Dune
{
  namespace Fem
  {

    template< class GridPart >
    struct HasBoundaryIntersection
    {
      using EntityType = typename GridPart::template Codim<0>::EntityType;
      static bool apply(const EntityType &entity)
      {
        return entity.hasBoundaryIntersections();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_HASBOUNDARYINTERSECTION_HH
