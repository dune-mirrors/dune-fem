#ifndef DUNE_FEM_LUMPING_QUADRATURE_HH
#define DUNE_FEM_LUMPING_QUADRATURE_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

namespace Fem {

/**Define a lumping quadrature for all geometries. Note, however, that
 * this may not make sense for anything else than simplices or maybe
 * hexagonal grids. For simplicial meshes the quadrature formula is
 * exact on linear polynomials and hence the quadrature error is
 * quadratic in the mesh-size. A mass-matrix assembled with the
 * caching quadrature will be diagonal in the context of Lagrange
 * spaces. Generally, it is a bad idea to use mass-lumping for
 * anything else than linear (or maybe bilinear) finite elements.
 *
 * There are probably much more efficient ways to perform
 * mass-lumping. This "quadrature" approach is convenient, however, as
 * it can simply be plugged into existing code by replacing
 * the quadrature rules.
 */
template<class FieldImp, class Topology>
class LumpingQuadrature
  : public QuadratureImp<FieldImp, Topology::dimension>
{
 public:
  typedef FieldImp FieldType;
  typedef Topology TopologyType;
  enum { dimension = TopologyType::dimension };

 private:
  typedef LumpingQuadrature<FieldType, TopologyType> ThisType;
  typedef QuadratureImp<FieldType, dimension> BaseType;

 public:
  typedef typename BaseType::CoordinateType CoordinateType;

 protected:
  static const unsigned int topologyId = TopologyType::id;
  typedef Dune::GenericGeometry::ReferenceDomain<TopologyType> ReferenceDomain;

 public:
  /** \brief constructor filling the list of points and weights.
   *
   *  \param[in]  geometry  geometry type for which a quadrature is desired
   *  \param[in]  order     order, ignored
   *  \param[in]  id        unique identifier, ignored
   */
  LumpingQuadrature(const GeometryType& geometry, int order, int id)
    : BaseType(id)
  {
    for (unsigned i = 0; i < ReferenceDomain::numCorners; ++i)
    {
      CoordinateType pt;
      ReferenceDomain::corner(i, pt);
      this->addQuadraturePoint(pt, ReferenceDomain::template volume<FieldType>() / ReferenceDomain::numCorners);
    }
  }

  /** \copydoc QuadratureImp::geometry
   */
  virtual GeometryType geometryType() const { return GeometryType(topologyId, dimension); }
  /** \copydoc QuadratureImp::order
   */
  virtual int order () const { return 1; }

  //! maximal order of available quadratures
  static std::size_t maxOrder () { return 1; }
};

template<class FieldType, int dimension>
struct DefaultLumpingQuadratureTraits
{
  typedef QuadratureImp<FieldType, dimension> IntegrationPointListType;

  typedef LumpingQuadrature<FieldType, typename GenericGeometry::SimplexTopology<dimension>::type> SimplexQuadratureType;
  typedef LumpingQuadrature<FieldType, typename GenericGeometry::CubeTopology<dimension>::type> CubeQuadratureType;
  typedef LumpingQuadrature<FieldType, typename GenericGeometry::PrismTopology<dimension>::type> PrismQuadratureType;
  typedef LumpingQuadrature<FieldType, typename GenericGeometry::PyramidTopology<dimension>::type> PyramidQuadratureType;
  typedef SimplexQuadratureType PointQuadratureType;
  typedef SimplexQuadratureType LineQuadratureType;
};

template<class GridPartImp, int codim>
struct LumpingQuadratureTraits
{
  // type of a single coordinate
  typedef typename GridPartImp::ctype ctype;

  // dimension of the quadrature
  enum { dimension = GridPartImp::dimension };

  // codimension
  enum { codimension = codim };

  typedef Quadrature<ctype, dimension-codim, DefaultLumpingQuadratureTraits> IntegrationPointListType;

  // type of local coordinate
  typedef typename Quadrature<ctype, dimension, DefaultQuadratureTraits>::CoordinateType CoordinateType;
};

template<class GridPartImp, int codim>
class CachingLumpingQuadrature;

template<typename GridPart>
class CachingLumpingQuadrature<GridPart, 0>
  : public CachingPointList<GridPart, 0, LumpingQuadratureTraits<GridPart, 0> >
{
 public:
  //! type of grid partition
  typedef GridPart GridPartType;

  //! codimension of the element quadrature
  enum { codimension = 0 };

 private:
  typedef LumpingQuadratureTraits<GridPartType, codimension> IntegrationTraits;

  typedef CachingLumpingQuadrature<GridPartType, codimension> ThisType;
  typedef CachingPointList<GridPartType, codimension, IntegrationTraits>
  BaseType;

 public:
  //! dimension of the world
  enum { dimension = BaseType::dimension };

  //! just another name for double
  typedef typename BaseType::RealType RealType;
  //! type of the coordinates in the codim-0 reference element
  typedef typename BaseType::CoordinateType CoordinateType;
  //! type of the quadrature point
  typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
  //! type of iterator
  typedef QuadraturePointIterator< ThisType > IteratorType;

  // for compatibility
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

 protected:
  using BaseType::quadImp;

 public:
  using BaseType::nop;

  /** \brief constructor
   *
   *  \param[in]  entity  entity, on whose reference element the quadratre
   *                      lives
   *
   *  \param[in] order   desired minimal order of the quadrature,
   *                     which is of course fixed at 1. But we allow
   *                     for this parameter in order to plug the
   *                     LumpingQuadrature into generic code which
   *                     normally passes the quadrature order as
   *                     second parameter.
   */
  CachingLumpingQuadrature(const EntityType &entity, int order = 1)
    : BaseType(entity.type(), 1)
  {}

  /** \brief copy constructor
   *
   *  \param[in]  org  element quadrature to copy
   */
  CachingLumpingQuadrature(const ThisType &org) : BaseType(org) {}

  QuadraturePointWrapperType operator[] ( std::size_t i ) const
  {
    return QuadraturePointWrapperType( *this, i );
  }

  IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
  IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

  const RealType &weight (std::size_t ) const
  {
    // all weights have the same value
    return quadImp().weight(0);
  }
};

template<typename GridPart>
class CachingLumpingQuadrature<GridPart, 1>
  : public CachingPointList<GridPart, 1, LumpingQuadratureTraits<GridPart, 1> >
{
 public:
  //! type of grid partition
  typedef GridPart GridPartType;

  //! codimension of the element quadrature
  enum { codimension = 1 };

 private:
  typedef LumpingQuadratureTraits<GridPartType, codimension> IntegrationTraits;

  typedef CachingLumpingQuadrature<GridPartType, codimension> ThisType;
  typedef CachingPointList<GridPartType, codimension, IntegrationTraits> BaseType;

 public:
  //! dimension of the world
  enum { dimension = BaseType::dimension };

  //! just another name for double
  typedef typename BaseType::RealType RealType;

  //! type of the coordinates in the codim-0 reference element
  typedef typename BaseType::CoordinateType CoordinateType;

  //! type of the intersection iterator
  typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  //! type of the quadrature point
  typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
  //! type of iterator
  typedef QuadraturePointIterator< ThisType > IteratorType;

  using BaseType::nop;

 protected:
  using BaseType::quadImp;

 public:
  /** \brief constructor
   *
   *  \param[in]  gridPart      grid partition
   *  \param[in]  intersection  intersection
   *  \param[in]  ignored       desired order of the quadrature, which is ignored here
   *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for
   *                            which the ElementQuadrature shall be created
   */
  CachingLumpingQuadrature(const GridPartType& gridPart,
                           const IntersectionType& intersection,
                           int ignored,
                           typename BaseType::Side side)
    : BaseType(gridPart, intersection, 1, side)
  {}

  /** \brief copy constructor
   *
   *  \param[in]  org  element quadrature to copy
   */
  CachingLumpingQuadrature(const ThisType &org) : BaseType(org) {}

  QuadraturePointWrapperType operator[] ( std::size_t i ) const
  {
    return QuadraturePointWrapperType( *this, i );
  }

  IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
  IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

  const RealType &weight (std::size_t ) const
  {
    // all weights should have the same value
    return quadImp().weight(0);
  }
};

} // Fem

} // Dune


#endif // DUNE_FEM_LUMPING_QUADRATURE_HH
