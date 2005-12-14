#ifndef DUNE_CACHINGQUAD_HH
#define DUNE_CACHINGQUAD_HH

//- Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/utility/structureutility.hh>

//- Local includes
#include "quadrature.hh"
#include "../caching/pointmapper.hh"
#include "../caching/twistprovider.hh"
#include "../caching/cacheprovider.hh"

namespace Dune {

  template <typename GridImp, int codim>
  class CacheQuadrature 
  {
    typedef CompileTimeChecker<false> Only_specialisations_for_codim_0_and_1_so_far;
  };

  template <typename GridImp>
  class CacheQuadrature<GridImp, 0> 
  {
  public:
    enum { codim = 0 };

    typedef Quadrature<
      typename GridImp::ctype, GridImp::dimension> QuadratureType;
    typedef typename QuadratureType::CoordinateType CoordinateType;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ct;

  public:
    //! Constructor which initializes the quadrature with a codim-0 entity
    CacheQuadrature(Entity& en, int order) :
      quad_(en.geometry().type(), order) 
    {
      CacheProvider<GridImp, 0>::registerQuadrature(quad_);
    }

    //! The total number of quadrature points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(size_t i) const {
      return quad_.point(i);
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const ct& weight(size_t i) const {
      return quad_.weight(i);
    }

    //! A unique id per quadrature type.
    //! Quadratures are considered as distinct when they differ in the
    //! following points: geometry type, order, dimension and implementation.
    //! \note At the time of writing this, there is only one implementation
    //! per geometry type, order and dimension provided, but the concept is
    //! easily extendible beyond that.
    size_t id() const {
      return quad_.id();
    }

    //! The actual order of the quadrature.
    //! The actual order can be higher as the desired order when no 
    //! implementation for the desired order is found.
    int order() const {
      return quad_.order();
    }

    //! The geometry type the quadrature points belong to.
    GeometryType geo() const {
      return quad_.geo();
    }

    //! The geometry type of the attached element (codim-0 entity)
    GeometryType elementGeo() const {
      return geo();
    }

    //! Additional method to map quadrature points to caching points.
    //! For codim-0 entities, the quadrature points are the same as the
    //! caching points.
    int cachingPoint(int quadraturePoint) const {
      return quadraturePoint;
    }
    
  private:
    QuadratureType quad_;

  };

  //! Specialisation for codim 1
  template <class GridImp>
  class CacheQuadrature<GridImp, 1> 
  {
  private:
    enum { dim = GridImp::dimension };

  public:
    //! The codimension
    enum { codim = 1 };
    //! Attributes the face either to the inside or the outside according
    //! to the definition of the intersection iterator
    enum Side { INSIDE, OUTSIDE };

    typedef Quadrature<typename GridImp::ctype, dim-codim> QuadratureType;

    typedef typename QuadratureType::CoordinateType CoordinateType;
    typedef typename GridImp::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridImp::ctype ct;

  public:
    //! Constructor
    //! Constructs a quadrature on an intersection.
    CacheQuadrature(IntersectionIterator& it, 
                    int order, 
                    int twist, // Move inside once the twist is in interface
                    Side side) :
      quad_(it.intersectionGlobal().type(),
            order),
      faceTwist_(twist), // it.twistSelf(), it.twistNeighbor() later on
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      elementGeo_(side == INSIDE ? 
                  it.inside()->geometry().type() : 
                  it.outside()->geometry().type()),
      // ? replace by pointer?
      mapper_(CacheProvider<GridImp, codim>::
              getMapper(quad_, elementGeo_, faceNumber_, faceTwist_))
    {}

   //! The total number of quadrature points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(int i) const {
      assert(i >= 0 && i < nop());
      return quad_.point(i);
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const ct& weight(int i) const { 
      assert(i >= 0 && i < nop());
      return quad_.weight(i);
    }

    //! A unique id per quadrature type.
    //! Quadratures are considered as distinct when they differ in the
    //! following points: geometry type, order, dimension and implementation.
    //! \note At the time of writing this, there is only one implementation
    //! per geometry type, order and dimension provided, but the concept is
    //! easily extendible beyond that.
    size_t id() const {
      return quad_.id();
    }

    //! The actual order of the quadrature.
    //! The actual order can be higher as the desired order when no 
    //! implementation for the desired order is found.
    int order() const {
      return quad_.order();
    }

    //! The geometry type the face quadrature points belong to.
    GeometryType geo() const {
      return quad_.geo();
    }

    //! The geometry type of the attached element
    GeometryType elementGeo() const {
      return elementGeo_;
    }

    //! Additional method to map quadrature points to caching points.
    //! For codim-1 entities, the mapping consists of two stages: First,
    //! consider the twist to get the quadrature point number on the reference
    //! element face and then map it to the caching point.
    size_t cachingPoint(size_t quadraturePoint) const {
      assert(quadraturePoint >= 0 && quadraturePoint < (size_t)nop());
      return mapper_[quadraturePoint];
    }

  private:
    typedef typename CachingTraits<ct, dim>::MapperType MapperType;

  private:
    QuadratureType quad_;
    int faceTwist_;
    int faceNumber_;
    GeometryType elementGeo_;

    const MapperType& mapper_;
  };

}

#endif
