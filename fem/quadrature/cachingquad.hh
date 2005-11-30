#ifndef DUNE_CACHINGQUAD_HH
#define DUNE_CACHINGQUAD_HH

#include <dune/common/misc.hh>
#include "quadrature.hh"

namespace Dune {

  template <typename GridImp, int codim>
  class CacheQuadrature 
  {
    typedef CompileTimeChecker<false> Only_specialisations_for_codim_0_and_1_so_far;
  };

  template <typename GridImp>
  class CacheQuadrature<GridImp, 0> 
  {
  private:
    typedef Quadrature<
      typename GridImp::ctype, GridImp::dimension> QuadratureType;
    
  private:
    QuadratureType quad_;

  public:
    typedef typename QuadratureType::CoordinateType CoordinateType;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ct;

  public:
    CacheQuadrature(Entity& en, int order) :
      quad_(en.geometry().type(), order)
    {}

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

    //! Additional method to map quadrature points to caching points.
    //! For codim-0 entities, the quadrature points are the same as the
    //! caching points.
    int cachingPoint(int quadraturePoint) const {
      return quadraturePoint;
    }
  };

  template <class GridImp>
  class CacheQuadrature<GridImp, 1> 
  {
  private:
    enum { codim = 1 };
    typedef Quadrature<
      typename GridImp::ctype, GridImp::dimension-codim> QuadratureType;

  private:
    QuadratureType quad_;
    int faceTwist_;
    int faceNumber_;

  public:
    enum Side { INSIDE, OUTSIDE };

    typedef typename QuadratureType::CoordinateType CoordinateType;
    typedef typename GridImp::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridImp::ctype ct;

  public:
    CacheQuadrature(IntersectionIterator& it, 
                    int order, 
                    int twist, // Move inside once the twist is in interface
                    Side side) :
      quad_(side == INSIDE ? 
            it.inside()->geometry().type() : 
            it.outside()->geometry().type(),
            order),
      faceTwist_(twist), // it.twistSelf(), it.twistNeighbor() later on
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor())
    {}

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

    //! Additional method to map quadrature points to caching points.
    //! For codim-1 entities, the mapping consists of two stages: First,
    //! consider the twist to get the quadrature point number on the reference
    //! element face and then map it to the caching point.
    int cachingPoint(int quadraturePoint) const {
      assert(false); // More to do
      return 0;
    }
  };

  // Alternative: Spezialisierung fuer dim, codim
  /*
  template <typename ct, int dim>
  class CachingQuadrature {
  public:
    enum Side {INSIDE, OUTSIDE};
  public:
    template <class Entity>
    CachingQuadrature(Entity& en, int order) :
      quad_(en.geometry().type(), order),
      side_(INSIDE),
      twist_(0),
      caching_(quad_.nop())
    {}

    template <class IntersectionIterator>
    CachingQuadrature(IntersectionIterator& it, int order, 
                      Side side, int twist) :
      quad_(side == INSIDE ? 
            it.inside()->geometry().type() : it.outside()->geometry().type(),
            order),
      side_(side),
      twist_(twist),
      caching_(quad_.nop())
    {}

 
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

    //! Additional method to map quadrature points to caching points.
    int cachingPoint(int quadraturePoint) const {
      assert(quadraturePoint >= 0 && quadraturePoint < nop());
      return caching_[i];
    }

  private:
    Quadrature<ct, dim> quad_;
    Side side_;
    int twist_;

    // To be replaced (Use caching provider here as well)
    std::vector<int> caching_;
  };
  */
}

#endif
