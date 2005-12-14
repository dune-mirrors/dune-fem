#ifndef DUNE_ELEMENTQUADRATURE_HH
#define DUNE_ELEMENTQUADRATURE_HH

#include "quadrature.hh"

namespace Dune {
  template <typename GridImp, int codim>
  class ElementQuadrature {
    typedef CompileTimeChecker<false> Only_implementation_for_codim_1_exists; 
  };

  template <typename GridImp>
  class ElementQuadrature<GridImp, 0>
  {
  public:
    enum { dimension = GridImp::dimension };
    enum { codimension = 0 };
    
    typedef typename GridImp::ctype RealType;
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    typedef typename GridImp::template Codim<0>::Entity Entity;

  public:
    ElementQuadrature(const Entity& en, int order) :
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

    //! Access to the ith quadrature point.
    const CoordinateType& localPoint(size_t i) const {
      return quad_.point(i);
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const {
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
    GeometryType geometry() const {
      return quad_.geometry();
    }
    
    //! The geometry type of the codim 0 element (which is the same)
    GeometryType elementGeometry() const {
      return quad_.geometry();
    }

  private:
    Quadrature<RealType, dimension> quad_;
  };

  template <typename GridImp>
  class ElementQuadrature<GridImp, 1> 
  {
  public:
    enum { dimension = GridImp::dimension };
    enum { codimension = 1 };
    
    enum Side { INSIDE, OUTSIDE };

    typedef typename GridImp::ctype RealType;
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    typedef typename Quadrature<
      RealType, dimension-codimension>::CoordinateType LocalCoordinateType;
    typedef typename GridImp::Traits::IntersectionIterator IntersectionIterator;

  public:
    ElementQuadrature(const IntersectionIterator& it, int order, Side side) :
      quad_(it.intersectionGlobal().type(), order),
      referenceGeometry_(side == INSIDE ?
                         it.intersectionSelfLocal() : 
                         it.intersectionNeighborLocal()),
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      dummy_(0.)
    {}

    //! The total number of quadrature points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(size_t i) const {
      dummy_ = referenceGeometry_.global(quad_.point(i));
      return dummy_;
    }

    const LocalCoordinateType& localPoint(size_t i) const {
      return quad_.point(i);
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const {
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
    GeometryType geometry() const {
      return quad_.geo();
    }

    //! The geometry type of the codim 0 reference element.
    GeometryType elementGeometry() const {
      return referenceGeometry_.type();
    }

  protected:
    int faceNumber() const { return faceNumber_; }

    const Quadrature<RealType, dimension-codimension>& quadImp() const
    { 
      return quad_; 
    }

  private:
   typedef typename IntersectionIterator::LocalGeometry ReferenceGeometry;

  private:
    Quadrature<RealType, dimension-codimension> quad_;
    const ReferenceGeometry& referenceGeometry_;
    int faceNumber_;

    mutable CoordinateType dummy_;
  };


}

#endif
