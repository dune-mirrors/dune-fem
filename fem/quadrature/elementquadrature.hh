#ifndef DUNE_ELEMENTQUADRATURE_HH
#define DUNE_ELEMENTQUADRATURE_HH

#include "quadrature.hh"

namespace Dune {
  //! \brief Quadrature on codim-0 reference element.
  //! Dune quadratures are defined on a certain geometry, using local 
  //! coordinates for the quadrature points. Often (especially in face
  //! integrals), the quadrature points need to be transformed to local
  //! coordinates in the codim-0 reference element, for instance in the 
  //! evaluation of base functions. To this end, additional information
  //! about the context, i.e. the face number or other information about the
  //! intersection needs to be known. The element quadratures store exactly
  //! this identical information. As a consequence, ElementQuadrature is
  //! not a generic quadrature anymore, but depends on the situation in
  //! which it is used.
  template <typename GridImp, int codim>
  class ElementQuadrature {
    typedef CompileTimeChecker<false> Only_implementation_for_codim_1_exists; 
  };

  //! \brief Element quadrature on codim-0 entities.
  //! For codim-0 element quadratures, there is no additional information
  //! from the context needes, in consequence, the quadrature behaves like
  //! a generic quadrature class, independent from the situation in the grid.
  template <typename GridImp>
  class ElementQuadrature<GridImp, 0>
  {
  public:
    //! Dimension of the world
    enum { dimension = GridImp::dimension };
    //! Codimension is zero by definition
    enum { codimension = 0 };
    
    enum Side { INSIDE, OUTSIDE };

    //! The type for reals (mostly double)
    typedef typename GridImp::ctype RealType;
    //! Type for coordinates in the codim-0 reference element 
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    //! Type of the codim-0 entity
    typedef typename GridImp::template Codim<0>::Entity Entity;

  public:
    //! Constructor
    //! \param en Entity the quadrature lives on (respectively on its reference element).
    //! \param order Desired minimal order of the quadrature.
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

    //! returns quadraturePoint, to behave like a cahcing qaud without
    //! caching, only works for codim 0
    size_t cachingPoint(size_t quadraturePoint) const {
      return quadraturePoint; 
    }
    
  protected:
    const Quadrature<RealType, dimension>& quadImp() const {
      return quad_;
    }
    
  private:
    Quadrature<RealType, dimension> quad_;
  };


  //! \brief Element quadrature on codim-1 entities.
  //! For codimension 1, the quadrature needs information about the 
  //! intersection. Plus, the user must decide if the quadrature shall live
  //! on the reference element of the outside or inside element of the 
  //! intersection.
  template <typename GridImp>
  class ElementQuadrature<GridImp, 1> 
  {
  public:
    //! Dimension of the world
    enum { dimension = GridImp::dimension };
    //! The codimension is one by definition
    enum { codimension = 1 };
    
    enum Side { INSIDE, OUTSIDE };

    //! Type of the reals (just a fancy name for a double...)
    typedef typename GridImp::ctype RealType;
    //! Type of coordinates in codim-0 reference element
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    //! Type of coordinate in codim-1 reference element
    typedef typename Quadrature<
      RealType, dimension-codimension>::CoordinateType LocalCoordinateType;
    //! Type of the intersection iterator
    typedef typename GridImp::Traits::IntersectionIterator IntersectionIterator;

  public:
    //! Constructor
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param twist the twist of the codim 1 entity
    //! \param side Is either INSIDE or OUTSIDE
    ElementQuadrature(const IntersectionIterator& it, int order, int twist,  Side side) :
      quad_(it.intersectionGlobal().type(), order),
      referenceGeometry_(side == INSIDE ?
                         it.intersectionSelfLocal() : 
                         it.intersectionNeighborLocal()),
      elementGeometry_(referenceGeometry_.type().basicType() ,dimension),
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      dummy_(0.)
    {
      /*
      assert( (side == INSIDE) ? 
          (it.inside ()->geometry().type() == elementGeometry_ ) : 
          (it.outside()->geometry().type() == elementGeometry_ ) );
      */
    }
    //! Constructor
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param side Is either INSIDE or OUTSIDE
    ElementQuadrature(const IntersectionIterator& it, int order, Side side) :
      quad_(it.intersectionGlobal().type(), order),
      referenceGeometry_(side == INSIDE ?
                         it.intersectionSelfLocal() : 
                         it.intersectionNeighborLocal()),
      elementGeometry_(referenceGeometry_.type().basicType() ,dimension),
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      dummy_(0.)
    {
      /*
      assert( (side == INSIDE) ? 
          (it.inside ()->geometry().type() == elementGeometry_ ) : 
          (it.outside()->geometry().type() == elementGeometry_ ) );
      */
    }
    
    //! The total number of quadrature points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(size_t i) const {
      dummy_ = referenceGeometry_.global(quad_.point(i));
      return dummy_;
    }

    //! Access to the ith quadrature point in local (codim-1 reference element)
    //! coordinates
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
      return elementGeometry_;
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
    GeometryType elementGeometry_;
    int faceNumber_;

    mutable CoordinateType dummy_;
  };

} // end namespace Dune

#endif
