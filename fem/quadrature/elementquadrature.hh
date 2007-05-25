#ifndef DUNE_ELEMENTQUADRATURE_HH
#define DUNE_ELEMENTQUADRATURE_HH

#include "quadrature.hh"
#include "elementpointlist.hh"

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
  template <typename GridPartImp, int codim>
  class ElementQuadrature {
    typedef CompileTimeChecker<false> Only_implementation_for_codim_0_and_1_exists; 
  };

  //! \brief Element quadrature on codim-0 entities.
  //! For codim-0 element quadratures, there is no additional information
  //! from the context needes, in consequence, the quadrature behaves like
  //! a generic quadrature class, independent from the situation in the grid.
  template <typename GridPartImp>
  class ElementQuadrature<GridPartImp, 0> : public
     ElementIntegrationPointList<GridPartImp,0,Quadrature>
  {
    typedef typename GridPartImp :: GridType GridType;
  public:
    //! Dimension of the world
    enum { dimension = GridType::dimension };
    //! Codimension is zero by definition
    enum { codimension = 0 };
    
    //! The type for reals (mostly double)
    typedef typename GridType::ctype RealType;
    
    typedef ElementIntegrationPointList<GridPartImp,0,Quadrature> BaseType;
  public:
    //! Type for coordinates in the codim-0 reference element 
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    
    //! Type of the codim-0 entity
    typedef typename GridType::template Codim<0>::Entity Entity;

  public:
    //! Constructor
    //! \param en Entity the quadrature lives on (respectively on its reference element).
    //! \param order Desired minimal order of the quadrature.
    ElementQuadrature(const Entity& en, int order) :
      BaseType(en, order)
    {}
    
    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const {
      return this->quadImp().weight(i);
    }
  };


  //! \brief Element quadrature on codim-1 entities.
  //! For codimension 1, the quadrature needs information about the 
  //! intersection. Plus, the user must decide if the quadrature shall live
  //! on the reference element of the outside or inside element of the 
  //! intersection.
  template <typename GridPartImp>
  class ElementQuadrature<GridPartImp, 1> : public 
     ElementIntegrationPointList<GridPartImp,1,Quadrature>
  {
    typedef ElementQuadrature<GridPartImp, 1> ThisType;
    typedef GridPartImp GridPartType;
    typedef typename GridPartType :: GridType GridType;

    typedef ElementIntegrationPointList<GridPartImp,1,Quadrature> BaseType;
  public:
    //! Dimension of the world
    enum { dimension = GridType::dimension };
    //! The codimension is one by definition
    enum { codimension = 1 };
    
    //! Type of the reals (just a fancy name for a double...)
    typedef typename GridType::ctype RealType;

    //! Type of coordinates in codim-0 reference element
    typedef typename Quadrature<RealType, dimension>::CoordinateType CoordinateType;
    //! Type of coordinate in codim-1 reference element
    typedef typename Quadrature<
      RealType, dimension-codimension>::CoordinateType LocalCoordinateType;

    //! Type of the intersection iterator
    typedef typename GridPartImp::IntersectionIteratorType IntersectionIterator;

    //! specify quadrature for use on conforming and non-conforming
    //! intersections 
    typedef ThisType NonConformingQuadratureType;
    
  public:
    //! Constructor
    //! \param gridPart s dummy parameter here 
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param side Is either INSIDE or OUTSIDE
    ElementQuadrature(const GridPartType & gridPart, 
                      const IntersectionIterator& it, 
                      int order, 
                      typename BaseType::Side side) :
      BaseType(gridPart,it,order,side)
    {
    }

    //! Constructor
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param side Is either INSIDE or OUTSIDE
    ElementQuadrature(const IntersectionIterator& it, 
                      int order, 
                      typename BaseType::Side side) :
      BaseType(it,order,side)
    {
    }
    
    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const {
      return this->quadImp().weight(i);
    }
  };

} // end namespace Dune
#endif
