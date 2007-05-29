#ifndef DUNE_CACHEQUAD_HH
#define DUNE_CACHEQUAD_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

#include "cachepointlist.hh"

namespace Dune {
  
  //! \brief Quadrature class which enables use of caching in base function
  //! sets.
  //! CachingQuadratures are a conceptual extension of an ElementQuadrature:
  //! they depend on the context (the actual entity, intersection etc.) and
  //! in addition they provide a mapping from local quadrature point numbers
  //! on a subentity's reference element to global quadrature point numbers
  //! on a codim-0 reference element. (Consider for instance a quadrature on
  //! a triangular face of a tetrahedron: it provides n local quadrature points
  //! which can lie on one of the four faces of the tetrahedron, resulting in
  //! 4*n global quadrature point. With the information from the intersection,
  //! for every face in the mesh, the mapping from those local quadrature
  //! points to the global quadrature points can be calculated and used to
  //! cache the evaluation of, say, a base function on those global quadrature
  //! points.)
  template <typename GridPartImp, int codim>
  class CachingQuadrature
  {
    typedef CompileTimeChecker<false> Only_specialisations_for_codim_0_and_1_so_far;
  };

  //! \brief Specialisation for codimension 0.
  template <typename GridPartImp>
  class CachingQuadrature<GridPartImp, 0> :
    public CachingPointList<GridPartImp, 0, 
                            ElementQuadratureTraits<GridPartImp,0> >
  {
    typedef ElementQuadratureTraits<GridPartImp,0> IntegrationTraits;
    typedef CachingPointList<GridPartImp, 0, IntegrationTraits> BaseType;

    typedef typename GridPartImp :: GridType GridType;

  public:
    //! Dimension of the world.
    enum { dimension = BaseType::dimension };
    //! The codimension is zero by definition.
    enum { codimension = 0 };

    //! Just another name for double...
    typedef typename BaseType::RealType RealType;
    //! The type of the coordinates in the codim-0 reference element.
    typedef typename BaseType::CoordinateType CoordinateType;
    //! The type of the codim-0 entity.
    typedef typename BaseType::Entity Entity;
    
  public:
    //! Constructor
    CachingQuadrature(const Entity& en, int order) 
      : BaseType(en.geometry().type(), order)
    {
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const 
    {
      return this->quadImp().weight(i);
    }
  };
  
  //! \brief Specialisation for codimension 1.
  //! Codimension one gets a little tricky, especially when face twists
  //! and non-symmetric quadrature rules are employed... But the details
  //! are safely hidden behind the interface and you don't need to bother.
  template <typename GridPartImp>
  class CachingQuadrature<GridPartImp, 1> :
   public CachingPointList<GridPartImp,1, 
                           ElementQuadratureTraits<GridPartImp,1> >
  {
    typedef ElementQuadratureTraits<GridPartImp,1> IntegrationTraits;

    // type of base class 
    typedef CachingPointList<GridPartImp,1,IntegrationTraits> BaseType;

    // type of grid part 
    typedef GridPartImp GridPartType; 
    // type of grid 
    typedef typename GridPartType :: GridType GridType;

  public:
    //! Dimeinsion of the world
    enum { dimension = BaseType::dimension };
    //! The codimension is one by definition
    enum { codimension = 1 };
    //! A double... or whatever your grid wants
    typedef typename BaseType::RealType RealType;
    //! The coordinates of the quadrature points in the codim-0 reference
    //! element
    typedef typename BaseType::CoordinateType CoordinateType;
    //! Type of the intersection iterator
    typedef typename BaseType::IntersectionIterator IntersectionIterator;

    //! type of quadrature used for non-conforming intersections  
    typedef ElementQuadrature<GridPartImp,codimension> NonConformingQuadratureType; 

    //! type of twist utility 
    typedef TwistUtility<GridType> TwistUtilityType;
  public:
    //! Constructor
    //! \param gridPart grid part to get twist from twist utility 
    //! \param it Intersection iterator.
    //! \param order The desired order of the quadrature.
    //! \param side Is either INSIDE or OUTSIDE
    CachingQuadrature(const GridPartType & gridPart, 
                      const IntersectionIterator& it, 
                      int order, 
                      typename BaseType::Side side) :
      BaseType(gridPart,it, order, side)
    {
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const RealType& weight(size_t i) const 
    {
      return this->quadImp().weight(i);
    }
  };
}

#endif
