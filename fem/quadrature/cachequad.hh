#ifndef DUNE_CACHEQUAD_HH
#define DUNE_CACHEQUAD_HH

//- Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/utility/structureutility.hh>

//- Local includes
#include "elementquadrature.hh"
#include "../caching/pointmapper.hh"
#include "../caching/cacheprovider.hh"


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
  template <typename GridImp, int codim>
  class CachingQuadrature
  {
    typedef CompileTimeChecker<false> Only_specialisations_for_codim_0_and_1_so_far;
  };
  
  //! \brief Specialisation for codimension 0.
  template <typename GridImp>
  class CachingQuadrature<GridImp, 0> : public ElementQuadrature<GridImp, 0>
  {
    typedef ElementQuadrature<GridImp, 0> BaseType;

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
    CachingQuadrature(Entity& en, int order) :
      BaseType(en, order)
    {}

    /* Probably not needed, though correct
    using BaseType::nop;
    using BaseType::point;
    using BaseType::weight;
    using BaseType::id;
    using BaseType::geometry;
    using BaseType::order;
    using BaseType::elementGeometry;
    */

    //! Additional method to map quadrature points to caching points.
    //! For codim-0 entities, the quadrature points are the same as the
    //! caching points.
    int cachingPoint(int quadraturePoint) const {
      return quadraturePoint;
    }
  };
  
  //! \brief Specialisation for codimension 1.
  //! Codimension one gets a little tricky, especially when face twists
  //! and non-symmetric quadrature rules are employed... But the details
  //! are safely hidden behind the interface and you don't need to bother.
  template <typename GridImp>
  class CachingQuadrature<GridImp, 1> : public ElementQuadrature<GridImp, 1>
  {
    typedef ElementQuadrature<GridImp, 1> BaseType;

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

  public:
    //! Constructor
    //! \param it Intersection iterator.
    //! \param order The desired order of the quadrature.
    //! \param twist Twist of the face (is 0 in structured grids)
    //! \param side Is either INSIDE or OUTSIDE
    CachingQuadrature(const IntersectionIterator& it, 
                      int order, 
                      int twist,
                      typename BaseType::Side side) :
      BaseType(it, order, side),
      mapper_(CacheProvider<GridImp, codimension>::
              getMapper(this->quadImp(), this->elementGeometry(), 
                        this->faceNumber(), twist))
    {}
      
    //! Additional method to map quadrature points to caching points.
    //! For codim-1 entities, the mapping consists of two stages: First,
    //! consider the twist to get the quadrature point number on the reference
    //! element face and then map it to the caching point.
    size_t cachingPoint(size_t quadraturePoint) const {
      assert(quadraturePoint >= 0 && quadraturePoint < (size_t)this->nop());
      return mapper_[quadraturePoint];
    }

  private:
    typedef typename CachingTraits<RealType, dimension>::MapperType MapperType;

  private:
    const MapperType& mapper_;
  };
}

#endif
