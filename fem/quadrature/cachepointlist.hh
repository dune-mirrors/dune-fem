#ifndef DUNE_CACHEPOINTLIST_HH
#define DUNE_CACHEPOINTLIST_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

namespace Dune {

  class CachingInterface 
  {
  protected:
    // do not create instances of this class 
    CachingInterface () {}
    
  public:
    //! Method to map quadrature points to caching points.
    //! For codim-1 entities, the mapping consists of two stages: First,
    //! consider the twist to get the quadrature point number on the reference
    //! element face and then map it to the caching point.
    size_t cachingPoint(const size_t quadraturePoint) const 
    {
      std::cerr << "CachingInterface::cachingPoint method not implemented! \n";
      DUNE_THROW(NotImplemented,"CachingInterface::cachingPoint must be overloaded!");
    }
  };
  
  //! \brief Quadrature class which enables use of caching in base function
  //! sets.
  //! CachingPointLists are a conceptual extension of an ElementQuadrature:
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
  template <typename GridPartImp, int codim, 
            class IntegrationTraits> 
  class CachingPointList
  {
    typedef CompileTimeChecker<false> Only_specialisations_for_codim_0_and_1_so_far;
  };
  
  //! \brief Specialisation for codimension 0.
  template <typename GridPartImp, class IntegrationTraits>  
  class CachingPointList<GridPartImp,0,IntegrationTraits> : 
    public ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>,
    public CachingInterface 
  {
    typedef ElementIntegrationPointList<GridPartImp,0,IntegrationTraits> BaseType;

    // type of grid 
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
    CachingPointList(const GeometryType& geo, int order) : BaseType(geo, order)
    {
      CacheProvider<GridType, codimension>::registerQuadrature(this->quadImp());
    }

    //! Additional method to map quadrature points to caching points.
    //! For codim-0 entities, the quadrature points are the same as the
    //! caching points.
    size_t cachingPoint(const size_t quadraturePoint) const {
      return quadraturePoint;
    }
  };
  
  //! \brief Specialisation for codimension 1.
  //! Codimension one gets a little tricky, especially when face twists
  //! and non-symmetric quadrature rules are employed... But the details
  //! are safely hidden behind the interface and you don't need to bother.
  template <typename GridPartImp, class IntegrationTraits>  
  class CachingPointList<GridPartImp, 1, IntegrationTraits> : 
    public ElementIntegrationPointList<GridPartImp,1,IntegrationTraits>, 
    public CachingInterface
  {
    typedef ElementIntegrationPointList<GridPartImp,1,IntegrationTraits> BaseType;

    typedef GridPartImp GridPartType; 
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
    typedef BaseType NonConformingQuadratureType; 

    //! type of twist utility 
    typedef TwistUtility<GridType> TwistUtilityType;
  public:
    //! Constructor
    //! \param gridPart grid part to get twist from twist utility 
    //! \param it Intersection iterator.
    //! \param order The desired order of the quadrature.
    //! \param side Is either INSIDE or OUTSIDE
    CachingPointList(const GridPartType & gridPart, 
                     const IntersectionIterator& it, 
                     int order, 
                     typename BaseType::Side side) :
      BaseType(it, order, side),
      mapper_(CacheProvider<GridType, codimension>::
              getMapper(this->quadImp(), this->elementGeometry(), 
                        this->faceNumber(), (side == BaseType :: INSIDE) ? 
                           TwistUtilityType::twistInSelf(gridPart.grid(),it) : 
                           TwistUtilityType::twistInNeighbor(gridPart.grid(),it)))
    {
      // make sure CachingPointList is only created for conforming
      // intersections 
      assert( TwistUtilityType::conforming(gridPart.grid(),it) );
    }

    //! Additional method to map quadrature points to caching points.
    //! For codim-1 entities, the mapping consists of two stages: First,
    //! consider the twist to get the quadrature point number on the reference
    //! element face and then map it to the caching point.
    size_t cachingPoint(const size_t quadraturePoint) const 
    {
      // this makes no sense for usigned ints ;)
      assert(quadraturePoint >= 0);
      assert(quadraturePoint < (size_t)this->nop());

      return mapper_[quadraturePoint];
    }

    //! return local caching point 
    //! for debugging issues only 
    size_t localCachingPoint(size_t quadraturePoint) const 
    {
      // this makes no sense for usigned ints ;)
      assert(quadraturePoint >= 0);
      assert(quadraturePoint < (size_t)this->nop());

      int faceIndex = this->faceNumber();
      int point = mapper_[quadraturePoint] - faceIndex*mapper_.size();
      assert( mapper_[quadraturePoint] >= 0 );

      assert( point < this->nop() );
      return point;
    }

  private:
    typedef typename CachingTraits<RealType, dimension>::MapperType MapperType;

  private:
    const MapperType& mapper_;
  };
}

#endif
