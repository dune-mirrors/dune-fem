#ifndef DUNE_QUADPROVIDER_HH
#define DUNE_QUADPROVIDER_HH

#include "quadrature.hh"
#include "idprovider.hh"

namespace Dune {

  // Forward declarations
  template <class ct, int dim>
  class QuadratureImp;
  template <class ct, int dim>
  class SimplexQuadrature;
  template <class ct, int dim>
  class CubeQuadrature;
  template <class ct>
  class LineQuadrature;
  template <class ct>
  class TriangleQuadrature;
  template <class ct>
  class TetraQuadrature;
  template <class ct>
  class QuadrilateralQuadrature;
  template <class ct>
  class HexaQuadrature;
  template <class ct>
  class PrismQuadrature;
  template <class ct>
  class PyramidQuadrature;

  //! Utility class that factors out the repetitive lookup process.
  class QuadCreator 
  {
    //! class holding vector with pointer to quadrature objects 
    template<class QuadImp>
    class QuadratureStorage  
    {
      // vector holding pointer to quadratures 
      std::vector<QuadImp*> vec_;
    public:   
      //! constructor creating empty vec of length maxOrder
      QuadratureStorage() : vec_(QuadImp::maxOrder()+1 , 0) {}
      
      //! deletes all quadratures
      ~QuadratureStorage() 
      {
        for(size_t i=0; i<vec_.size(); ++i) 
        {
          delete vec_[i];
        }
      }

      //! return reference to quadrature, if pointer not exists
      //! object is created 
      QuadImp& getQuadrature(size_t order) 
      {
        assert( order <= vec_.size() );
        // if not exists, create quadrature  
        if(!vec_[order]) 
        {
          vec_[order] = new QuadImp(order, IdProvider::instance().newId());
        }
        return *(vec_[order]);    
      }      
    }; // end class QuadratureStorage 
  
  public:
    //! provide quadrature, stores internal singleton holding list 
    //! with existing quadratures
    template <class QuadImp>
    static const QuadImp& provideQuad(int order)
    {
      // vector holding pointers 
      static QuadratureStorage<QuadImp> storage;
      return storage.getQuadrature(order);
    }
  };

  //! QuadratureProvider follows the monostate pattern. It provides a single
  //! point of access (and storage) for the actual implementation of
  //! quadratures so that the number of costly creations can be reduced to
  //! a minimum.
  template <typename ct, int dim, template <class,int> class QuadratureTraits >
  class QuadratureProvider {
  public:
    //! Access to the quadrature implementations.
    //! Implementation: if the desired hasn't been constructed, this is
    //! done during the call. In all subsequent calls, the stored object is
    //! returned.
    static const QuadratureImp<ct, dim>& getQuadrature(GeometryType geo, 
                                                       int order);
  };

  //! Specialisation for dimension 0
  template <typename ct, template <class,int> class QuadratureTraits >
  class QuadratureProvider<ct, 0, QuadratureTraits> 
  {
    typedef QuadratureTraits<ct,0> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: PointQuadratureType PointQuadratureType;
  public:
    //! Access to the quadrature implementations.
    static const QuadratureImp<ct, 0>& getQuadrature(GeometryType geo, 
                                                     int order) 
    {
      return QuadCreator::template provideQuad<PointQuadratureType>(order);
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  }; 

  
  //! Specialisation for dimension 1.
  template <class ct, template <class,int> class QuadratureTraits >
  class QuadratureProvider<ct, 1, QuadratureTraits> 
  {
    typedef QuadratureTraits<ct,1> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: LineQuadratureType LineQuadratureType;
  public:
    //! Access to the quadrature implementations.
    static const QuadratureImp<ct, 1>& getQuadrature(GeometryType geo, 
                                                     int order) {
      assert(geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      return QuadCreator::template provideQuad<LineQuadratureType> (order);
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  }; 

  //! Specialisation for dimension = 2
  template <class ct, template <class,int> class QuadratureTraits >
  class QuadratureProvider<ct, 2, QuadratureTraits> 
  {
    typedef QuadratureTraits<ct,2> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: SimplexQuadratureType SimplexQuadratureType;
    typedef typename QuadratureTraitsType :: CubeQuadratureType CubeQuadratureType;
  public:
    //! Access to the quadrature implemenations.
    static const QuadratureImp<ct, 2>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert( geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      if(geo.isTriangle()) 
        return QuadCreator::template provideQuad<SimplexQuadratureType>(order);
      if(geo.isQuadrilateral())
        return QuadCreator::template provideQuad<CubeQuadratureType>(order);

      DUNE_THROW(RangeError, "Element type not available for dim == 2");
      // dummy return
      return QuadCreator::template provideQuad<SimplexQuadratureType>(0);
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  };

  //! Specialisation for dimension = 3.
  template <class ct, template <class,int> class QuadratureTraits >
  class QuadratureProvider<ct, 3, QuadratureTraits> 
  {
    typedef QuadratureTraits<ct,3> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: SimplexQuadratureType SimplexQuadratureType;
    typedef typename QuadratureTraitsType :: CubeQuadratureType    CubeQuadratureType;
    typedef typename QuadratureTraitsType :: PrismQuadratureType   PrismQuadratureType;
    typedef typename QuadratureTraitsType :: PyramidQuadratureType PyramidQuadratureType;
  public:
    //! Access to the quadrature implementation.
    static const QuadratureImp<ct, 3>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert( geo.isTetrahedron() || geo.isHexahedron() || 
              geo.isPrism() || geo.isPyramid() );
      assert(order >= 0);

      if(geo.isTetrahedron()) 
        return QuadCreator::template provideQuad<SimplexQuadratureType> (order);
      if(geo.isHexahedron())
        return QuadCreator::template provideQuad<CubeQuadratureType> (order);
      if(geo.isPrism())
        return QuadCreator::template provideQuad<PrismQuadratureType> (order);
      if(geo.isPyramid())
        return QuadCreator::template provideQuad<PyramidQuadratureType> (order);

      DUNE_THROW(RangeError, "Element type not available for dim == 3");
      // dummy return
      return QuadCreator::template provideQuad<CubeQuadrature<ct, 3> > (0);
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  };

} // end namespace Dune 
#endif
