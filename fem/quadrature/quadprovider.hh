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
        if(vec_.size() <= order) 
        {
          std::cerr << "WARNING: couldn't create quadrature of order=" << order << " in: " __FILE__ << " line: "<< __LINE__ <<"\n";
          order = vec_.size()-1;
        }
        
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
  template <typename ct, int dim>
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
  template <typename ct>
  class QuadratureProvider<ct, 0> 
  {
  public:
    //! Access to the quadrature implementations.
    static const QuadratureImp<ct, 0>& getQuadrature(GeometryType geo, 
                                                     int order) 
    {
      return QuadCreator::template 
        provideQuad<CubeQuadrature<ct, 0> >(order);
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  }; 

  
  //! Specialisation for dimension 1.
  template <typename ct>
  class QuadratureProvider<ct, 1> 
  {
  public:
    //! Access to the quadrature implementations.
    static const QuadratureImp<ct, 1>& getQuadrature(GeometryType geo, 
                                                     int order) {
      assert(geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      return QuadCreator::template provideQuad<CubeQuadrature<ct, 1> > (order);
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  }; 

  //! Specialisation for dimension = 2
  template <typename ct>
  class QuadratureProvider<ct, 2> 
  {
  public:
    //! Access to the quadrature implemenations.
    static const QuadratureImp<ct, 2>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert( geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      if(geo.isTriangle()) 
        return QuadCreator::template provideQuad<SimplexQuadrature<ct, 2> >(order);
      if(geo.isQuadrilateral())
        return QuadCreator::template provideQuad<CubeQuadrature<ct, 2> >(order);

      DUNE_THROW(RangeError, "Element type not available for dim == 2");
      // dummy return
      return QuadCreator::template provideQuad<SimplexQuadrature<ct, 2> >(0);
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  };

  //! Specialisation for dimension = 3.
  template <class ct>
  class QuadratureProvider<ct, 3> 
  {
  public:
    //! Access to the quadrature implementation.
    static const QuadratureImp<ct, 3>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert( geo.isTetrahedron() || geo.isHexahedron() || 
              geo.isPrism() || geo.isPyramid() );
      assert(order >= 0);

      if(geo.isTetrahedron()) 
        return QuadCreator::template provideQuad<SimplexQuadrature<ct, 3> > (order);
      if(geo.isHexahedron())
        return QuadCreator::template provideQuad<CubeQuadrature<ct, 3> > (order);
      if(geo.isPrism())
        return QuadCreator::template provideQuad<PrismQuadrature<ct> > (order);
      if(geo.isPyramid())
        return QuadCreator::template provideQuad<PyramidQuadrature<ct> > (order);

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
