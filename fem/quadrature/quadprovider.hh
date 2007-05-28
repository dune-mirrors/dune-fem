#ifndef DUNE_QUADPROVIDER_HH
#define DUNE_QUADPROVIDER_HH

#include "quadrature.hh"
#include "idprovider.hh"

namespace Dune {

  template <class ct, int dim>
  class QuadratureImp;
  
  //! Utility class that factors out the repetitive lookup process.
  //! the template parameter is only to distinguish the classes for
  //! different geometry types if the QuadImp are the same for different
  //! geometry classes 
  template <int dummy> 
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
      QuadratureStorage() : vec_(QuadImp::maxOrder()+1 ,(QuadImp*) 0) {}
      
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
      QuadImp& getQuadrature(const GeometryType& geo, size_t order) 
      {
        assert( order >= 0 );
        assert( order < vec_.size() );
        // if not exists, create quadrature  
        if(!vec_[order]) 
        {
          vec_[order] = new QuadImp(geo, order, IdProvider::instance().newId());
        }
        return *(vec_[order]);    
      }      
    }; // end class QuadratureStorage 
  
  public:
    //! provide quadrature, stores internal singleton holding list 
    //! with existing quadratures
    template <class QuadImp>
    static const QuadImp& provideQuad(const GeometryType& geo, int order)
    {
      // vector holding pointers 
      static QuadratureStorage<QuadImp> storage;
      return storage.getQuadrature(geo,order);
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
    typedef QuadratureTraits<ct,dim> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: CubeQuadratureType CubeQuadratureType;
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;

    //! Access to the quadrature implementations.
    static const IntegrationPointListType& getQuadrature(const GeometryType& geo,
                                                         int order)
    {
      assert(geo.isCube());
      return QuadCreator<0>::template provideQuad<CubeQuadratureType>(geo,order);
    } 
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
    
  };

  //! Specialisation for dimension 0
  template <typename ct, template <class,int> class QuadratureTraits >
  class QuadratureProvider<ct, 0, QuadratureTraits> 
  {
    typedef QuadratureTraits<ct,0> QuadratureTraitsType;
    typedef typename QuadratureTraitsType :: PointQuadratureType PointQuadratureType;
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;
  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType& getQuadrature(const GeometryType& geo,
                                                         int order) 
    {
      return QuadCreator<0>::template provideQuad<PointQuadratureType>(geo,order);
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
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;
  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType& getQuadrature(const GeometryType& geo,
                                                         int order) {
      assert(geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      return QuadCreator<0>::template provideQuad<LineQuadratureType> (geo,order);
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
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;
  public:
    //! Access to the quadrature implemenations.
    static const IntegrationPointListType& getQuadrature(const GeometryType& geo,
                                                         int order) {
      assert( geo.isCube() || geo.isSimplex() );
      assert(order >= 0);

      if(geo.isQuadrilateral())
        return QuadCreator<0>::template provideQuad<CubeQuadratureType>(geo,order);
      if(geo.isTriangle()) 
        return QuadCreator<1>::template provideQuad<SimplexQuadratureType>(geo,order);

      DUNE_THROW(RangeError, "Element type not available for dim == 2");
      // dummy return
      return QuadCreator<1>::template provideQuad<SimplexQuadratureType>(geo,0);
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
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;
  public:
    //! Access to the quadrature implementation.
    static const IntegrationPointListType& getQuadrature(const GeometryType& geo,
                                                         int order) {
      assert( geo.isTetrahedron() || geo.isHexahedron() || 
              geo.isPrism() || geo.isPyramid() );
      assert(order >= 0);

      if(geo.isHexahedron())
        return QuadCreator<0>::template provideQuad<CubeQuadratureType> (geo,order);
      if(geo.isTetrahedron()) 
        return QuadCreator<1>::template provideQuad<SimplexQuadratureType> (geo,order);
      if(geo.isPrism())
        return QuadCreator<2>::template provideQuad<PrismQuadratureType> (geo,order);
      if(geo.isPyramid())
        return QuadCreator<3>::template provideQuad<PyramidQuadratureType> (geo,order);

      DUNE_THROW(RangeError, "Element type not available for dim == 3");
      // dummy return
      return QuadCreator<1>::template provideQuad<SimplexQuadratureType> (geo,0);
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);
  };

} // end namespace Dune 
#endif
