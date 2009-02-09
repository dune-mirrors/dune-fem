#ifndef DUNE_FEM_LDG_NUMERICALFLUX_HH
#define DUNE_FEM_LDG_NUMERICALFLUX_HH

// *** Numerical fluxes ***

// Dune includes
namespace Dune {

  // *** the Local Lax-Friedrichs flux *** //
  template <class Model>
  class LLFFlux {
  public:
    typedef typename Model::Traits Traits;
    enum { dimRange = Model::dimRange };
    typedef typename Model::RangeType RangeType;
    typedef typename Model::FluxRangeType FluxRangeType;
  public:
    LLFFlux(const Model& mod) : model_(mod) {}
    const Model& model() const {return model_;}
    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    inline double numericalFlux(typename Traits::IntersectionIterator& it,
               double time, 
               const typename Traits::FaceDomainType& x,
               const RangeType& uLeft, 
               const RangeType& uRight,
               RangeType& gLeft,
               RangeType& gRight) const {
      const typename Traits::DomainType normal = it->integrationOuterNormal(x);  
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;

      model_.analyticalFlux(*(it->inside()), time,
                            it->intersectionSelfLocal().global(x),
                            uLeft, anaflux);
      gLeft*=0;
      anaflux.umv(normal, gLeft);
      if (it->neighbor())
         model_.analyticalFlux(*(it->outside()), time,
                               it->intersectionNeighborLocal().global(x),
                               uRight, anaflux);
      else
         model_.analyticalFlux(*(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               uRight, anaflux);
                               anaflux.umv(normal,gLeft);

      double maxspeedl,maxspeedr,maxspeed;
      double viscparal,viscparar,viscpara;

      model_.maxSpeed(normal,time,it->intersectionGlobal().global(x),
                      uLeft,viscparal,maxspeedl);
      model_.maxSpeed(normal,time,it->intersectionGlobal().global(x),
                      uRight,viscparar,maxspeedr);

      maxspeed=(maxspeedl>maxspeedr)?maxspeedl:maxspeedr;
      viscpara=(viscparal>viscparar)?viscparal:viscparar;

      visc=uRight;
      visc-=uLeft;
      visc*=2.*viscpara;
      gLeft-=visc;
      
      gLeft*=0.5;
      gRight=gLeft;

      return maxspeed;
    }
  private:
    const Model& model_;
  };

} // end namespace Dune

#endif

