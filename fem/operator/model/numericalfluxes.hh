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
  }; // end of LLFFlux class


  // *** the Local Lax-Friedrichs flux, modified for the wetting-drying treatment *** //
  template <class Model>
  class WDLLFFlux
  {
    public:
      typedef typename Model::Traits Traits;
      enum { dimRange = Model::dimRange };
      typedef typename Model::RangeType RangeType;
      typedef typename Model::FluxRangeType FluxRangeType;
    public:
      WDLLFFlux(Model& mod) 
        : model_(mod) {}
          const Model& model() const
          {
            return model_;
          }
          Model& model() 
          {
            return model_;
          }
    
    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    inline double numericalFlux(typename Traits::IntersectionIterator& it,
                                double time, 
                                const typename Traits::FaceDomainType& x,
                                const RangeType& uLeft, 
                                const RangeType& uRight,
                                bool reflectionLeft, bool reflectionRight,
                                RangeType& gLeft,
                                RangeType& gRight) const
    {
      const typename Traits::DomainType normal = it->integrationOuterNormal(x);  
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;

      gLeft = 0;
      gRight = 0;

      if (reflectionRight == false) // (no reflection) or (reflection in entity)
      {
        // anaflux = F( uLeft )
        std::cout << " reflectionRight = false, compute g(uLeft) .. "  << std::endl;
        model_.analyticalFlux( *(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               uLeft, anaflux);
        anaflux.umv( normal, gLeft );
      }
      else // (reflection in neighbor)
      {
        // anaflux = F( uRight )
        std::cout << " reflectionRight = true, compute g(uRight) .. "  << std::endl;
        model_.analyticalFlux( *(it->outside()), time,
                               it->intersectionNeighborLocal().global(x),
                               uRight, anaflux);
        anaflux.umv( normal, gRight );
      }

      if (reflectionLeft == true) // (reflection in entity)
      {
        // anaflux = F( uLeftRef )
        std::cout << " reflectionLeft = true, compute g(uLeftRef) .. "  << std::endl;
        RangeType uLeftRef;
        model_.analyticalFlux( *(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               model_.reflectU( uLeft, normal, uLeftRef ), anaflux);

        /*reflectU ( uLeft, normal, uLeftRef );
        model_.analyticalFlux( *(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               uLeftRef, anaflux);*/

        /*model_.analyticalFlux( *(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               reflectU( uLeft, normal ), anaflux); */

//                               uLeft, anaflux);
        anaflux.umv( normal, gLeft );
      }
      else
      {
        if (reflectionRight == true) // (reflection in neigbor)
        {
          // anaflux = F( uRightRef )
          std::cout << " reflectionRight = true, compute g(uRightRef) ..  "  << std::endl;
          RangeType uRightRef;
          model_.analyticalFlux( *(it->outside()), time,
                                 it->intersectionNeighborLocal().global(x),
                                 model_.reflectU( uRight, normal, uRightRef ), anaflux);
          anaflux.umv( normal, gRight );
        }
        else // (no reflection)
        {
          std::cout << " reflectionLeft = false, reflectionRight = false, compute g(uRight) ..  "  << std::endl;
          if (it->neighbor())
             model_.analyticalFlux( *(it->outside()), time,
                                    it->intersectionNeighborLocal().global(x),
                                    uRight, anaflux );
          else
          {
             model_.analyticalFlux( *(it->inside()), time,
                                    it->intersectionSelfLocal().global(x),
                                    uRight, anaflux );
          }
          anaflux.umv( normal,gLeft );
        }
      }

/*      model_.analyticalFlux( *(it->inside()), time,
                             it->intersectionSelfLocal().global(x),
                             uLeft, anaflux);
      gLeft = 0;
      anaflux.umv( normal, gLeft );
      if (it->neighbor())
         model_.analyticalFlux( *(it->outside()), time,
                                it->intersectionNeighborLocal().global(x),
                                uRight, anaflux );
      else
         model_.analyticalFlux( *(it->inside()), time,
                                it->intersectionSelfLocal().global(x),
                                uRight, anaflux );
      anaflux.umv( normal,gLeft ); */

      double maxspeedl, maxspeedr, maxspeed;
      double viscparal, viscparar, viscpara;

// nsh: the version of the method maxSpeed is used with the additional argument for wetting-drying treatment
      model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->inside()),
                       uLeft, viscparal, maxspeedl);
      model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->outside()),
                       uRight, viscparar, maxspeedr);

      maxspeed = (maxspeedl>maxspeedr) ? maxspeedl : maxspeedr;
      viscpara = (viscparal>viscparar) ? viscparal : viscparar;

      visc  = uRight;
      visc -= uLeft;
      visc *= 2.*viscpara;
      gLeft -= visc;
      
      gLeft *= 0.5;

      gRight = gLeft;

      return maxspeed;
    }
  private:
    Model& model_;
  }; // end of WDLLFFlux class

} // end namespace Dune

#endif

