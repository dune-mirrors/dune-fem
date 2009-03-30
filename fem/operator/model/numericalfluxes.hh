#ifndef DUNE_FEM_LDG_NUMERICALFLUX_HH
#define DUNE_FEM_LDG_NUMERICALFLUX_HH

#include <dune/common/version.hh>

// *** Numerical fluxes ***

// Dune includes
namespace Dune
{

  // *** the Local Lax-Friedrichs flux *** //
  template< class Model >
  struct LLFFlux
  {
    typedef typename Model::Traits Traits;
    enum { dimRange = Model::dimRange };
    typedef typename Model::RangeType RangeType;
    typedef typename Model::FluxRangeType FluxRangeType;

    typedef typename Traits::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    LLFFlux ( const Model &model )
    : model_( model )
    {}

    const Model &model () const
    {
      return model_;
    }

    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    double numericalFlux ( const Intersection &intersection,
                           double time, 
                           const typename Traits::FaceDomainType &x,
                           const RangeType &uLeft, 
                           const RangeType &uRight,
                           RangeType &gLeft,
                           RangeType &gRight) const
    {
      typedef typename Intersection::Geometry IntersectionGeometry;
      typedef typename Intersection::LocalGeometry IntersectionLocalGeometry;

      const typename Traits::DomainType normal = intersection.integrationOuterNormal(x);
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;

#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
      const IntersectionLocalGeometry &geoInInside = intersection.geometryInInside();
#else
      const IntersectionLocalGeometry &geoInInside = intersection.intersectionSelfLocal();
#endif

      model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global( x ),
                             uLeft, anaflux);
      gLeft*=0;
      anaflux.umv(normal, gLeft);
      if( intersection.neighbor() )
      {
#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
        const IntersectionLocalGeometry &geoInOutside = intersection.geometryInOutside();
#else
        const IntersectionLocalGeometry &geoInOutside = intersection.intersectionNeighborLocal();
#endif
         model_.analyticalFlux( *(intersection.outside()), time, geoInOutside.global( x ),
                                uRight, anaflux);
      }
      else
      {
         model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global( x ),
                                uRight, anaflux);
      }
      anaflux.umv( normal, gLeft );

      double maxspeedl,maxspeedr,maxspeed;
      double viscparal,viscparar,viscpara;

#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
      const IntersectionGeometry &geo = intersection.geometry();
#else      
      const IntersectionGeometry &geo = intersection.intersectionGlobal();
#endif
      model_.maxSpeed( normal, time, geometry.global( x ), uLeft, viscparal, maxspeedl );
      model_.maxSpeed( normal, time, geometry.global( x ), uRight, viscparar, maxspeedr );

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

    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    double numericalFlux ( const IntersectionIterator &it,
                           double time, 
                           const typename Traits::FaceDomainType &x,
                           const RangeType &uLeft, 
                           const RangeType &uRight,
                           RangeType &gLeft,
                           RangeType &gRight) const
    {
      return numericalFlux( *it, time, x, uLeft, uRight, gLeft, gRight );
    }

  private:
    const Model &model_;
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

      //std::cout << " reflectionLeft is " << reflectionLeft << std::endl;
      //std::cout << " reflectionRight is " << reflectionRight << std::endl;

      gLeft = 0;
      gRight = 0;
      RangeType uLeftRef;
      RangeType uRightRef;

      if (reflectionRight == false) // (no reflection) or (reflection in entity)
      {
        // anaflux = F( uLeft )
        //std::cout << " reflectionRight = false, compute g(uLeft) .. "  << std::endl;
        model_.analyticalFlux( *(it->inside()), time,
                               it->intersectionSelfLocal().global(x),
                               uLeft, anaflux);
        anaflux.umv( normal, gLeft );
      }
      else // (reflection in neighbor)
      {
        // anaflux = F( uRight )
        //std::cout << " reflectionRight = true, compute g(uRight) .. "  << std::endl;
        model_.analyticalFlux( *(it->outside()), time,
                               it->intersectionNeighborLocal().global(x),
                               uRight, anaflux);
        anaflux.umv( normal, gRight );
      }

      if (reflectionLeft == true) // (reflection in entity)
      {
        // anaflux = F( uLeftRef )
        //std::cout << " reflectionLeft = true, compute g(uLeftRef) .. "  << std::endl;
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
          //std::cout << " reflectionRight = true, compute g(uRightRef) ..  "  << std::endl;
          model_.analyticalFlux( *(it->outside()), time,
                                 it->intersectionNeighborLocal().global(x),
                                 model_.reflectU( uRight, normal, uRightRef ), anaflux);
          anaflux.umv( normal, gRight );
        }
        else // (no reflection)
        {
          //std::cout << " reflectionLeft = false, reflectionRight = false, compute g(uRight) ..  "  << std::endl;
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

      double maxspeedl, maxspeedr, maxspeed;
      double viscparal, viscparar, viscpara;

// the version of the method maxSpeed is used with the additional argument for wetting-drying treatment

      if (reflectionRight == false) // (no reflection) or (reflection in entity)
      {
        model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->inside()),
                         uLeft, viscparal, maxspeedl);
      }
      else // (reflection in neighbor)
      {
        model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->outside()),
                         uRight, viscparar, maxspeedr);
      }

      if (reflectionLeft == true) // (reflection in entity)
      {
        model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->inside()),
                         uLeftRef, viscparar, maxspeedr);
      }
      else
      {
        if (reflectionRight == true) // (reflection in neigbor)
        {
          model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->outside()),
                           uRightRef, viscparal, maxspeedl );
        }
        else // (no reflection)
        {
          model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->outside()),
                           uRight, viscparar, maxspeedr );
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
      anaflux.umv( normal,gLeft );

      model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->inside()),
                       uLeft, viscparal, maxspeedl);

      model_.maxSpeed( normal, time, it->intersectionGlobal().global(x), *(it->outside()),
                       uRight, viscparar, maxspeedr); */

      maxspeed = (maxspeedl>maxspeedr) ? maxspeedl : maxspeedr;
      viscpara = (viscparal>viscparar) ? viscparal : viscparar;

      //std::cout << " reflectionRight = true, compute g(uRightRef) ..  "  << std::endl;

      if (reflectionRight == true)
      {
        visc  = uRightRef - uRight;  
      }
      if (reflectionLeft == true)
      {
        visc  = uLeftRef - uLeft;  
      }
      else
      {
        visc  = uRight - uLeft;
      }

      visc *= 2.*viscpara;

      if (reflectionRight == true)
      {
        gRight -= visc;
        gRight *= 0.5;        
        gLeft = gRight;  
      }
      else
      {
        gLeft -= visc;
        gLeft *= 0.5;
        gRight = gLeft;
      }

      return maxspeed;
    }
  private:
    Model& model_;
  }; // end of WDLLFFlux class

} // end namespace Dune

#endif

