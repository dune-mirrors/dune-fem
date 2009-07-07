#ifndef DUNE_FEM_LDG_NUMERICALFLUX_HH
#define DUNE_FEM_LDG_NUMERICALFLUX_HH

// *** Numerical fluxes ***

namespace Dune
{

  // local Lax-Friedrichs flux
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

    // return value: maximum wavespeed * length of integrationOuterNormal
    // gLeft, gRight are fluxes * length of integrationOuterNormal
    double numericalFlux ( const Intersection &intersection,
                           double time, 
                           const typename Traits::FaceDomainType &x,
                           const RangeType &uLeft, 
                           const RangeType &uRight,
                           RangeType &gLeft,
                           RangeType &gRight ) const
    {
      typedef typename Intersection::Geometry IntersectionGeometry;
      typedef typename Intersection::LocalGeometry IntersectionLocalGeometry;

      const typename Traits::DomainType normal = intersection.integrationOuterNormal(x);
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;

      const IntersectionLocalGeometry &geoInInside = intersection.geometryInInside();

      model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global( x ),
                             uLeft, anaflux);
      gLeft*=0;
      anaflux.umv(normal, gLeft);
      if( intersection.neighbor() )
      {
        const IntersectionLocalGeometry &geoInOutside = intersection.geometryInOutside();
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

      const IntersectionGeometry &geo = intersection.geometry();
      model_.maxSpeed( normal, time, geo.global( x ), uLeft, viscparal, maxspeedl );
      model_.maxSpeed( normal, time, geo.global( x ), uRight, viscparar, maxspeedr );

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

    // return value: maximum wavespeed * length of integrationOuterNormal
    // gLeft, gRight are fluxes * length of integrationOuterNormal
    double numericalFlux ( const IntersectionIterator &it,
                           double time, 
                           const typename Traits::FaceDomainType &x,
                           const RangeType &uLeft, 
                           const RangeType &uRight,
                           RangeType &gLeft,
                           RangeType &gRight ) const
    {
      return numericalFlux( *it, time, x, uLeft, uRight, gLeft, gRight );
    }

  private:
    const Model &model_;
  }; // end of LLFFlux


  // local Lax-Friedrichs flux, modified for the wetting-drying treatment for shallow water equations
  template <class Model>
  class WDLLFFlux
  {
    public:
      typedef typename Model::Traits Traits;
      enum { dimRange = Model::dimRange };
      typedef typename Model::RangeType RangeType;
      typedef typename Model::FluxRangeType FluxRangeType;

      typedef typename Traits::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;

    public:
      WDLLFFlux(Model& mod) 
        : model_(mod) {}
          const Model& model() const
          {
            return model_;
          }
          // in LLFFlux model_ is const, but here model contains wet-or-dry status flags, which are updated, therefore, const can't be used here
          Model& model()
          {
            return model_;
          }
    
    // return value: maximum wavespeed * length of integrationOuterNormal
    // gLeft, gRight are fluxes * length of integrationOuterNormal
    inline double numericalFlux( const Intersection &intersection,
                                 double time, 
                                 const typename Traits::FaceDomainType& x,
                                 const RangeType& uLeft, 
                                 const RangeType& uRight,
                                 bool reflectionLeft, bool reflectionRight,
                                 RangeType& gLeft,
                                 RangeType& gRight ) const
    {
      typedef typename Intersection::Geometry IntersectionGeometry;
      typedef typename Intersection::LocalGeometry IntersectionLocalGeometry;

      const typename Traits::DomainType normal = intersection.integrationOuterNormal(x);  
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;

#if WRITENUMFLUX == 1
      std::cout << std::endl;
      std::cout << " .. uLeft is " << uLeft << std::endl;
      std::cout << " .. uRight is " << uRight << std::endl;
      std::cout << " .. reflectionLeft is " << reflectionLeft << std::endl;
      std::cout << " .. reflectionRight is " << reflectionRight << std::endl;
#endif

      const IntersectionLocalGeometry &geoInInside = intersection.geometryInInside();
      const IntersectionLocalGeometry &geoInOutside = intersection.geometryInOutside();

      gLeft = 0;
      gRight = 0;
      RangeType uLeftRef;
      RangeType uRightRef;


      if ( reflectionLeft == false && reflectionRight == false )
      // (no reflection)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. no refelction: reflectionLeft = false & reflectionRight = false"  << std::endl;
        std::cout << " ..   evaluate f(a) = F( en, time, x_en, uLeft )"  << std::endl;
#endif
        // anaflux = F( en, time, x_en, uLeft )
        model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x),
                               uLeft, 
                               anaflux );
        // gLeft = F( en, time, x_en, uLeft ) * n
        anaflux.umv( normal, gLeft );

        if (intersection.neighbor())
        {
#if WRITENUMFLUX == 1
          std::cout << " ..   neighbor: evaluate f(b) = F( nb, time, x_nb, uRight)"  << std::endl;
#endif
          // anaflux = F( nb, time, x_nb, uRight)
          model_.analyticalFlux( *(intersection.outside()), time, geoInOutside.global(x),
                                 uRight, 
                                 anaflux );
        }
        else
        {
#if WRITENUMFLUX == 1
          std::cout << " ..   no neighbor: evaluate f(b) = F( en, time, x_en, uRight)"  << std::endl;
#endif
          // anaflux = F( en, time, x_en, uRight)
          model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x), 
                                 uRight, 
                                 anaflux );
        }
        
        // gLeft = F( en, time, x_en, uLeft ) * n + F( nb, time, x_nb, uRight) * n (neighbor)
        //   or
        // gLeft = F( en, time, x_en, uLeft ) * n + F( en, time, x_en, uRight) * n (no neighbor)
        anaflux.umv( normal,gLeft );

      } // end of (no reflection)


      if ( reflectionLeft == true )
      // (reflection in entity) or (reflection in both entity and neighbor)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. reflection in entity: reflectionLeft = true, reflectionRight can be true or false"  << std::endl;
        std::cout << " ..   evaluate f(a) = F( en, time, x_en, uLeft )"  << std::endl;
#endif
        // anaflux = F( en, time, x_en, uLeft )
        model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x),
                               uLeft, 
                               anaflux );
        // gLeft = F( en, time, x_en, uLeft ) * n
        anaflux.umv( normal, gLeft );
#if WRITENUMFLUX == 1
        std::cout << " ..   evaluate f(b) = F( en, time, x_en, uLeftRef )"  << std::endl;
#endif
        // anaflux = F( en, time, x_en, uLeftRef )
        model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x),
                               model_.reflectU( uLeft, normal, uLeftRef ), 
                               anaflux );

        // different notation ???
        //reflectU ( uLeft, normal, uLeftRef );
        //model_.analyticalFlux( *(intersection.inside()), time,
        //                       geoInInside.global(x),
        //                       uLeftRef, anaflux);

        // (compile with -O3) see /swe/model/swe.hh ???
        //model_.analyticalFlux( *(intersection.inside()), time,
        //                       geoInInside.global(x),
        //                       reflectU( uLeft, normal ), anaflux);

        // gLeft = F( en, time, x_en, uLeft ) * n + F( en, time, x_en, uLeftRef ) * n
        anaflux.umv( normal, gLeft );
      } // end of (reflection in entity) or (reflection in both entity and neighbor)
      else
      {
      if ( reflectionLeft == false && reflectionRight == true )
      // (reflection in neighbor)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. reflection in neighbor: reflectionLeft = false && reflectionRight = true "  << std::endl;
        std::cout << " ..   evaluate f(a) = F( nb, time, x_nb, uRightRef )"  << std::endl;
#endif
          // anaflux = F( nb, time, x_nb, uRightRef )
          model_.analyticalFlux( *(intersection.outside()), time,
                                 geoInOutside.global(x),
                                 model_.reflectU( uRight, normal, uRightRef ), 
                                 anaflux);
          // gLeft = F( nb, time, x_nb, uRightRef ) * n
          anaflux.umv( normal, gRight );

#if WRITENUMFLUX == 1
        std::cout << " ..   evaluate f(b) = F( nb, time, x_nb, uRight )"  << std::endl;
#endif
        // anaflux = F( nb, time, x_nb, uRight )
        model_.analyticalFlux( *(intersection.outside()), time, geoInOutside.global(x),
                               uRight, 
                               anaflux );
        // gLeft = F( nb, time, x_nb, uRightRef ) * n + F( nb, time, x_nb, uRight ) * n
        anaflux.umv( normal, gRight );
      } // end of (reflection in neighbor)
      } // end of else

      double maxspeedl, maxspeedr, maxspeed;
      double viscparal, viscparar, viscpara;

      const IntersectionGeometry &geo = intersection.geometry();

      // the version of the method maxSpeed is used with the additional argument for wetting-drying treatment
      if ( reflectionRight == false )
      // (no reflection) or (reflection in entity)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.inside()),
                         uLeft, viscparal, maxspeedl );
      }
      else
      // (reflection in neighbor)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                         uRight, viscparar, maxspeedr );
      }

      if ( reflectionLeft == true )
      // (reflection in entity)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.inside()),
                         uLeftRef, viscparar, maxspeedr );
      }
      else
      {
        if ( reflectionRight == true )
        // (reflection in neigbor)
        {
          model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                           uRightRef, viscparal, maxspeedl );
        }
        else
        // (no reflection)
        {
          model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                           uRight, viscparar, maxspeedr );
        }
      }

/*
      if ( reflectionRight == false )
      // (no reflection) or (reflection in entity)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. reflectionRight = false, evaluate F( en, time, x_en, uLeft )"  << std::endl;
#endif
        // anaflux = F( uLeft )
        model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x),
                               uLeft, 
                               anaflux );
        anaflux.umv( normal, gLeft );
      }
      else
      // (reflection in neighbor)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. reflectionRight = true, evaluate F( nb, time, x_nb, uRight )"  << std::endl;
#endif
        // anaflux = F( uRight )
        model_.analyticalFlux( *(intersection.outside()), time, geoInOutside.global(x),
                               uRight, 
                               anaflux );
        anaflux.umv( normal, gRight );
      }

      if ( reflectionLeft == true )
      // (reflection in entity)
      {
#if WRITENUMFLUX == 1
        std::cout << " .. reflectionLeft = true, evaluate F( en, time, x_en, uLeftRef )"  << std::endl;
#endif
        // anaflux = F( uLeftRef )
        model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x),
                               model_.reflectU( uLeft, normal, uLeftRef ), 
                               anaflux );

        //reflectU ( uLeft, normal, uLeftRef );
        //model_.analyticalFlux( *(intersection.inside()), time,
        //                       geoInInside.global(x),
        //                       uLeftRef, anaflux);
        //model_.analyticalFlux( *(intersection.inside()), time,
        //                       geoInInside.global(x),
        //                       reflectU( uLeft, normal ), anaflux);

        anaflux.umv( normal, gLeft );
      }
      else
      {
        if ( reflectionRight == true )
        // (reflection in neigbor)
        {
#if WRITENUMFLUX == 1
          std::cout << " .. reflectionRight = true, evaluate F ( nb, time, x_nb, uRightRef )"  << std::endl;
#endif
          // anaflux = F( uRightRef )
          model_.analyticalFlux( *(intersection.outside()), time,
                                 geoInOutside.global(x),
                                 model_.reflectU( uRight, normal, uRightRef ), 
                                 anaflux);
          anaflux.umv( normal, gRight );
        }
        else
        // (no reflection)
        {
#if WRITENUMFLUX == 1
      std::cout << " .. reflectionLeft = false, reflectionRight = false, evaluate F( nb, time, uRight)"  << std::endl;
#endif
          if (intersection.neighbor())
             model_.analyticalFlux( *(intersection.outside()), time, geoInOutside.global(x),
                                    uRight, 
                                    anaflux );
          else
          {
             model_.analyticalFlux( *(intersection.inside()), time, geoInInside.global(x), 
                                    uRight, 
                                    anaflux );
          }
          anaflux.umv( normal,gLeft );
        }
      }

      double maxspeedl, maxspeedr, maxspeed;
      double viscparal, viscparar, viscpara;

      const IntersectionGeometry &geo = intersection.geometry();

      // the version of the method maxSpeed is used with the additional argument for wetting-drying treatment
      if ( reflectionRight == false )
      // (no reflection) or (reflection in entity)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.inside()),
                         uLeft, viscparal, maxspeedl );
      }
      else
      // (reflection in neighbor)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                         uRight, viscparar, maxspeedr );
      }

      if ( reflectionLeft == true )
      // (reflection in entity)
      {
        model_.maxSpeed( normal, time, geo.global(x), *(intersection.inside()),
                         uLeftRef, viscparar, maxspeedr );
      }
      else
      {
        if ( reflectionRight == true )
        // (reflection in neigbor)
        {
          model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                           uRightRef, viscparal, maxspeedl );
        }
        else
        // (no reflection)
        {
          model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                           uRight, viscparar, maxspeedr );
        }
      }
*/

/*      model_.analyticalFlux( *(intersection.inside()), time,
                             geoInInside.global(x),
                             uLeft, anaflux );
      gLeft = 0;
      anaflux.umv( normal, gLeft );
      if (it->neighbor())
         model_.analyticalFlux( *(intersection.outside()), time,
                                geoInOutside.global(x),
                                uRight, anaflux );
      else
         model_.analyticalFlux( *(intersection.inside()), time,
                                geoInInside.global(x),
                                uRight, anaflux );
      anaflux.umv( normal,gLeft );

      model_.maxSpeed( normal, time, geo.global(x), *(intersection.inside()),
                       uLeft, viscparal, maxspeedl );

      model_.maxSpeed( normal, time, geo.global(x), *(intersection.outside()),
                       uRight, viscparar, maxspeedr ); */

      maxspeed = (maxspeedl>maxspeedr) ? maxspeedl : maxspeedr;
      viscpara = (viscparal>viscparar) ? viscparal : viscparar;

#if WRITENUMFLUX == 1
      std::cout << " .. reflectionRight = true, evaluate g(uRightRef)"  << std::endl;
#endif

      if ( reflectionRight == true )
      {
        visc  = uRightRef - uRight;  
      }
      if ( reflectionLeft == true )
      {
        visc  = uLeftRef - uLeft;  
      }
      else
      {
        visc  = uRight - uLeft;
      }

      visc *= 2.*viscpara;

      if ( reflectionRight == true )
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

    // return value: maximum wavespeed * length of integrationOuterNormal
    // gLeft, gRight are fluxes * length of integrationOuterNormal
    double numericalFlux ( const IntersectionIterator &it,
                           double time, 
                           const typename Traits::FaceDomainType &x,
                           const RangeType &uLeft, 
                           const RangeType &uRight,
                           bool reflectionLeft, bool reflectionRight,
                           RangeType &gLeft,
                           RangeType &gRight ) const
    {
      return numericalFlux( *it, time, x, uLeft, uRight, reflectionLeft, reflectionRight, gLeft, gRight );
    }

  private:
    Model& model_;
  }; // end of WDLLFFlux

} // end namespace Dune

#endif

