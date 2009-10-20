/**************************************************************
 **** Structure of model class:
 **************************************************************
 **** Traits class: parameters are dimRange and dimRange1
   DomainType: vector in world-space FV<dimDomain>
   RangeType: vector in phase-space FV<dimRange>
   FluxRangeType: matrix for analytical flux FM<Range,Domain>
   GradientType: vector in gradient phase-space FV<dimRange1>
   DiffusionRangeType: matrix for diffusion flux FM<Range1,Domain>
 **** Description class for    V + div a(U)          = 0
                            dt U + div (F(U)+A(U,V)) = 0
      template <class GridType,class ProblemType>
      class ModelInterface {
      public:
        enum { dimDomain = GridType::dimensionworld };
        enum { dimRange = 1};
        typedef ModelTraits<GridType,dimRange> Traits;
        typedef typename Traits::RangeType RangeType;
        typedef typename Traits::DomainType DomainType;
        typedef typename Traits::FluxRangeType FluxRangeType;
        typedef typename Traits::GradientType GradientType;
        typedef typename Traits::DiffusionRangeType DiffusionRangeType;
      public:
        ModelInterface(Problem& prob) {}  // set parameters from Problem class
        // F
        inline void analyticalFlux(typename Traits::EntityType& en,
           double time,  
           const typename Traits::DomainType& x,
           const RangeType& u, 
           FluxRangeType& f) const {}
        // a
        inline void diffusion(typename Traits::EntityType& en,
                  double time, 
            const DomainType& x,
            const RangeType& u, 
            DiffusionRangeType& a) const {}
        // A
        inline double diffusion(typename Traits::EntityType& en,
              double time, 
              const DomainType& x,
              const RangeType& u, 
              const GradientType& v,
              FluxRangeType& A) const {}
        // **** Method for boundary treatment:
        inline bool hasBoundaryValue(typename Traits::IntersectionIterator& it,
                double time, 
                const typename Traits::FaceDomainType& x) const {
        // return true:  numerical flux is used with uRight given by boundaryValue
        // return false: boundaryFlux is used for flux computation
        inline void boundaryValue(typename Traits::IntersectionIterator& it,
          double time, 
                const typename Traits::FaceDomainType& x,
                const RangeType& uLeft, 
                RangeType& uRight) const {}
        // in the case that no boundary data is available:
        // return total flux on the boundary, i.e., advecton+diffusion part
        // Version 1: for advection-diffusion operator
        inline double boundaryFlux(typename Traits::IntersectionIterator& it,
             double time, 
                 const typename Traits::FaceDomainType& x,
                 const RangeType& uLeft, 
                 const GradientType& vLeft, 
                                   RangeType& gLeft) {} 
        // Version 2: for advection-diffusion operator
        inline double boundaryFlux(typename Traits::IntersectionIterator& it,
             double time, 
                 const typename Traits::FaceDomainType& x,
                 const RangeType& uLeft, 
                                   RangeType& gLeft) {} 
        // for time-step control: advspeed for advection-timestep, 
        //                        totalspeed for timestep with diffusion
        inline void maxSpeed(const typename Traits::DomainType& normal,
                 double time,  
                 const typename Traits::DomainType& x,
                 const RangeType& u,
                             double& advspeed,double& totalpeed) const {}
      };
      template <class Model>
      class FluxInterface {
      public:
        typedef typename Model::Traits Traits;
        enum { dimRange = Model::dimRange };
        typedef typename Model::RangeType RangeType;
        typedef typename Model::FluxRangeType FluxRangeType;
        typedef typename Model::DiffusionRangeType DiffusionRangeType;
      public:
        FluxInterface(const Model& mod) : model_(mod) {}
        const Model& model() const {return model_;}
        // Return value: maximum wavespeed*length of integrationOuterNormal
        // gLeft,gRight are fluxed * length of integrationOuterNormal
        inline double numericalFlux(typename Traits::IntersectionIterator& it,
            double time, 
            const typename Traits::FaceDomainType& x,
            const RangeType& uLeft, 
            const RangeType& uRight,
            RangeType& gLeft,
            RangeType& gRight) const {}
      private:
        const Model& model_;
      };

**************************************************************/
