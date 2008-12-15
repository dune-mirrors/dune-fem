/**************************************************************
 **** Structure of numerical class:
 **************************************************************
 **** Traits class: parameters are dimRange and dimRange1
   DomainType: vector in world-space FV<dimDomain>
   RangeType: vector in phase-space FV<dimRange>
   FluxRangeType: matrix for analytical flux FM<Range,Domain>
   GradientType: vector in gradient phase-space FV<dimRange1>
   DiffusionRangeType: matrix for diffusion flux FM<Range1,Domain>
 **** Description class for    V + div a(U)          = 0
                            dt U + div (F(U)+A(U,V)) = 0

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
