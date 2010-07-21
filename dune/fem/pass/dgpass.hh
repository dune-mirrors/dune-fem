#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

//- system includes 
#include <dune/fem/misc/utility.hh>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
// #include "discretemodel.hh"
#include <dune/fem/pass/dgmodelcaller.hh>

#include <dune/fem/solver/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/fem/space/common/allgeomtypes.hh> 
#include <dune/fem/space/common/arrays.hh> 
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/selection.hh>

#include <dune/fem/quadrature/intersectionquadrature.hh>

namespace Dune {
/*! @addtogroup PassHyp
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   v + div(f(x,u)) + A(x,u)\nabla u &=& S(x,u)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ u \f$ is the argument and \f$ v \f$ is computed.
** Weak formulation on a cell T: 
** \f[
** \int_T v \phi = -\int_{\partial T} g \phi + \int_T f \cdot \nabla \phi + \int_T Q \phi 
** \f]
** with \f$ g \approx f \cdot n + \tilde{A}[u] \cdot n \f$ and \f$ Q \approx S - A \nabla u \f$
** where \f$ \tilde{A} \f$ denotes the arithmetic average and \f$ [u] \f$ the jump of 
** \f$ u \f$ over the cell interface.\\
** The discrete model provides the \b analyticalFlux f, the \b source Q and the \b numericalFlux g.
** @{
**************************************************************************/

  //! Concrete implementation of Pass for first hyperbolic systems using
  //! LDG
  template <class DiscreteModelImp, class PreviousPassImp , int passIdImp = -1 >
  class LocalDGPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp > 
  {
    typedef LocalDGPass< DiscreteModelImp , PreviousPassImp , passIdImp > ThisType;
  public:
    
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp > BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType:: BaseFunctionSetType
      BaseFunctionSetType; 

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef DGDiscreteModelCaller< DiscreteModelType 
                                   , ArgumentType 
                                   , CombinedSelectorType
                                 > DiscreteModelCallerType;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };

    // type of local id set 
    typedef typename GridPartType::IndexSetType IndexSetType; 
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType > TemporaryLocalFunctionType;

    //! type of local mass matrix 
    typedef LocalDGMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
    int volumeQuadOrd =-1,int faceQuadOrd=-1) :
      BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      indexSet_(gridPart_.indexSet()),
      visited_(0),
      updEn_(spc_),
      updNb_(spc_),
      fMatVec_( 20 ),
      valEnVec_( 20 ),
      valNbVec_( 20 ),
      dtMin_(std::numeric_limits<double>::max()),
      minLimit_(2.0*std::numeric_limits<double>::min()),
      //volumeQuadOrd_( 2 * spc_.order() ),
      //faceQuadOrd_( 2 * spc_.order() + 1),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? 
          ( 2 * spc_.order()) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? 
        ( 2 * spc_.order()+1) : faceQuadOrd ),
      localMassMatrix_( spc_ , volumeQuadOrd_ ) 
    {
      fMatVec_.setMemoryFactor( 1.1 );
      valEnVec_.setMemoryFactor( 1.1 );
      valNbVec_.setMemoryFactor( 1.1 );

      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );
    }
   
    //! Destructor
    virtual ~LocalDGPass() {
    }

    //! print tex info
    void printTexInfo(std::ostream& out) const {
      BaseType::printTexInfo(out);
      out << "LocalDGPass: "
          << "\\\\ \n";
    }
    
    //! Estimate for the timestep size 
    double timeStepEstimateImpl() const 
    {
      // factor for LDG  Discretization 
      const double p = 2 * spc_.order() + 1;
      return dtMin_ / p;
    }

  protected:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      // clear destination 
      dest_->clear();

      // set arguments to caller 
      caller_.setArgument(*arg_);

      // resize indicator function 
      visited_.resize( indexSet_.size(0) );
      // set all values to false 
      const int indSize = visited_.size();
      for(int i=0; i<indSize; ++i) visited_[i] = false;

      // time initialisation to max value 
      dtMin_ = std::numeric_limits<double>::max();

      // time is member of pass 
      caller_.setTime( this->time() );
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // communicate calculated function 
      spc_.communicate( dest );
      
      // call finalize 
      caller_.finalize();
    }

    void applyLocalMass( const EntityType& en) const 
    {
      // get local function and add update 
      LocalFunctionType function = dest_->localFunction(en);

      // apply local inverse mass matrix 
      localMassMatrix_.applyInverse( caller_, en, function );
    }

    void applyLocal( const EntityType& en) const 
    {
      // init local function 
      initLocalFunction( en , updEn_ );

      // call real apply local 
      applyLocal(en, updEn_);
      
      // add update to real function 
      updateFunctionAndApplyMass(en, updEn_ );
    }

    //! local integration 
    void applyLocal( const EntityType& en, TemporaryLocalFunctionType& updEn) const
    {
      // only call geometry once, who know what is done in this function 
      const Geometry & geo = en.geometry();

      // get volume of element 
      const double vol = geo.volume(); 

      // only apply volumetric integral if order > 0 
      // otherwise this contribution is zero 
      
      if( (spc_.order() > 0) || problem_.hasSource() ) 
      {
        // create quadrature 
        VolumeQuadratureType volQuad( en, volumeQuadOrd_ );

        // set entity and evaluate local functions for quadrature 
        caller_.setEntity( en , volQuad );

        if( problem_.hasSource() )
        {
          // evaluate flux and source
          evalVolumetricPartBoth( en, geo, volQuad, updEn );
        }
        else if ( problem_.hasFlux() ) 
        {
          // if only flux, evaluate only flux 
          evalVolumetricPartFlux( en, geo, volQuad, updEn );
        }
      }
      else 
      {
        // only set entity here without evaluation 
        caller_.setEntity( en );
      }

      /////////////////////////////
      // Surface integral part
      /////////////////////////////
      if ( problem_.hasFlux() ) 
      {
        IntersectionIteratorType endnit = gridPart_.iend(en);
        for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
        {
          // get intersection from intersection iterator
          const IntersectionType& intersection = *nit;

          //double nbvol;
          double nbvol = vol;
          double wspeedS = 0.0;
          if (intersection.neighbor()) 
          {
            // get neighbor 
            EntityPointerType ep = intersection.outside();
            EntityType & nb = *ep;
      
            if ( ! visited_[ indexSet_.index( nb ) ] ) 
            {
              // init local function 
              initLocalFunction( nb, updNb_ );

              // for conforming situations apply Quadrature given
              if( ! GridPartType :: conforming && ! intersection.conforming() )
              {
                // apply neighbor part, return is volume of neighbor which is
                // needed below 
                nbvol = applyLocalNeighbor< false > 
                                    (intersection,en,nb,
                                     updEn, updNb_ , 
                                     wspeedS);
              }
              else
              { 
                // apply neighbor part, return is volume of neighbor which is
                // needed below 
                nbvol = applyLocalNeighbor< true > 
                                    (intersection,en,nb,
                                     updEn, updNb_ , 
                                     wspeedS);
              }

              // add update to real function 
              updateFunction(nb, updNb_ );

            } // end if do something 
              
          } // end if neighbor
          else if( intersection.boundary() )
          {
            FaceQuadratureType faceQuadInner(gridPart_, intersection, faceQuadOrd_, 
                                             FaceQuadratureType::INSIDE);

            // set neighbor entity to inside entity 
            caller_.setBoundary(en, faceQuadInner);

            // cache number of quadrature points 
            const int faceQuadInner_nop = faceQuadInner.nop();

            if( valEnVec_.size() < faceQuadInner_nop )
              valEnVec_.resize( faceQuadInner_nop );

            // loop over quadrature points 
            for (int l = 0; l < faceQuadInner_nop; ++l) 
            {
              RangeType& flux = valEnVec_[ l ];

              // eval boundary Flux  
              wspeedS += caller_.boundaryFlux( *nit, faceQuadInner, l, flux )
                       * faceQuadInner.weight(l);
              
              // apply weights 
              flux *= -faceQuadInner.weight(l);
            }

            // add factor 
            updEn.axpyQuadrature ( faceQuadInner, valEnVec_ );

          } // end if boundary
          
          if (wspeedS > minLimit_ ) 
          {
            double minvolS = std::min(vol,nbvol);
            dtMin_ = std::min(dtMin_,minvolS/wspeedS);
          }
        } // end intersection loop

      } // end if problem_.hasFlux()

      // this entity is finised by now 
      visited_[ indexSet_.index( en ) ] = true ;
    }

    // initialize local update function 
    template <class LocalFunctionImp>
    void initLocalFunction( const EntityType& en, LocalFunctionImp& update) const 
    {
      // init local function  
      update.init( en );
      // clear dof values 
      update.clear();
    }

    //! add update to destination 
    template <class LocalFunctionImp>
    void updateFunction( const EntityType& en, 
                         LocalFunctionImp& update) const 
    {
      // get local function and add update 
      LocalFunctionType function = dest_->localFunction(en);
      function += update;
    }

    //! add update to destination 
    template <class LocalFunctionImp>
    void updateFunctionAndApplyMass(
                        const EntityType& en, 
                        LocalFunctionImp& update) const
    {
      // get local function and add update 
      LocalFunctionType function = dest_->localFunction(en);
      function += update;

      // apply local inverse mass matrix 
      localMassMatrix_.applyInverse( caller_, en, function );
    }

    //////////////////////////////////////////
    // Volumetric integral part only flux 
    //////////////////////////////////////////
    template <class LocalFunctionImp>
    void evalVolumetricPartFlux( const EntityType& en, 
                                 const Geometry& geo , 
                                 const VolumeQuadratureType& volQuad,
                                 LocalFunctionImp& updEn ) const
    {
      const int volQuad_nop = volQuad.nop();
      if( fMatVec_.size() < volQuad_nop ) 
      {
        fMatVec_.resize( volQuad_nop );
      }

      for (int l = 0; l < volQuad_nop; ++l) 
      {
        JacobianRangeType& flux = fMatVec_[ l ];

        // evaluate analytical flux and source 
        caller_.analyticalFlux(en, volQuad, l, flux );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);
        
        // apply integration weights 
        flux *= intel;
      }

      // add values to local function 
      updEn.axpyQuadrature (volQuad, fMatVec_ );
    }
    
    //////////////////////////////////////////
    // Volumetric integral part only flux 
    //////////////////////////////////////////
    template <class LocalFunctionImp>
    void evalVolumetricPartBoth( const EntityType& en, 
                                 const Geometry& geo, 
                                 const VolumeQuadratureType& volQuad,
                                 LocalFunctionImp& updEn ) const
    {
      const int volQuad_nop = volQuad.nop();

      if( fMatVec_.size() < volQuad_nop ) 
      {
        fMatVec_.resize( volQuad_nop );
      }

      if( valEnVec_.size() < volQuad_nop ) 
      {
        valEnVec_.resize( volQuad_nop );
      }

      for (int l = 0; l < volQuad_nop; ++l) 
      {
        JacobianRangeType& flux = fMatVec_[ l ];
        RangeType& source = valEnVec_[ l ];

        // evaluate analytical flux and source 
        const double dtEst =
          caller_.analyticalFluxAndSource(en, volQuad, l, flux, source );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);
        
        // apply integration weights 
        source *= intel;
        flux   *= intel;

        if( dtEst > minLimit_ )
          dtMin_ = std::min(dtMin_, dtEst);
      }

      // add values to local function 
      updEn.axpyQuadrature (volQuad, valEnVec_, fMatVec_ );
    }
    
    template <bool conforming, class LocalFunctionImp>  
    double applyLocalNeighbor ( const IntersectionType &intersection,
                                const EntityType &en, 
                                const EntityType &nb, 
                                LocalFunctionImp &updEn,
                                LocalFunctionImp &updNb,
                                double &wspeedS ) const 
    {
      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      // use IntersectionQuadrature to create appropriate face quadratures 
      typedef IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType; 
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

      // create intersection quadrature 
      IntersectionQuadratureType interQuad( gridPart_, intersection, faceQuadOrd_ );

      // get appropriate references 
      const QuadratureImp &faceQuadInner = interQuad.inside();
      const QuadratureImp &faceQuadOuter = interQuad.outside();

      // make Entity known in caller  
      caller_.setNeighbor(nb, faceQuadInner, faceQuadOuter);
     
      // get goemetry of neighbor 
      const Geometry & nbGeo = nb.geometry();

      // get neighbors volume 
      const double nbvol = nbGeo.volume(); 
      
      const int faceQuadInner_nop = faceQuadInner.nop();
      assert( faceQuadInner.nop() == faceQuadOuter.nop() );

      // check  valNbVec_ here, valEnVec_ might have been resized  
      if( valNbVec_.size() < faceQuadInner_nop )
      {
        valEnVec_.resize( faceQuadInner_nop );
        valNbVec_.resize( faceQuadInner_nop );
      }

      for (int l = 0; l < faceQuadInner_nop; ++l) 
      {
        RangeType& fluxEn = valEnVec_[ l ];
        RangeType& fluxNb = valNbVec_[ l ];

        wspeedS += caller_.numericalFlux(intersection, 
                                         faceQuadInner, 
                                         faceQuadOuter,
                                         l, 
                                         fluxEn, fluxNb )
                 * faceQuadInner.weight(l);
        
        // apply weights 
        fluxEn *= -faceQuadInner.weight(l);
        fluxNb *=  faceQuadOuter.weight(l);
      }
        
      // add values to local functions 
      updEn.axpyQuadrature( faceQuadInner, valEnVec_ );
      updNb.axpyQuadrature( faceQuadOuter, valNbVec_ );

      return nbvol;
    }
                         
  private:
    LocalDGPass();
    LocalDGPass(const LocalDGPass&);
    LocalDGPass& operator=(const LocalDGPass&);

  protected:
    mutable DiscreteModelCallerType caller_;
    const DiscreteModelType& problem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const IndexSetType& indexSet_;

    // indicator for grid walk 
    mutable MutableArray<bool> visited_;

    mutable TemporaryLocalFunctionType updEn_;
    mutable TemporaryLocalFunctionType updNb_;

    //! Some helper variables
    mutable MutableArray< JacobianRangeType > fMatVec_;
    mutable MutableArray< RangeType > valEnVec_;
    mutable MutableArray< RangeType > valNbVec_;

    mutable double dtMin_;
    const double minLimit_;

    const int volumeQuadOrd_, faceQuadOrd_;
    LocalMassMatrixType localMassMatrix_;
  };
//! @}  
} // end namespace Dune

#endif
