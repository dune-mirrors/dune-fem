#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

//- system includes 
#include <map>

#include <dune/fem/misc/utility.hh>

#include "pass.hh"
#include "selection.hh"
#include "discretemodel.hh"
#include "modelcaller.hh"

#include <dune/fem/solver/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/fem/space/common/allgeomtypes.hh> 
#include <dune/fem/space/common/arrays.hh> 
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/selection.hh>

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

    //! map that stores the volume of the reference element  
    typedef std::map<const Dune::GeometryType, double> RefVolumeMapType;

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
    typedef DiscreteModelCaller< DiscreteModelType 
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
      refVolMap_ (),
      updEn_(spc_),
      updNeigh_(spc_),
      dtMin_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      grads_(0.0),
      minLimit_(2.0*std::numeric_limits<double>::min()),
      diffVar_(),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? 
          (2*spc_.order()) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? 
        (2*spc_.order()+1) : faceQuadOrd ),
      localMassMatrix_( spc_ , volumeQuadOrd_ ) 
    {
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

    void applyLocal(EntityType& en) const 
    {
      // init local function 
      initLocalFunction( en , updEn_ );

      // call real apply local 
      applyLocal(en, updEn_);
      
      // add update to real function 
      updateFunctionAndApplyMass(en, updEn_ );
    }

    //! local integration 
    void applyLocal(EntityType& en, TemporaryLocalFunctionType& updEn) const
    {
      //- statements
      caller_.setEntity(en);

      // only call geometry once, who know what is done in this function 
      const Geometry & geo = en.geometry();

      // get volume of element 
      const double vol = geo.volume(); 

      // only apply volumetric integral if order > 0 
      // otherwise this contribution is zero 
      
      if( (spc_.order() > 0) || problem_.hasSource()) 
      {
        // if only flux, evaluate only flux 
        if ( problem_.hasFlux() && !problem_.hasSource() ) 
        {
          evalVolumetricPartFlux(en, geo, updEn);
        }
        else 
        {
          // evaluate flux and source 
          evalVolumetricPartBoth(en, geo, updEn);
        }
      }

      /////////////////////////////
      // Surface integral part
      /////////////////////////////
      if ( problem_.hasFlux() ) 
      {
        IntersectionIteratorType endnit = gridPart_.iend(en);
        for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
        {
          const IntersectionType& inter=*nit;
          //double nbvol;
          double nbvol = vol;
          double wspeedS = 0.0;
          if (inter.neighbor()) 
          {
            // get neighbor 
            EntityPointerType ep = inter.outside();
            EntityType & nb = *ep;
      
            if ( ! visited_[ indexSet_.index( nb ) ] ) 
            {
              // init local function 
              initLocalFunction( nb, updNeigh_ );

              typedef TwistUtility<GridType> TwistUtilityType;
              // for conforming situations apply Quadrature given
              if( TwistUtilityType::conforming(gridPart_.grid(),inter) )
              {
                FaceQuadratureType faceQuadInner(gridPart_, inter, faceQuadOrd_,
                                                 FaceQuadratureType::INSIDE);
          
                FaceQuadratureType faceQuadOuter(gridPart_, inter, faceQuadOrd_,
                                                 FaceQuadratureType::OUTSIDE);
                // apply neighbor part, return is volume of neighbor which is
                // needed below 
                nbvol = applyLocalNeighbor(nit,en,nb,
                          faceQuadInner,faceQuadOuter,
                          updEn, updNeigh_ , 
                          wspeedS);
              }
              else
              { 
                // for non-conforming situations apply the non-conforming 
                // type of the qaudrature 
                 
                // we only should get here whne a non-conforming situation 
                // occurs in a non-conforming grid 
                assert( GridPartType :: conforming == false );

                typedef typename FaceQuadratureType :: NonConformingQuadratureType
                  NonConformingFaceQuadratureType;

                NonConformingFaceQuadratureType 
                    nonConformingFaceQuadInner(gridPart_, inter, faceQuadOrd_,
                                               NonConformingFaceQuadratureType::INSIDE);

                NonConformingFaceQuadratureType 
                    nonConformingFaceQuadOuter(gridPart_,inter, faceQuadOrd_,
                                               NonConformingFaceQuadratureType::OUTSIDE);

                // apply neighbor part, return is volume of neighbor which is
                // needed below 
                nbvol = applyLocalNeighbor(nit,en,nb,
                          nonConformingFaceQuadInner,
                          nonConformingFaceQuadOuter,
                          updEn, updNeigh_ , 
                          wspeedS);
              }

              // add update to real function 
              updateFunction(nb, updNeigh_ );

            } // end if do something 
              
          } // end if neighbor
          else if( inter.boundary() )
          {
            FaceQuadratureType faceQuadInner(gridPart_, inter, faceQuadOrd_, 
                                             FaceQuadratureType::INSIDE);
            //nbvol = vol;
            caller_.setNeighbor(en);
            const int faceQuadInner_nop = faceQuadInner.nop();
            for (int l = 0; l < faceQuadInner_nop; ++l) 
            {
              // eval boundary Flux  
              wspeedS += caller_.boundaryFlux(nit, faceQuadInner, l, valEn_ )
                       * faceQuadInner.weight(l);
              
              // apply weights 
              valEn_ *= -faceQuadInner.weight(l);

              // add factor 
              updEn.axpy( faceQuadInner[l], valEn_  );
            }
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
    void initLocalFunction(EntityType& en, LocalFunctionImp& update) const 
    {
      // init local function  
      update.init( en );
      // clear dof values 
      update.clear();
    }

    //! add update to destination 
    template <class LocalFunctionImp>
    void updateFunction(EntityType& en, 
                        LocalFunctionImp& update) const 
    {
      // get local function and add update 
      LocalFunctionType function = dest_->localFunction(en);
      function += update;
    }

    //! add update to destination 
    template <class LocalFunctionImp>
    void updateFunctionAndApplyMass(
                        EntityType& en, 
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
    void evalVolumetricPartFlux(EntityType& en, 
                                const Geometry& geo , 
                                LocalFunctionImp& updEn) const
    {
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        // evaluate analytical flux and source 
        caller_.analyticalFlux(en, volQuad, l, fMat_ );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);
        
        // apply integration weights 
        fMat_ *= intel;

        // add values to local function 
        updEn.axpy( volQuad[l], fMat_);
      }
    }
    
    //////////////////////////////////////////
    // Volumetric integral part only flux 
    //////////////////////////////////////////
    template <class LocalFunctionImp>
    void evalVolumetricPartBoth(EntityType& en, 
                                const Geometry& geo, 
                               LocalFunctionImp& updEn) const
    {
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        // evaluate analytical flux and source 
        caller_.analyticalFluxAndSource(en, volQuad, l, fMat_, source_ );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);
        
        // apply integration weights 
        source_ *= intel;
        fMat_   *= intel;
        
        // add values to local function 
        updEn.axpy(volQuad[l],source_,fMat_);
      }
    }
    
    template <class QuadratureImp, class LocalFunctionImp >  
    double applyLocalNeighbor(IntersectionIteratorType & nit, 
            EntityType & en, EntityType & nb, 
            const QuadratureImp & faceQuadInner, 
            const QuadratureImp & faceQuadOuter,
            LocalFunctionImp & updEn,
            LocalFunctionImp & updNeigh,
            double & wspeedS) const 
    {
      // make Entity known in caller  
      caller_.setNeighbor(nb);
      
      // get goemetry of neighbor 
      const Geometry & nbGeo = nb.geometry();

      // get neighbors volume 
      const double nbvol = nbGeo.volume(); 
      
      const int faceQuadInner_nop = faceQuadInner.nop();
      for (int l = 0; l < faceQuadInner_nop; ++l) 
      {
        wspeedS += caller_.numericalFlux(nit, 
                                         faceQuadInner, 
                                         faceQuadOuter,
                                         l, 
                                         valEn_, valNeigh_)
                 * faceQuadInner.weight(l);
        
        // apply weights 
        valEn_    *= -faceQuadInner.weight(l);
        valNeigh_ *=  faceQuadOuter.weight(l);
        
        // add values to local functions 
        updEn.axpy   ( faceQuadInner[l], valEn_ );
        updNeigh.axpy( faceQuadOuter[l], valNeigh_ );
      }
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

    mutable RefVolumeMapType refVolMap_;

    mutable TemporaryLocalFunctionType updEn_;
    mutable TemporaryLocalFunctionType updNeigh_;

    mutable double dtMin_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_,fMatTmp;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable DomainType grads_;
    const double minLimit_;
    FieldVector<int, 0> diffVar_;

    const int volumeQuadOrd_, faceQuadOrd_;
    LocalMassMatrixType localMassMatrix_;
  };
//! @}  
} // end namespace Dune

#endif
