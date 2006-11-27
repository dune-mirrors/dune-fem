/**************************************************************
 * Known problems:
 * 1) works only with ortho-normal basefunction set
 *    we should have multiplication with inverse mass matrix on baseset
 * 2) Caching does not work with non-conforming grids
 *    here a switch is required
 * 3) would be good to pass quadrature to discrete model
 *    prehaps also a setEntity method on the discrete model
 * 4) hexaedrons with non-linear mapping will definitly not work
 *    (-> diplom thesis!)
*****************************************************************/


#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

#include "pass.hh"
#include "selection.hh"
#include "discretemodel.hh"
#include "modelcaller.hh"

// * needs to move
// #include "../misc/timenew.hh"
#include "../misc/timeutility.hh"

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/utility/twistutility.hh>

#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune {
/*! @defgroup PassHyp Local Discontinous Galerkin for first order hyperbolic equations
 *  @ingroup Pass
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   v + div(f(x,u)) + A(x,u)\nabla u &=& S(x,u)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ u \f$ is the argument and \f$ v \f$ is computed.
** @{
**************************************************************************/

  //! Concrete implementation of Pass for first hyperbolic systems using
  //! LDG
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

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
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // type of Communication Manager 
    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
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
                DiscreteFunctionSpaceType& spc,
    int volumeQuadOrd =-1,int faceQuadOrd=-1) :
      BaseType(pass, spc),
      caller_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      communicationManager_(spc_),
      dtMin_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      grads_(0.0),
      time_(0),
      diffVar_(),
      twistUtil_(spc.grid()),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? 
          (2*spc_.order()) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? 
        (2*spc_.order()+1) : faceQuadOrd )
    {
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );
    }
   
    //! Destructor
    virtual ~LocalDGPass() {
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void processTimeProvider(TimeProvider* time) {
      time_ = time;
    }

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return dtMin_;
    }

  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      dest_->clear();

      caller_.setArgument(*arg_);

      // time initialisation
      dtMin_ = std::numeric_limits<double>::max();
      if (time_) {
        caller_.setTime(time_->time());
      }
      else {
        caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // communicate calculated function 
      communicationManager_.exchange( dest );
      
      if (time_) {
        time_->provideTimeStepEstimate(dtMin_);
      }
      caller_.finalize();
    }

    void applyLocal(EntityType& en) const
    {
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;

      const IndexSetType& iset = spc_.indexSet();
      
      //- statements
      caller_.setEntity(en);
      LocalFunctionType updEn = dest_->localFunction(en);
      const int updEn_numDofs = updEn.numDofs();
      const BaseFunctionSetType& bsetEn = updEn.baseFunctionSet(); 
      

      // only call geometry once, who know what is done in this function 
      const GeometryType & geo = en.geometry();

      double massVolElinv;
      const double vol = volumeElement(geo, massVolElinv);

      ///////////////////////////////
      // Volumetric integral part
      ///////////////////////////////
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        // evaluate analytical flux and source 
        caller_.analyticalFluxAndSource(en, volQuad, l, fMat_, source_ );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                             * massVolElinv*volQuad.weight(l);
        
        for (int i = 0; i < updEn_numDofs; ++i) 
        {
          updEn[i] += 
            (bsetEn.evaluateGradientSingle(i, en, volQuad, l, fMat_) +
             bsetEn.evaluateSingle(i, volQuad, l, source_))*intel;
        }
      }

      /////////////////////////////
      // Surface integral part
      /////////////////////////////
      double dtLocal = 0.0;
      double nbvol;

      IntersectionIteratorType endnit = gridPart_.iend(en);
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      {
        double wspeedS = 0.0;
        if (nit.neighbor()) 
        {
          // get neighbor 
          EntityPointerType ep = nit.outside();
          EntityType & nb = *ep;
    
          if ((iset.index(nb) > iset.index(en) && en.level()==nb.level())
              || en.level() > nb.level()
              || nb.partitionType() != InteriorEntity) 
          {
            // for conforming situations apply Quadrature given
            if( twistUtil_.conforming(nit) )
            {
              const int twistSelf = twistUtil_.twistInSelf(nit); 
              FaceQuadratureType faceQuadInner(nit, faceQuadOrd_, twistSelf, 
                                             FaceQuadratureType::INSIDE);
        
              const int twistNeighbor = twistUtil_.twistInNeighbor(nit);
              FaceQuadratureType faceQuadOuter(nit, faceQuadOrd_, twistNeighbor,
                                             FaceQuadratureType::OUTSIDE);

              // apply neighbor part, return is volume of neighbor which is
              // needed below 
              nbvol = applyLocalNeighbor(nit,en,nb,massVolElinv,
                        faceQuadInner,faceQuadOuter,
                        bsetEn,updEn_numDofs,updEn,
                        dtLocal,wspeedS);
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

              const int twistSelf = twistUtil_.twistInSelf(nit); 
              NonConformingFaceQuadratureType ncFaceQuadInner(nit, faceQuadOrd_, twistSelf,
                                               NonConformingFaceQuadratureType::INSIDE);

              const int twistNeighbor = twistUtil_.twistInNeighbor(nit);
              NonConformingFaceQuadratureType ncFaceQuadOuter(nit, faceQuadOrd_, twistNeighbor,
                                               NonConformingFaceQuadratureType::OUTSIDE);

              // apply neighbor part, return is volume of neighbor which is
              // needed below 
              nbvol = applyLocalNeighbor(nit,en,nb,massVolElinv,
                        ncFaceQuadInner,ncFaceQuadOuter,
                        bsetEn,updEn_numDofs,updEn,
                        dtLocal,wspeedS);
            }

          }
            
        } // end if neighbor

        if (nit.boundary()) 
        {
          const int twistSelf = twistUtil_.twistInSelf(nit); 
          FaceQuadratureType faceQuadInner(nit, faceQuadOrd_, twistSelf, 
                                           FaceQuadratureType::INSIDE);
          nbvol = vol;
          caller_.setNeighbor(en);
          const int faceQuadInner_nop = faceQuadInner.nop();
          for (int l = 0; l < faceQuadInner_nop; ++l) 
          {
            double dtLocalS = 
              caller_.boundaryFlux(nit, faceQuadInner, l, source_);
            
            dtLocal += dtLocalS*faceQuadInner.weight(l);
            wspeedS += dtLocalS*faceQuadInner.weight(l);
                    
            for (int i = 0; i < updEn_numDofs; ++i) 
            {
              updEn[i] -= bsetEn.evaluateSingle(i, faceQuadInner, l, source_)
                *faceQuadInner.weight(l)*massVolElinv;
            }
          }
        } // end if boundary
        
        if (wspeedS>2.*std::numeric_limits<double>::min()) 
        {
          double minvolS = std::min(vol,nbvol);
          dtMin_ = std::min(dtMin_,minvolS/wspeedS);
        }
      }
    }

    template <class QuadratureImp>  
    double applyLocalNeighbor(IntersectionIteratorType & nit, 
            EntityType & en, EntityType & nb, 
            const double massVolElinv, 
            const QuadratureImp & faceQuadInner, 
            const QuadratureImp & faceQuadOuter,
            const BaseFunctionSetType & bsetEn, 
            const int updEn_numDofs, 
            LocalFunctionType & updEn,
            double & dtLocal, 
            double & wspeedS) const 
    {
      // make Entity known in caller  
      caller_.setNeighbor(nb);
      
      // get local function  
      LocalFunctionType updNeigh = dest_->localFunction(nb);
      const BaseFunctionSetType& bsetNeigh = updNeigh.baseFunctionSet();

      // get goemetry of neighbor 
      const GeometryType & nbGeo = nb.geometry();
      double massVolNbinv;
      double nbvol = volumeElement(nbGeo, massVolNbinv);
      
      const int faceQuadInner_nop = faceQuadInner.nop();
      for (int l = 0; l < faceQuadInner_nop; ++l) 
      {
        double dtLocalS = 
          caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,
                                l, valEn_, valNeigh_);
        
        dtLocal += dtLocalS*faceQuadInner.weight(l);
        wspeedS += dtLocalS*faceQuadInner.weight(l);

        for (int i = 0; i < updEn_numDofs; ++i) 
        {
          updEn[i] -= 
            bsetEn.evaluateSingle(i, faceQuadInner, l, valEn_)*
            faceQuadInner.weight(l)*massVolElinv;
          updNeigh[i] += 
            bsetNeigh.evaluateSingle(i, faceQuadOuter, l, valNeigh_)*
            faceQuadOuter.weight(l)*massVolNbinv;
        }
      }
      return nbvol;
    }
                         
  private:
    LocalDGPass();
    LocalDGPass(const LocalDGPass&);
    LocalDGPass& operator=(const LocalDGPass&);

  private:
    double volumeElement(const GeometryType& geo,
                         double& massVolinv) const 
    {
      double volume = geo.volume(); 

      typedef typename GeometryType :: ctype coordType; 
      enum { dim = GridType :: dimension };
      const ReferenceElement< coordType, dim > & refElem =
             ReferenceElements< coordType, dim >::general(geo.type());
      
      double volRef = refElem.volume();

      massVolinv = volRef/volume;
      return volume;
    }
    
  private:
    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    mutable CommunicationManagerType communicationManager_;

    mutable double dtMin_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable DomainType grads_;
    TimeProvider* time_;
    FieldVector<int, 0> diffVar_;

    TwistUtility<GridType> twistUtil_;

    int volumeQuadOrd_,faceQuadOrd_;
  };
//! @}  
} // end namespace Dune

#endif
