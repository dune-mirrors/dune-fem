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

//- system includes 
#include <map>

#include <dune/fem/misc/utility.hh>

#include "pass.hh"
#include "selection.hh"
#include "discretemodel.hh"
#include "modelcaller.hh"

// * needs to move
#include <dune/fem/misc/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/allgeomtypes.hh> 

namespace Dune {
/*! @addtogroup PassHyp
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

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // type of Communication Manager 
    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };

    // type of local id set 
    typedef typename GridType::Traits::LocalIdSet LocalIdSetType; 

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
      localIdSet_(spc_.grid().localIdSet()),
      communicationManager_(spc_),
      refVolMap_ (),
      dtMin_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      grads_(0.0),
      time_(0),
      minLimit_(2.0*std::numeric_limits<double>::min()),
      diffVar_(),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? 
          (2*spc_.order()) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? 
        (2*spc_.order()+1) : faceQuadOrd )
    {
      {
        typedef AllGeomTypes< typename GridPartType :: IndexSetType,
                              GridType> AllGeomTypesType;
        AllGeomTypesType allGeomTypes( gridPart_.indexSet() );
        const std::vector<Dune::GeometryType>& geomTypes =
          allGeomTypes.geomTypes(0);

        for(size_t i=0; i<geomTypes.size(); ++i)
        {
          typedef typename Geometry :: ctype coordType; 
          enum { dim = GridType :: dimension };
          const ReferenceElement< coordType, dim > & refElem =
                 ReferenceElements< coordType, dim >::general( geomTypes[i] );
      
          // store volume of reference element
          refVolMap_[ geomTypes[i] ] = refElem.volume();
        }
      }

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

  protected:
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
      if (time_) 
      {
        caller_.setTime(time_->time());
      }
      else 
      {
        caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // communicate calculated function 
      communicationManager_.exchange( dest );
      
      // if time provider exists, check dtMin_
      if (time_) 
      {
        time_->provideTimeStepEstimate(dtMin_);
      }
      
      // call finalize 
      caller_.finalize();
    }

    //! local integration 
    void applyLocal(EntityType& en) const
    {
      //- statements
      caller_.setEntity(en);
      LocalFunctionType updEn = dest_->localFunction(en);

      // only call geometry once, who know what is done in this function 
      const Geometry & geo = en.geometry();

      const double vol = geo.volume(); 
      assert( refVolMap_.find (  geo.type() ) != refVolMap_.end() );
      const double massVolElinv = refVolMap_[ geo.type() ] / vol;

      // only apply volumetric integral if order > 0 
      // otherwise this contribution is zero 
      
      if( (spc_.order() > 0) || problem_.hasSource()) 
      {
        // if only flux, evaluate only flux 
        if ( problem_.hasFlux() && !problem_.hasSource() ) 
        {
          evalVolumetricPartFlux(en, geo, updEn , massVolElinv);
        }
        else 
        {
          // evaluate flux and source 
          evalVolumetricPartBoth(en, geo, updEn , massVolElinv);
        }
      }

      /////////////////////////////
      // Surface integral part
      /////////////////////////////
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
    
          if (localIdSet_.id(en) < localIdSet_.id(nb) && en.level()==nb.level() 
              || en.level() > nb.level()
#if HAVE_MPI 
              || nb.partitionType() != InteriorEntity
#endif
              ) 
          {
            typedef TwistUtility<GridType> TwistUtilityType;
            // for conforming situations apply Quadrature given
            if( TwistUtilityType::conforming(gridPart_.grid(),nit) )
            {
              FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                               FaceQuadratureType::INSIDE);
        
              FaceQuadratureType faceQuadOuter(gridPart_, nit, faceQuadOrd_,
                                               FaceQuadratureType::OUTSIDE);
              // apply neighbor part, return is volume of neighbor which is
              // needed below 
              nbvol = applyLocalNeighbor(nit,en,nb,massVolElinv,
                        faceQuadInner,faceQuadOuter,
                        updEn,
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
                  nonConformingFaceQuadInner(gridPart_, nit, faceQuadOrd_,
                                             NonConformingFaceQuadratureType::INSIDE);

              NonConformingFaceQuadratureType 
                  nonConformingFaceQuadOuter(gridPart_,nit, faceQuadOrd_,
                                             NonConformingFaceQuadratureType::OUTSIDE);

              // apply neighbor part, return is volume of neighbor which is
              // needed below 
              nbvol = applyLocalNeighbor(nit,en,nb,massVolElinv,
                        nonConformingFaceQuadInner,
                        nonConformingFaceQuadOuter,
                        updEn,
                        wspeedS);
            }

          }
            
        } // end if neighbor

        if (nit.boundary()) 
        {
          FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_, 
                                           FaceQuadratureType::INSIDE);
          nbvol = vol;
          caller_.setNeighbor(en);
          const int faceQuadInner_nop = faceQuadInner.nop();
          for (int l = 0; l < faceQuadInner_nop; ++l) 
          {
            // eval boundary Flux  
            wspeedS += caller_.boundaryFlux(nit, faceQuadInner, l, valEn_ )
                     * faceQuadInner.weight(l);
            
            // apply weights 
            valEn_ *= -faceQuadInner.weight(l) * massVolElinv;

            // add factor 
            updEn.axpy( faceQuadInner, l, valEn_  );
          }
        } // end if boundary
        
        if (wspeedS > minLimit_ ) 
        {
          double minvolS = std::min(vol,nbvol);
          dtMin_ = std::min(dtMin_,minvolS/wspeedS);
        }
      }
    }

    //////////////////////////////////////////
    // Volumetric integral part only flux 
    //////////////////////////////////////////
    void evalVolumetricPartFlux(EntityType& en , const Geometry& geo , 
        LocalFunctionType& updEn , const double massVolElinv ) const
    {
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        // evaluate analytical flux and source 
        caller_.analyticalFlux(en, volQuad, l, fMat_ );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                             * massVolElinv * volQuad.weight(l);
        
        fMat_ *= intel;

        // add fMat 
        updEn.axpy( volQuad, l, fMat_);
      }
    }
    
    //////////////////////////////////////////
    // Volumetric integral part only flux 
    //////////////////////////////////////////
    void evalVolumetricPartBoth(EntityType& en , const Geometry& geo , 
        LocalFunctionType& updEn , const double massVolElinv ) const
    {
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        // evaluate analytical flux and source 
        caller_.analyticalFluxAndSource(en, volQuad, l, fMat_, source_ );
        
        const double intel = geo.integrationElement(volQuad.point(l))
                             * massVolElinv*volQuad.weight(l);
        
        source_ *= intel;
        fMat_   *= intel;
        
        updEn.axpy(volQuad,l,source_,fMat_);
      }
    }
    
    template <class QuadratureImp>  
    double applyLocalNeighbor(IntersectionIteratorType & nit, 
            EntityType & en, EntityType & nb, 
            const double massVolElinv, 
            const QuadratureImp & faceQuadInner, 
            const QuadratureImp & faceQuadOuter,
            LocalFunctionType & updEn,
            double & wspeedS) const 
    {
      // make Entity known in caller  
      caller_.setNeighbor(nb);
      
      // get local function  
      LocalFunctionType updNeigh = dest_->localFunction(nb);

      // get goemetry of neighbor 
      const Geometry & nbGeo = nb.geometry();

      const double nbvol = nbGeo.volume(); 
      assert( refVolMap_.find (  nbGeo.type() ) != refVolMap_.end() );
      const double massVolNbinv = refVolMap_[ nbGeo.type() ] / nbvol;
      
      const int faceQuadInner_nop = faceQuadInner.nop();
      for (int l = 0; l < faceQuadInner_nop; ++l) 
      {
        wspeedS += caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,l, valEn_, valNeigh_)
                 * faceQuadInner.weight(l);
        
        // apply weights 
        valEn_    *= -faceQuadInner.weight(l)*massVolElinv;
        valNeigh_ *=  faceQuadOuter.weight(l)*massVolNbinv;
        
        updEn.axpy( faceQuadInner, l, valEn_ );
        updNeigh.axpy( faceQuadOuter, l, valNeigh_ );
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
    const LocalIdSetType& localIdSet_;
    mutable CommunicationManagerType communicationManager_;
    mutable RefVolumeMapType refVolMap_;

    mutable double dtMin_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_,fMatTmp;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable DomainType grads_;
    TimeProvider* time_;
    const double minLimit_;
    FieldVector<int, 0> diffVar_;

    const int volumeQuadOrd_, faceQuadOrd_;
  };
//! @}  
} // end namespace Dune

#endif
