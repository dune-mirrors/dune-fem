#ifndef DUNE_DGMASSPASS_HH
#define DUNE_DGMASSPASS_HH

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

//- Dune-fem includes 
#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune
{

/*! @ingroup GradientOperator 
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   u &=& A(x)\nabla p \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ p \f$ is the argument and \f$ u \f$ is computed.
** @{
**************************************************************************/

  /** \brief Implementation of operator to calculate gradient of 
      a given discrete function using the pass concept.
  */
  template< class DiscreteModelImp, class PreviousPassImp >
  class LocalDGMassPass
  : public LocalDGPass< DiscreteModelImp, PreviousPassImp >
  {
    typedef LocalDGMassPass< DiscreteModelImp, PreviousPassImp > ThisType;

  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalDGPass<DiscreteModelImp, PreviousPassImp> BaseType;

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
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType, SelectorType > CombinedSelectorType;
    typedef DiscreteModelCaller< DiscreteModelType, ArgumentType, CombinedSelectorType >
      DiscreteModelCallerType;

    // type of Communication Manager 
    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType;
    
    typedef typename DiscreteModelType :: MassFactorType MassFactorType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param factor mass matrix factor, default is -1.0
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGMassPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                double factor = -1.0,
                int volumeQuadOrd = -1, int faceQuadOrd=-1) :
      BaseType(problem, pass, spc,volumeQuadOrd,faceQuadOrd),
      problem_(problem),
      spc_(spc),
      communicationManager_(spc_),
      tau_(0.0),
      tauTmp_(0.0),
      massVal_(0.0),
      factor_(factor)
    {
    }
   
    //! Destructor
    virtual ~LocalDGMassPass() {
    }

  private:
    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      if(problem_.hasMass())
      {
        IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
          applyLocalMass(*it);
        }
      }

      // call finalize of dg pass (i.e. data communication)
      BaseType :: finalize(arg,dest);
    }

  private:    
    // apply mass matrix multiplication
    // here matrix free implementation due to memory savings
    void applyLocalMass(EntityType& en) const
    {
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;

      //- statements
      this->caller_.setEntity(en);
      LocalFunctionType updEn = this->dest_->localFunction(en);
      const int updEn_numDofs = updEn.numDofs();
      const BaseFunctionSetType& bsetEn = updEn.baseFunctionSet(); 
      
      // only call geometry once, who know what is done in this function 
      const GeometryType & geo = en.geometry();

      const double massVolElinv = massVolInv(geo);

      if((int)massMatrix_.size() != updEn_numDofs)
      {
        massMatrix_.resize(updEn_numDofs);
      }

      // clear mass entries 
      for (int i = 0; i < updEn_numDofs; ++i) 
      {
        massMatrix_[i] = 0.0;
      }
      
      ///////////////////////////////
      // Volumetric integral part
      ///////////////////////////////
      VolumeQuadratureType volQuad(en, this->volumeQuadOrd_);
      const int volQuad_nop = volQuad.nop();
      for (int l = 0; l < volQuad_nop; ++l) 
      {
        const double intel = geo.integrationElement(volQuad.point(l))
                             * massVolElinv * volQuad.weight(l);

        // evaluate mass factor 
        this->caller_.mass(en, volQuad, l, massVal_ );
        
        for (int i = 0; i < updEn_numDofs; ++i) 
        {
          // eval tau_k 
          bsetEn.evaluate(i, volQuad[l], tau_ );

          // apply mass factor 
          massVal_.mv(tau_,tauTmp_);
            
          massMatrix_[i] += bsetEn.evaluateSingle(i, volQuad[l], tauTmp_ ) * intel;
        }
      }

      // multiply with mass matrix 
      for(int i=0; i<updEn_numDofs; ++i)
      {
        updEn[i] *= factor_ * massMatrix_[i];
      }
    }

  private:
    LocalDGMassPass();
    LocalDGMassPass(const LocalDGMassPass&);
    LocalDGMassPass& operator=(const LocalDGMassPass&);

  private:
    double massVolInv(const GeometryType& geo) const 
    {
      double volume = geo.volume(); 

      typedef typename GeometryType :: ctype coordType; 
      enum { dim = GridType :: dimension };
      const ReferenceElement< coordType, dim > & refElem =
             ReferenceElements< coordType, dim >::general(geo.type());
      
      double volRef = refElem.volume();
      return volRef/volume;
    }
    
  private:
    DiscreteModelType& problem_;     

    const DiscreteFunctionSpaceType& spc_;
    mutable CommunicationManagerType communicationManager_;

    mutable RangeType tau_;
    mutable RangeType tauTmp_;
    mutable MassFactorType massVal_;
    const double factor_;

    mutable std::vector<RangeFieldType> massMatrix_;
  };

//! @}  
} // end namespace Dune

#endif
