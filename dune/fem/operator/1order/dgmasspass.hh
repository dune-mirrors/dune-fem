#ifndef DUNE_FEM_DGMASSPASS_HH
#define DUNE_FEM_DGMASSPASS_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/common/memory.hh>
#include <dune/fem/pass/localdg.hh>
#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune
{

  namespace Fem
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
        typedef LocalDGPass<DiscreteModelImp, PreviousPassImp> BaseType;

      public:
        typedef typename BaseType::DiscreteModelType DiscreteModelType;
        typedef typename BaseType::PreviousPassType PreviousPassType;

        typedef typename BaseType::ArgumentType ArgumentType;

        typedef typename BaseType::DestinationType DestinationType;
        typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
        typedef typename BaseType::LocalFunctionType LocalFunctionType;

        typedef typename BaseType::GridPartType GridPartType;
        typedef typename BaseType::GridType GridType;
        typedef typename BaseType::IteratorType IteratorType;
        typedef typename BaseType::Entity EntityType;
        typedef typename BaseType::Geometry GeometryType;
        typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;


        typedef typename BaseType::RangeType RangeType;
        typedef typename RangeType::value_type RangeFieldType;

        typedef typename BaseType::VolumeQuadratureType VolumeQuadratureType;

        typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType;

        typedef typename DiscreteModelType::MassFactorType MassFactorType;

        enum { dimRange = BaseType::dimRange };

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
          spc_( referenceToSharedPtr( spc ) ),
          communicationManager_( space() ),
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
            IteratorType endit = space().end();
            for (IteratorType it = space().begin(); it != endit; ++it)
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
          //- statements
          this->caller_.setEntity(en);
          LocalFunctionType updEn = this->dest_->localFunction(en);
          const int updEn_numDofs = updEn.numDofs();
          const BasisFunctionSetType& bsetEn = updEn.basisFunctionSet();

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

        const DiscreteFunctionSpaceType &space () const { return *spc_; }

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
          const Dune::ReferenceElement< coordType, dim > & refElem =
                 Dune::ReferenceElements< coordType, dim >::general(geo.type());

          double volRef = refElem.volume();
          return volRef/volume;
        }

      private:
        DiscreteModelType& problem_;

        std::shared_ptr< const DiscreteFunctionSpaceType > spc_;
        mutable CommunicationManagerType communicationManager_;

        mutable RangeType tau_;
        mutable RangeType tauTmp_;
        mutable MassFactorType massVal_;
        const double factor_;

        mutable std::vector<RangeFieldType> massMatrix_;
      };

  //! @}
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DGMASSPASS_HH
