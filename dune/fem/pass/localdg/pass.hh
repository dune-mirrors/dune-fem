#ifndef DUNE_FEM_PASS_LOCALDG_HH
#define DUNE_FEM_PASS_LOCALDG_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/dgpass.hh> instead!"
#endif

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/common/local.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/storage/dynamicarray.hh>

#include "modelcaller.hh"

namespace Dune
{

  namespace Fem
  {

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

      //! pass ids up to here (tuple of integral constants)
      typedef typename BaseType::PassIds PassIds;

      //! Repetition of template arguments
      typedef DiscreteModelImp DiscreteModelType;
      //! Repetition of template arguments
      typedef PreviousPassImp PreviousPassType;

      // Types from the base class
      typedef typename BaseType::EntityType EntityType;
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
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

      // Types extracted from the underlying grids
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType::Intersection IntersectionType;
      typedef typename GridType::template Codim<0>::Geometry Geometry;


      // Various other types
      typedef typename DestinationType::LocalFunctionType LocalFunctionType;

      typedef DGDiscreteModelCaller< DiscreteModelType, ArgumentType, PassIds > DiscreteModelCallerType;

      // Range of the destination
      enum { dimRange = DiscreteFunctionSpaceType::dimRange };

      // type of local id set
      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef TemporaryLocalFunction< DiscreteFunctionSpaceType > TemporaryLocalFunctionType;

      //! type of local mass matrix
      typedef LocalMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;

    public:
      //- Public methods
      //! Constructor
      //! \param problem Actual problem definition (see problem.hh)
      //! \param pass Previous pass
      //! \param spc Space belonging to the discrete function local to this pass
      //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
      //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
      //! \param notThreadParallel  true if pass is used in single thread mode
      LocalDGPass(DiscreteModelType& problem,
                  PreviousPassType& pass,
                  const DiscreteFunctionSpaceType& spc,
                  const int volumeQuadOrd =-1,
                  const int faceQuadOrd=-1,
                  const bool notThreadParallel = true ) :
        BaseType(pass, spc),
        caller_(0),
        problem_(problem),
        arg_(0),
        dest_(0),
        gridPart_( space().gridPart()),
        indexSet_(gridPart_.indexSet()),
        visited_(0),
        updEn_( space() ),
        updNb_( space() ),
        fMatVec_( 20 ),
        valEnVec_( 20 ),
        valNbVec_( 20 ),
        dtMin_(std::numeric_limits<double>::max()),
        minLimit_(2.0*std::numeric_limits<double>::min()),
        volumeQuadOrd_( (volumeQuadOrd < 0) ?  ( 2 * space().order()) : volumeQuadOrd ),
        faceQuadOrd_( (faceQuadOrd < 0) ?  ( 2 * space().order()+1) : faceQuadOrd ),
        localMassMatrix_( space() , volumeQuadOrd_ ),
        notThreadParallel_( notThreadParallel )
      {
        fMatVec_.setMemoryFactor( 1.1 );
        valEnVec_.setMemoryFactor( 1.1 );
        valNbVec_.setMemoryFactor( 1.1 );

        assert( volumeQuadOrd_ >= 0 );
        assert( faceQuadOrd_ >= 0 );
      }

      LocalDGPass ( const LocalDGPass & ) = delete;

      //! Destructor
      virtual ~LocalDGPass() = default;

      LocalDGPass &operator= ( const LocalDGPass & ) = delete;

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
        const double p = 2 * space().order() + 1;
        return dtMin_ / p;
      }

    public:
      //! In the preparations, store pointers to the actual arguments and
      //! destinations. Filter out the "right" arguments for this pass.
      virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
      {
        arg_ = const_cast<ArgumentType*>(&arg);
        dest_ = &dest;

        if( notThreadParallel_ )
        {
          // clear destination
          dest_->clear();
        }

        caller_ = new DiscreteModelCallerType( *arg_, problem_ );
        assert( caller_ );
        caller_->setTime( this->time() );

        // resize indicator function
        visited_.resize( indexSet_.size(0) );
        // set all values to false
        const int indSize = visited_.size();
        for(int i=0; i<indSize; ++i) visited_[i] = false;

        // time initialisation to max value
        dtMin_ = std::numeric_limits<double>::max();
      }

      //! Some timestep size management.
      virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
      {
        if( notThreadParallel_ )
        {
          // communicate calculated function
          space().communicate( dest );
        }

        if( caller_ )
          delete caller_;
        caller_ = 0;
      }

      size_t numberOfElements() const { return 0; }

      struct DefaultNBChecker
      {
        bool operator ()(const EntityType& , const EntityType& ) const
        {
          return true ;
        }
      };

      void applyLocal( const EntityType& en) const
      {
        applyLocal( en, DefaultNBChecker() );
      }

      void applyLocalMass( const EntityType& en) const
      {
        // get local function and add update
        LocalFunctionType function = dest_->localFunction(en);

        // apply local inverse mass matrix
        localMassMatrix_.applyInverse( caller(), en, function );
      }

      //! local integration
      template <class NeighborChecker>
      void applyLocal( const EntityType& en, const NeighborChecker& nbChecker ) const
      {
        // init local function
        initLocalFunction( en , updEn_ );

        // call real apply local
        applyLocal(en, updEn_, nbChecker );

        // add update to real function
        updateFunctionAndApplyMass(en, updEn_ );
      }

      using BaseType::space;

    protected:
      //! local integration
      template <class NeighborChecker>
      void applyLocal( const EntityType& en,
                       TemporaryLocalFunctionType& updEn,
                       const NeighborChecker& nbChecker ) const
      {
        // only call geometry once, who know what is done in this function
        const Geometry & geo = en.geometry();

        // get volume of element
        const double vol = geo.volume();

        // only apply volumetric integral if order > 0
        // otherwise this contribution is zero

        if( (space().order() > 0) || problem_.hasSource() )
        {
          // create quadrature
          VolumeQuadratureType volQuad( en, volumeQuadOrd_ );

          // set entity and evaluate local functions for quadrature
          caller().setEntity( en , volQuad );

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
          caller().setEntity( en );
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
              EntityType nb = intersection.outside();

              const bool canUpdateNeighbor = nbChecker( en, nb );

              if ( ! visited_[ indexSet_.index( nb ) ] )
              {
                // for conforming situations apply Quadrature given
                if( !GridPartCapabilities::isConforming< GridPartType >::v
                    && !intersection.conforming() )
                {
                  // apply neighbor part, return is volume of neighbor which is
                  // needed below
                  nbvol = applyLocalNeighbor< false >
                                      (intersection,en,nb,
                                       updEn, updNb_ ,
                                       wspeedS,
                                       canUpdateNeighbor );
                }
                else
                {
                  // apply neighbor part, return is volume of neighbor which is
                  // needed below
                  nbvol = applyLocalNeighbor< true >
                                      (intersection,en,nb,
                                       updEn, updNb_ ,
                                       wspeedS,
                                       canUpdateNeighbor );
                }

              } // end if do something

            } // end if neighbor
            else if( intersection.boundary() )
            {
              FaceQuadratureType faceQuadInner(gridPart_, intersection, faceQuadOrd_,
                                               FaceQuadratureType::INSIDE);

              // set neighbor entity to inside entity
              caller().setBoundary(en, faceQuadInner);

              // cache number of quadrature points
              const size_t faceQuadInner_nop = faceQuadInner.nop();

              if( valEnVec_.size() < faceQuadInner_nop )
                valEnVec_.resize( faceQuadInner_nop );

              // loop over quadrature points
              for (size_t l = 0; l < faceQuadInner_nop; ++l)
              {
                RangeType& flux = valEnVec_[ l ];

                // eval boundary Flux
                wspeedS += caller().boundaryFlux( intersection, faceQuadInner, l, flux )
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
        localMassMatrix_.applyInverse( caller(), en, function );
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
        const size_t volQuad_nop = volQuad.nop();
        if( fMatVec_.size() < volQuad_nop )
        {
          fMatVec_.resize( volQuad_nop );
        }

        for (size_t l = 0; l < volQuad_nop; ++l)
        {
          JacobianRangeType& flux = fMatVec_[ l ];

          // evaluate analytical flux and source
          caller().analyticalFlux(en, volQuad, l, flux );

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
        const size_t volQuad_nop = volQuad.nop();

        if( fMatVec_.size() < volQuad_nop )
        {
          fMatVec_.resize( volQuad_nop );
        }

        if( valEnVec_.size() < volQuad_nop )
        {
          valEnVec_.resize( volQuad_nop );
        }

        for (size_t l = 0; l < volQuad_nop; ++l)
        {
          JacobianRangeType& flux = fMatVec_[ l ];
          RangeType& source = valEnVec_[ l ];

          // evaluate analytical flux and source
          const double dtEst =
            caller().analyticalFluxAndSource(en, volQuad, l, flux, source );

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
                                  double &wspeedS,
                                  const bool canUpdateNeighbor ) const
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
        caller().setNeighbor(nb, faceQuadInner, faceQuadOuter);

        // get goemetry of neighbor
        const Geometry & nbGeo = nb.geometry();

        // get neighbors volume
        const double nbvol = nbGeo.volume();

        const size_t faceQuadInner_nop = faceQuadInner.nop();
        assert( faceQuadInner.nop() == faceQuadOuter.nop() );

        // check  valNbVec_ here, valEnVec_ might have been resized
        if( valNbVec_.size() < faceQuadInner_nop )
        {
          valEnVec_.resize( faceQuadInner_nop );
          valNbVec_.resize( faceQuadInner_nop );
        }

        for (size_t l = 0; l < faceQuadInner_nop; ++l)
        {
          RangeType& fluxEn = valEnVec_[ l ];
          RangeType& fluxNb = valNbVec_[ l ];

          wspeedS += caller().numericalFlux(intersection,
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

        // update on neighbor
        if( canUpdateNeighbor )
        {
          // init local function
          initLocalFunction( nb, updNb );
          // add fluxes
          updNb.axpyQuadrature( faceQuadOuter, valNbVec_ );
          // add update to real function
          updateFunction(nb, updNb );
        }

        return nbvol;
      }

    protected:
      DiscreteModelCallerType &caller () const
      {
        assert( caller_ );
        return *caller_;
      }

    protected:
      mutable DiscreteModelCallerType *caller_;
      DiscreteModelType& problem_;

      mutable ArgumentType* arg_;
      mutable DestinationType* dest_;

      const GridPartType & gridPart_;
      const IndexSetType& indexSet_;

      // indicator for grid walk
      mutable DynamicArray<bool> visited_;

      mutable TemporaryLocalFunctionType updEn_;
      mutable TemporaryLocalFunctionType updNb_;

      //! Some helper variables
      mutable DynamicArray< JacobianRangeType > fMatVec_;
      mutable DynamicArray< RangeType > valEnVec_;
      mutable DynamicArray< RangeType > valNbVec_;

      mutable double dtMin_;
      const double minLimit_;

      const int volumeQuadOrd_, faceQuadOrd_;
      LocalMassMatrixType localMassMatrix_;
      const bool notThreadParallel_;
    };
  //! @}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_LOCALDG_HH
