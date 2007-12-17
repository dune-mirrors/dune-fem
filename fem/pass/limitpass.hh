#ifndef DUNE_LIMITERPASS_HH
#define DUNE_LIMITERPASS_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/common/utility.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>


#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/misc/timeprovider.hh>

#include <dune/fem/pass/ldgflux.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>

#include <dune/fem/operator/artificialdiffusion.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>

//*************************************************************
namespace Dune {  
/*! @addtogroup PassLimit
*/
  template <class GlobalTraitsImp, class Model>
  class LimiterDefaultDiscreteModel;
  
  template <class GlobalTraitsImp, class Model>
  struct LimiterDefaultTraits 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef GlobalTraitsImp Traits; 
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef LimiterDefaultDiscreteModel<GlobalTraitsImp,Model> DiscreteModelType;
  };
  
  
  // **********************************************
  // **********************************************
  // **********************************************
  template <class GlobalTraitsImp, class Model>
  class LimiterDefaultDiscreteModel :
    public DiscreteModelDefaultWithInsideOutSide< LimiterDefaultTraits<GlobalTraitsImp,Model> > 
  {
  public:
    typedef LimiterDefaultTraits<GlobalTraitsImp,Model> Traits;
    
    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    
  public:
    LimiterDefaultDiscreteModel(const Model& mod) 
      : model_(mod) , velocity_(0) {}

    bool hasSource() const { return false; }
    bool hasFlux() const   { return true;  }
    
    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight) const
    { 
      
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
 
      if( checkDirection(it,time,x,uLeft,uRight,gLeft,gRight) )
      {
        gLeft  = argULeft;
        gLeft -= argURight;
        gRight = gLeft;
        return it.integrationOuterNormal( x ).two_norm();
      }
      else 
      {
        gLeft = gRight = 0.0;
        return 0.0;
      }
    }

  protected:
    template <class ArgumentTuple>
    bool checkDirection(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft, 
                        const ArgumentTuple& uRight,
                        RangeType& gLeft,
                        RangeType& gRight) const 
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      model_.velocity(this->inside(),time,it.intersectionSelfLocal().global(x),
                      argULeft,velocity_);
      return ((it.outerNormal(x) * velocity_) < 0);
    }

  protected:
    const Model& model_;
    mutable DomainType velocity_;
  };


  /** \brief Concrete implementation of Pass for Limiting.
      The implemented Shock detection is described in detail in: 
        L. Krivodonova and J. Xin and J.-F. Remacle and N. Chevaugeon and J. E. Flaherty
        Shock detection and limiting with discontinuous Galerkin methods for hyperbolic conservation laws.
        Appl. Numer. Math., 48(3-4), pages 323-338, 2004.
    
      Link to paper:
        http://www.scorec.rpi.edu/REPORTS/2003-3.pdf

      Limiting is done by simply setting the polynomial order to zero.
  */
  template <class DiscreteModelImp, class PreviousPassImp>
  class LimitDGPass :
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
    typedef typename BaseType::ArgumentType ArgumentType;
    
    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType 
    JacobianRangeType;
    
    // Types extracted from the underlying grids
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;
    
    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    typedef DofConversionUtility< PointBased > DofConversionUtilityType;

    typedef LeafGridPart< GridType > ContGridPartType; 

    typedef LagrangeDiscreteFunctionSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , ContGridPartType, 1 > LagrangeSpaceType;

    typedef DiscontinuousGalerkinSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , GridPartType, 0 > DG0SpaceType;

    typedef AdaptiveDiscreteFunction< LagrangeSpaceType > P1FunctionType; 
    typedef AdaptiveDiscreteFunction< DG0SpaceType > DG0FunctionType; 

    typedef ArtificialDiffusion<DestinationType> ArtificialDiffusionType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange,
     dimDomain = DiscreteFunctionSpaceType::DimDomain};
    enum { dimension = GridType :: dimension };
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LimitDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc) :
      BaseType(pass, spc),
      caller_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      orderPower_( -((spc_.order()+1.0)/2.0)),
      dofConversion_(dimRange),
      faceQuadOrd_(spc_.order()),
      jump_(0),
      jump2_(0),
      gradientFlux_(1.0,10.0),
      artDiff_(0),
      comboSet_()  
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

          // store local coordinates of barycenter 
          baryCenterMap_[ geomTypes[i] ] = refElem.position(0,0);
          //std::cout << "Local barycenter = " << baryCenterMap_[ geomTypes[i] ] << "\n";
        }
      }

      // we need the flux here 
      assert(problem.hasFlux());
    }
    
    //! Destructor
    virtual ~LimitDGPass() {
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
    }
    
    //! Some management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      caller_.finalize();
      DestinationType* U = const_cast<DestinationType*>(Element<0>::get(*arg_));
      U->assign(dest);
    }

    //! Perform the limitation on all elements.
    virtual void applyLocal(EntityType& en) const
    {
      enum { dim = EntityType :: dimension };
      // check argument is not zero
      assert( arg_ );
      
      //- statements
      // set entity to caller 
      caller_.setEntity(en);

      // get function to limit 
      const DestinationType& U = *(Element<0>::get(*arg_));

      // get U on entity
      const LocalFunctionType uEn = U.localFunction(en);

      // get local funnction for limited values
      LocalFunctionType limitEn = dest_->localFunction(en);

      // get geometry 
      const Geometry& geo = en.geometry();
      
      // number of scalar base functions
      const int numBasis = limitEn.numDofs()/dimRange;
      
      double radius = 0.0;
      const GeometryType geomType = geo.type();
      const DomainType enBary = geo.global( baryCenterMap_[geomType] );
      
      if ( geomType.isSimplex() )
      {
        const int numcorners = geo.corners();
        DomainType diff;
        for(int n=0; n<numcorners; ++n)
        {
          diff = geo[(n+1)%numcorners] - geo[n];
          radius = std::max(diff.two_norm(),radius);
        }
      }
      else if ( geomType.isCube() )
      {
        // map of used geometry points to determ side length (see
        // reference elements)
        const double hypo = SQR((geo[1][0] - geo[0][0])*0.5) + SQR((geo[2][1] - geo[0][1])*0.5);

        if( dim == 3 )
        {
          radius = sqrt(hypo + SQR((geo[4][2] - geo[0][2])* 0.5));
        } 
        else if ( dim == 2 )
        {
          // radium of circum 
          radius = sqrt(hypo);
        }
        else 
        {
          radius = std::abs((geo[1][0] - geo[0][0])*0.5);
        }
      }
      else 
      {
        DUNE_THROW(NotImplemented,"Unsupported geometry type!");
      }
            
      const double hPowPolOrder = (1.0/(geo.volume())) * pow(radius, orderPower_);// -((order+1.0)/2.0) );
      //const double hPowPolOrder = pow(4.0*M_SQRT2*radius, orderPower_);// -((order+1.0)/2.0) );
      // get value of U in barycenter

      FieldVector<bool,dimRange> limit(false);

      RangeType totaljump(0);

      // calculate circume during neighbor check 
      double circume = 0.0;

      JacobianRangeType enGrad;
      JacobianRangeType nbGrad;

      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
           nit != endnit; ++nit) 
      {
        // check all neighbors 
        if (nit.neighbor()) 
        {
          // get neighbor entity
          EntityPointerType ep = nit.outside();
          EntityType& nb = *ep; 

          // set neighbor to caller 
          caller_.setNeighbor(nb);
          
          typedef TwistUtility<GridType> TwistUtilityType;
          // conforming case 
          if( TwistUtilityType::conforming(gridPart_.grid(),nit) )
          {
            FaceQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
            FaceQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

            applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump,circume);
          }
          else  
          { // non-conforming case 
            typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
            NonConformingQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
            NonConformingQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

            applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump,circume);
          }
        }

        // check all neighbors 
        if (nit.boundary()) 
        {
          FaceQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
          applyBoundary(nit,faceQuadInner,totaljump,circume);
        }
      } // end intersection iterator 
       
      circume = (circume > 0.0) ? 1.0/circume : 0.0;

      // determ whether limitation is necessary  
      bool limiter = false;
      for (int r=0; r<dimRange; ++r) 
      {
        double jumpr = std::abs(totaljump[r]);
        if ( jumpr*hPowPolOrder*circume > 1 ) 
        {
          limit[r] = true;
          limiter = true;
        }
        else
          limit[r] = false;
      }
       
      // prepare limitEn 
      for (int r=0; r<dimRange; ++r) 
      {
        if (limit[r]) 
        {
          // set to zero
          {
            const int dofIdx = dofConversion_.combinedDof(0,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
          
          // set to zero
          for (int i=1; i<numBasis; ++i) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = 0.;
          }
        }
        else 
        {
          // copy function since we do not limit anything
          for (int i=0; i<numBasis; ++i) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
        }
      }

      if( limiter )
      {
        // default is zero 
        typedef FieldVector< DomainType , dimRange > DeoModType; 
        DeoModType deoMod;

        std::vector< DeoModType > deoMods;
        deoMods.reserve( 20 );
        std::vector< CheckType > comboVec;
        comboVec.reserve( 20 );

        std::vector< DomainType > barys;
        std::vector< RangeType >  nbVals;

        barys.reserve( 6 );
        nbVals.reserve( 6 );

        // get quadrature for barycenter 
        VolumeQuadratureType volQuad( en, 0 );

        RangeType enVal; 
        // evalaute function  
        uEn.evaluate( volQuad[0] , enVal );

        IntersectionIteratorType endnit = gridPart_.iend(en); 
        for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
             nit != endnit; ++nit) 
        {
          // check all neighbors 
          if (nit.neighbor()) 
          {
            EntityPointerType ep = nit.outside();
            EntityType& nb = *ep;
            
            // get U on entity
            const LocalFunctionType uNb = U.localFunction(nb);

            // get quadrature for barycenter 
            VolumeQuadratureType nbQuad( nb, 0 );

            RangeType nbVal;
            // evaluate function  
            uNb.evaluate( nbQuad[0] , nbVal);

            // calculate difference 
            nbVal -= enVal;

            // store value 
            nbVals.push_back(nbVal);

            const Geometry& nbGeo = nb.geometry();
            
            DomainType nbBary = nbGeo.global( baryCenterMap_[ nbGeo.type() ] );
            // calculate difference 
            nbBary -= enBary;

            // store value 
            barys.push_back(nbBary);

          } // end neighbor 

          // use ghost cell for limiting 
          if( nit.boundary() )
          {
            // assume same value on ghost element
            nbVals.push_back(RangeType(0));

            DomainType tmp (nit.intersectionGlobal()[0]);
            tmp -= enBary;
            tmp *= 2.0;

            FieldVector<double,dim-1> faceMid(0.5);

            DomainType nbBary = nit.unitOuterNormal(faceMid);
            const double factor = nbBary * tmp;

            // calculate new value 
            nbBary *= factor;

            // store value 
            barys.push_back(nbBary);
          }
          
        } // end intersection iterator 

        // get number of found neighbors 
        const size_t neighbors = nbVals.size();

        {
          bool allNegative = true;
          bool allPositive = true;
          // check whether all differences are positive or negative 
          for(size_t n=0; n<neighbors; ++n) 
          {
            for(int r=0; r<dimRange; ++r)
            {
              if(limit[r] ) 
              {
                if( nbVals[n][r] >= 0 ) allNegative = false;
                if( nbVals[n][r] <= 0 ) allPositive = false;
              }
            }
          }
          
          // stay with the p-zero version  
          if( allNegative || allPositive ) return ;
        }

        typedef FieldMatrix<double, dim, dim> MatrixType;
        MatrixType matrix(0);
        MatrixType inverse(0);
        DomainType rhs(0);
        
        // create combination set 
        setupComboSet( neighbors );

        // calculate linear functions 
        // D(x) = U_i + D_i * (x - w_i)
        typedef typename ComboSetType :: const_iterator iterator; 
        const iterator endit = comboSet_.end();
        for(iterator it = comboSet_.begin(); it != endit; ++it) 
        {
          // get tuple of numbers 
          const KeyType& v = (*it).first; 

          // setup matrix 
          for(int i=0; i<dim; ++i) 
          {
            matrix[i] = barys[ v[i] ];
          }

          // create new instance of limiter coefficients  
          DeoModType dM;

          const double det = FMatrixHelp :: invertMatrix(matrix,inverse);
          if( std::abs( det ) > 0 )
          {
            // calculate D
            for(int r=0; r<dimRange; ++r)
            {
              for(int i=0; i<dim; ++i) 
              {
                rhs[i] = nbVals[ v[i] ][r];
              }

              DomainType& D = dM[r];

              // set to zero 
              D = 0;
              // get solution 
              inverse.umv( rhs, D );
            }

            // store linear function 
            deoMods.push_back( dM );
            comboVec.push_back( (*it).second );
          }
#if 1
          else 
          {
            // apply least square by adding another point 
            // this should make the linear system solvable 
             
            std::vector<int> nV( dim );
            for(int i=0; i<dim; ++i) nV[i] = v[i];

            // take first point of list of pionts to check 
            CheckType check ( (*it).second );
            assert( check.size() > 0 );

            bool matrixSingular = true ;
            int count = 1; 
            while ( matrixSingular && check.size() > 0 ) 
            {
              // take new entry
              nV.push_back(check[0]);
              // remove point from points to check list 
              check.erase( check.begin () );

              if( count == 1 ) 
              {
                matrixSingular = 
                  applyLeastSquare<dim+1>( barys,
                                    nbVals,
                                    nV,
                                    dM,
                                    matrix,
                                    inverse,
                                    rhs );
              }
              else if( count == 2 ) 
              {
                matrixSingular = 
                  applyLeastSquare<dim+2>( barys,
                                    nbVals,
                                    nV,
                                    dM,
                                    matrix,
                                    inverse,
                                    rhs );
              }
              else 
              {
                // should work for 1 or 2 otherwise error 
                DUNE_THROW(InvalidStateException,"Matrix singular in Limiter");
              }

              ++count;
              assert( count <= dim );
            }

            if( ! matrixSingular ) 
            {
              // store linear function
              deoMods.push_back( dM );
              // store vector with points to check 
              comboVec.push_back( check );
            }
          }
#endif
        }

        // Limiting 
        {
          const size_t numFunctions = deoMods.size();
          // for all functions check with all values 
          for(size_t j=0; j<numFunctions; ++j) 
          {
            const std::vector<int> & v = comboVec[j];

            // loop over dimRange 
            for(int r=0; r<dimRange; ++r) 
            {
              double mini = 1.0;
              DomainType& D = deoMods[j][r];
              
              const size_t vSize = v.size();
              for(size_t m=0; m<vSize; ++m)
              {
                const size_t k = v[m];
                const DomainType& omega = barys[k];

                const double g = D * omega;
                const double d = nbVals[k][r];
                const double gd = g * d;

                double m_l = 1.0;
                
                // if product smaller than limit take it as zero 
                if( gd <= 1e-12 )
                {
                  m_l = 0.0;
                }
                else if( (gd > 0.) && (std::abs(g) > std::abs(d)) ) 
                {
                  m_l = (d/g);
                }

                mini = std::min( m_l , mini );
              }

              // scale linear function 
              D *= mini;
            }
          }
          
          // take maximum of limited functions 
          {
            RangeType max (0);
            std::vector< size_t > number(dimRange,0);
            for(int r=0; r<dimRange; ++r) 
            {
              max[r] = deoMods[0][r].two_norm();
            }

            for(size_t l=1; l<numFunctions; ++l) 
            {
              for(int r=0; r<dimRange; ++r) 
              {
                double D_abs = deoMods[l][r].two_norm();
                if( D_abs > max[r] ) 
                {
                  number[r] = l;
                  max[r] = D_abs;
                }
              }
            }
            
            for(int r=0; r<dimRange; ++r) 
            {
              deoMod[r] = deoMods[ number[r] ][r];
            }
          }

          // L2 Projection 
          {
            for(int r=0; r<dimRange; ++r) 
            {
              if( limit[r] ) 
              {
                {
                  const int dofIdx = dofConversion_.combinedDof(0,r);
                  limitEn[dofIdx] = 0;
                }
              }
            }
            
            RangeType ret, tmp;
            // get quadrature for barycenter 
            VolumeQuadratureType quad( en, 2 );

            typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
            //! Note: BaseFunctions must be ortho-normal!!!!
            const BaseFunctionSetType& baseset = limitEn.baseFunctionSet();
            const int quadNop = quad.nop();
            for(int qP = 0; qP < quadNop ; ++qP)
            {
              // get global coordinate 
              DomainType point = geo.global( quad.point( qP ) );
              point -= enBary;

              // get quadrature weight 
              const double intel = quad.weight(qP);

              // evaluate linear function 
              for(int r=0; r<dimRange; ++r) 
              {
                ret[r] = enVal[r] + (deoMod[r] * point);
              }

              for(int r=0; r<dimRange; ++r) 
              {
                if( limit[r] )
                {
                  // only take linear functions components 
                  for(int i=0; i<dim+1; ++i)
                  {
                    const int dofIdx = dofConversion_.combinedDof(i,r);
                    // here evaluateScalar could be used 
                    baseset.evaluate(dofIdx, quad[qP] , tmp);
                    limitEn[dofIdx] += intel * (ret * tmp) ;
                  }
                }
              }
            }
          }
        }
      } //end if limiter 
    }
    
  private:
    // matrix = A^T * A 
    template <class NewMatrixType, class MatrixType> 
    void multiply_AT_A(const NewMatrixType& A, MatrixType& matrix) const
    {
      assert( (int) MatrixType :: rows == (int) NewMatrixType :: cols );
      typedef typename MatrixType :: field_type value_type;

      for(int row=0; row< MatrixType :: rows; ++row)
      {
        for(int col=0; col< MatrixType :: cols; ++col)
        {
          value_type sum = 0;
          for(int k=0; k<NewMatrixType :: rows;  ++k)
          {
            sum += A[k][row] * A[k][col];
          }
          matrix[row][col] = sum;
        }
      }
    }

    template <int newDim, 
              class BaryVectorType, class NbValVectorType,
              class FunctionType, class MatrixType, class VectorType> 
    bool applyLeastSquare(const BaryVectorType& barys,
                          const NbValVectorType& nbVals,
                          const std::vector<int> &nV,
                          FunctionType& dM, 
                          MatrixType& matrix, 
                          MatrixType& inverse, 
                          VectorType& rhs) const 
    {
      enum { dim = dimension };

      // apply least square by adding another point 
      // this should make the linear system solvable 
       
      // need new matrix type containing one row more  
      typedef FieldMatrix<double, newDim , dim> NewMatrixType;
      typedef FieldVector<double, newDim > NewVectorType;

      // new matrix 
      NewMatrixType A ;

      // create matrix 
      for(int k=0; k<newDim; ++k) 
      {
        A[k] = barys[ nV[k] ];
      }

      // matrix = A^T * A 
      multiply_AT_A(A, matrix);

      // invert matrix 
      const double det = FMatrixHelp :: invertMatrix(matrix,inverse);

      if( std::abs( det ) > 0 )
      {
        // now matrix has to be invertable 
        //assert( std::abs( det ) > 0 );
        
        // need new right hand side 
        NewVectorType newRhs;
        
        // calculate D
        for(int r=0; r<dimRange; ++r)
        {
          // get right hand side 
          for(int i=0; i<newDim; ++i) 
          {
            newRhs[i] = nbVals[ nV[i] ][r];
          }
          
          // convert newRhs to matrix 
          rhs = 0;
          A.umtv(newRhs, rhs);
          
          DomainType& D = dM[r];

          // set to zero 
          D = 0;
          // get solution 
          inverse.umv( rhs, D );

        }
        
        // return false because matrix is invertable 
        return false;
      }

      // matrix is still singular 
      return true;
    }
    
    // setup set storing combinations of linear functions 
    void setupComboSet(const int neighbors) const 
    {
      if( spc_.multipleGeometryTypes() || (comboSet_.size() == 0) )
      {
        // clear set 
        comboSet_.clear();

        // maximal number of neighbors 
        std::vector<int> v(dimension,0);
        std::vector<int> null;

        for(int n=0; n<neighbors; ++n)
        {
          v[0] = n;
          if( dimension > 1 )
          {
            for(int j=n+1; j<neighbors; ++j) 
            {
              v[1] = j;
              if( dimension > 2 )
              {
                for(int l=j+1; l<neighbors; ++l) 
                {
                  v[2] = l;
                  comboSet_.insert( VectorCompType( v, null ) );
                }
              }
              else 
              {
                comboSet_.insert( VectorCompType( v, null ) );
              }
            }
          }
          else 
          {
            comboSet_.insert( VectorCompType( v, null ) );
          }
        }

        typedef typename ComboSetType :: iterator iterator; 
        const iterator endit = comboSet_.end();
        for(iterator it = comboSet_.begin(); it != endit; ++it) 
        {
          const KeyType& v = (*it).first;
          CheckType& check = const_cast<CheckType&> ((*it).second);

          // insert all j that are not equal to one of the used in v 
          for(int j=0; j<neighbors; ++j) 
          {
            bool c = true; 
            for(int i=0; i<dimension; ++i) 
            {
              if( v[i] == (int) j ) c = false;
            }
            if( c ) check.push_back(j);
          }
        }
      }
    }

    template <class QuadratureImp>
    void applyLocalNeighbor(IntersectionIteratorType & nit,
            const QuadratureImp & faceQuadInner,
            const QuadratureImp & faceQuadOuter,
            RangeType& totaljump,
            double& umfang) const
    {
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        double umf = caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,l,
                              jump_, jump2_);

        // if jump is not zero 
        if( umf > 0.0 )
        {
          umf *= faceQuadInner.weight(l);
          jump_ *= umf;
          umfang += umf;
          totaljump += jump_;
        }
      }
    }

    template <class QuadratureImp>
    void applyBoundary(IntersectionIteratorType & nit,
                      const QuadratureImp & faceQuadInner,
                      RangeType& totaljump,
                      double& umfang) const
    {
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        double umf = caller_.boundaryFlux(nit, faceQuadInner,l,jump_);

        // if jump is not zero 
        if( umf > 0.0 )
        {
          umf *= faceQuadInner.weight(l);
          jump_ *= umf;
          totaljump += jump_;
          umfang += umf; 
        }
      }
    }

    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // limitation only necessary if order > 0
      if( spc_.order() > 0 )
      {
        // prepare, i.e. set argument and destination 
        prepare(arg, dest);

        // dod limitation 
        IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) {
          applyLocal(*it);
        }
        // finalize
        finalize(arg, dest);
      }
      else 
      {
        // otherwise just copy 
        const DestinationType& U = *(Element<0>::get(arg));
        dest.assign(U);
      }
    }
    
    // make private 
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);
    
  private:
    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType& gridPart_;
    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
    const int faceQuadOrd_;
    mutable RangeType jump_; 
    mutable RangeType jump2_;
    GradientFlux gradientFlux_;

    mutable std::map< const GeometryType, DomainType > baryCenterMap_; 

    mutable ArtificialDiffusionType* artDiff_;

    typedef DGFEntityKey<int> KeyType;
    typedef std::vector<int> CheckType;
    typedef std::pair< KeyType, CheckType > VectorCompType;
    typedef std::set< VectorCompType > ComboSetType;
    mutable ComboSetType comboSet_;
  }; // end DGLimitPass 

} // end namespace Dune 
#endif
