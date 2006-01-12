#ifndef DUNE_FVFEMALT_CC
#define DUNE_FVFEMALT_CC

#define MYABS(a) (( a > 0.0 ) ? (a) : -(a)) 
#define INFTY 1.0/1.0E-16

#include <dune/quadrature/fixedorder.hh>
#include <dune/quadrature/barycenter.hh>
#include <dune/fem/common/boundary.hh>


namespace Dune 
{

  static double frac[11] = 
    {INFTY,
     1.0,
     0.5,
     1.0/3.0,
     0.25,
     0.2,
     1.0/6.0,
     1.0/7.0,
     0.125,
     1.0/9.0,
     0.1
    };

  inline static double fraction ( int num )
  {
    assert(num >= 0);
    assert(num <= 10);

    return frac[num];
  }

  
  class L1Norm
  {
  public:
    template <int polOrd, class DiscreteFunctionType>
    double norm (int level, DiscreteFunctionType &discFunc)
    {
      typedef typename DiscreteFunctionType::FunctionSpace FunctionSpaceType;
    
      const typename DiscreteFunctionType::FunctionSpace
        & functionSpace_= discFunc.getFunctionSpace();

      typedef typename FunctionSpaceType::GridType GridType;
      typedef typename FunctionSpaceType::IteratorType IteratorType;
      typedef typename DiscreteFunctionType::LocalFunctionType 
        LocalFunctionType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 

      const GridType & grid = functionSpace_.getGrid();

      IteratorType endit = functionSpace_.end();
      IteratorType it    = functionSpace_.begin(); 

      FieldVector<RangeFieldType,GridType::dimension> tmp;
    
      FixedOrderQuad < RangeFieldType,
        typename FunctionSpaceType::DomainType , polOrd > quad(*it);
      double volRefelem = 0.0;
      for(int i=0; i<quad.nop(); i++)
        volRefelem += quad.weight(i);
      double l1norm = 0.0;

      for(it; it != endit ; ++it)
        {
          LocalFunctionType lf = discFunc.localFunction( *it );
          double vol = volRefelem * (*it).geometry().integrationElement(tmp);
          l1norm += MYABS(vol * lf[0]);
        }
      return l1norm;
    }

  };

  class L1Error
  {
  public:
    template <int polOrd, class FunctionType, class DiscreteFunctionType>
    double error (int level, FunctionType &f, DiscreteFunctionType &discFunc, double time)
    {
      typedef typename DiscreteFunctionType::FunctionSpace FunctionSpaceType;
    
      const typename DiscreteFunctionType::FunctionSpace
        & functionSpace_= discFunc.getFunctionSpace();

      typedef typename FunctionSpaceType::GridType GridType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
      typedef typename GridType::template Traits<0>::LevelIterator LevelIterator;
      typedef typename DiscreteFunctionType::LocalFunctionType 
        LocalFunctionType;

      enum { dim = GridType::dimension };

      GridType & grid = functionSpace_.getGrid();

      LevelIterator endit = grid.template lend<0>   ( level );
      LevelIterator it    = grid.template lbegin<0> ( level ); 
  
      FieldVector<RangeFieldType,dim> tmp;
    
      FixedOrderQuad < RangeFieldType ,    
        typename FunctionSpaceType::DomainType , polOrd > quad(*it);
    
      double error = 0.0;
      typename FunctionSpaceType::RangeType ret (0.0);
      typename FunctionSpaceType::RangeType phi (0.0);

      FieldVector<RangeFieldType,dim> point;
      LocalFunctionType lf = discFunc.newLocalFunction();

      double intWeight;
    
      double val = 0.0;
      for( ; it != endit ; ++it)
        {
          discFunc.localFunction ( *it , lf );

          val = 0.0;
          for(int qP = 0; qP < quad.nop(); qP++)
            {
              point = (*it).geometry().global(quad.point(qP));
              intWeight = quad.weight(qP) * 
                (*it).geometry().integrationElement(quad.point(qP));
              f.evaluate(point,time,ret);
        
              lf.evaluate(*it,point,phi);
              val += intWeight * std::abs((ret(0) - phi(0)));
            }
          error += std::abs(val);
        }
      return error;
    }

  };


  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //**************************************************************************
  //
  //  --ScalarFV
  //
  //**************************************************************************
  template <class NumericalFluxFunction, class DiscFuncType>
  class ScalarFV : 
    public LocalOperatorDefault <DiscFuncType,
                                 DiscFuncType, 
                                 typename DiscFuncType::RangeFieldType , 
                                 ScalarFV < NumericalFluxFunction , DiscFuncType > >
  {
    typedef ScalarFV <NumericalFluxFunction, DiscFuncType > MyType;  

    typedef DiscFuncType DFType;
    typedef typename DiscFuncType::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::GridType GridType; 
    typedef typename DiscFuncType::RangeFieldType RangeFieldType;

    typedef typename NumericalFluxFunction::NormalType NormalType;
    typedef typename NumericalFluxFunction::ValueType ValueType;
    typedef BoundaryManager<FunctionSpaceType> BoundaryManagerType;

    enum { dim = GridType::dimension };
 
    //! the numerical flux function  
    const NumericalFluxFunction fluxFcn_;
    
    //! Boundary manager
    BoundaryManagerType bcManager_;

    //! function space of grid function 
    FunctionSpaceType& fSpace_;
 
    //! true if prepare and finalize should be done
    bool doPrepare_;
  
    //! true if simulation is adaptive
    bool adaptive_;

    //! level of interest 
    // int level_;

    //! timestep size 
    mutable double dtMin_;
 
    //! timestep size 
    mutable double dtOld_;
 
    typedef typename DiscFuncType::FunctionSpaceType::RangeType  DRangeType;
    typedef typename DiscFuncType::LocalFunctionType LocalFuncType;
    typedef RangeFieldType Real;
    // help variables 
    mutable FieldVector<RangeFieldType,dim-1> mid;
    mutable FieldVector<RangeFieldType,dim> midPoint;
    mutable FieldVector<RangeFieldType,dim> centerPoint_;

    // local storage of evaluation of the grid functions 
    mutable DRangeType valEl_;
    mutable DRangeType valNeigh_;

    // local storage of numerical fluxes 
    mutable DRangeType flux_;
    mutable DRangeType fluxEl_;

    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<1>::Entity FaceType;

    //! quadrature for evaluation 
    //! the quadrature points should be the mid points of the faces 
    typedef BaryCenterQuad < typename FunctionSpaceType::RangeFieldType ,
                             typename FunctionSpaceType::DomainType , 0 > QuadratureType;
    QuadratureType quad_;
  
    double volRefelem_;

    const DiscFuncType * arg_;
    DiscFuncType * dest_;

    DiscFuncType* errorFunc_;

    typedef typename GridType::Traits::LocalIdSet LocalIdSetType; 
    const LocalIdSetType & localIdSet_; 
  
  public:
    // Constructor 
    ScalarFV ( const NumericalFluxFunction& fluxFcn,
               const BoundaryManagerType& bc,
               FunctionSpaceType &f,
               DiscFuncType* error,
               bool adaptive,
               bool doPrepare) :
      fluxFcn_(fluxFcn) ,  
      bcManager_(bc),
      fSpace_ ( f ),
      doPrepare_ ( doPrepare ), 
      adaptive_(adaptive), 
      mid(0.5) ,
      quad_ ( *(f.begin()) ),
      errorFunc_ ( error ),
      localIdSet_( f.grid().localIdSet() )  
    { 
      //level_ = f.getGrid().maxlevel();
      typedef FieldVector<RangeFieldType,dim-1> MidType;
      typedef FixedOrderQuad < typename FunctionSpaceType::RangeFieldType,MidType,1> MidQuadType;
      IntersectionIterator nit = (f.begin())->ibegin();
      GeometryType eltype = nit.intersectionGlobal().type();
      MidQuadType midquad( eltype );
    
      mid = midquad.point(0);

      typedef FieldVector<RangeFieldType,dim> CenterType;
      typedef FixedOrderQuad <RangeFieldType, CenterType, 1> CenterQuadType;
      CenterQuadType centerQuad(*(f.begin()));
      centerPoint_ = centerQuad.point(0);

      volRefelem_ = 0.0;
      for(int qp=0; qp < quad_.nop(); qp++ )
        volRefelem_ += quad_.weight(qp);

      std::cout << "Constructor ScalarFV finished! \n";
    }

    // Destructor 
    ~ScalarFV () 
    {
      if( oldLf_ ) delete oldLf_; 
      if( neighLf_ ) delete neighLf_; 
      if( up_ ) delete up_;
      if( upNeigh_ ) delete upNeigh_;
    }

    void setTimeStep(double dt) {
      dtOld_ = dt;
    }

    double getTimeStep() const {
      return dtMin_;
    }

    // set update function to zero 
    void prepareGlobal(const DiscFuncType &arg, DiscFuncType &dest)
    {
      arg_  = &arg; 
      dest_ = &dest;    
      assert(arg_ != NULL); assert(dest_ != NULL);
      DiscFuncType & argTemp = const_cast<DiscFuncType &> (*arg_);

      dtMin_  = 1.0E+308;
      if(dtOld_ > 1.0 ) dtOld_ = 1.0;

      errorFunc_->clear();
 
      if(doPrepare_)
        {
          dest.clear();
        }
    }
 
    // return calculated time step size 
    void finalizeGlobal ()
    {
      dtOld_ = dtMin_;
    }

    // apply operator locally
    template <class EntityType>
    void applyLocal (EntityType &en)
    {
      DiscFuncType & argTemp = const_cast<DiscFuncType &> (*arg_);
      onCell( en , argTemp,  *dest_);
    }
  
  private:
    template <class EntityType>
    void onCell(EntityType &en, DiscFuncType &oldSol, DiscFuncType &update) 
    {
      //********************************************************************  
      // necessary typedefs 
      typedef typename GridType::Traits::
        IntersectionIterator NeighIt;
      //********************************************************************  
    
      LocalFuncType oldLf = oldSol.localFunction ( en ); 
      LocalFuncType up    = update.localFunction ( en ); 

      LocalFuncType enError = errorFunc_->localFunction( en );

      // get volume of actual element 
      double vol = 
        (volRefelem_ * en.geometry().integrationElement(centerPoint_)); 
      double vol_1 = 1.0/vol;

      double dtLocal;
      double nvol;

      NeighIt endnit = en.iend();
      for(NeighIt nit = en.ibegin(); nit != endnit; ++nit)
        {
          NormalType normal = nit.integrationOuterNormal(mid);
          double h = normal.two_norm();
          normal *= 1.0/h;
              
          if( nit.neighbor() ) {
            typedef EntityType :: EntityPointer EntityPointerType;
            EntityPointerType ep = nit.outside();
            EntityType & neighbour = *ep;

            // if intersection is with neighbor and neighbors number is bigger
            // then ours ==> calculate flux 
            if( (localIdSet_.id(neighbour) > localIdSet_.id(en) )
                || neighbour.partitionType() == GhostEntity) 
            {

              // edge or face volume
                         
              // volume of neighbor 
              nvol = 
                volRefelem_ * neighbour.geometry().integrationElement(centerPoint_);
              double nvol_1 = 1/nvol;
          
              // get access to local functions  
              LocalFuncType neighLf = oldSol.localFunction ( neighbour );
              LocalFuncType upNeigh = update.localFunction ( neighbour );
        
              // evaluate the local function on the middle point of the actual face
              // * Why numberInSelf, numberInNeighbor here
              //oldLf_->evaluate(en, quad_, nit.numberInSelf(), valEl_);
              //neighLf_->evaluate(*nit.outside(), quad_, nit.numberInNeighbor(), valNeigh_);

              oldLf.evaluate(en, quad_, 0, valEl_);
              neighLf.evaluate(neighbour , quad_, 0, valNeigh_);
        
              // calculate numerical flux , g(u,v)
              double dtLocal = fluxFcn_(valEl_, valNeigh_, normal, normal, flux_);
              dtLocal = 
                (dtLocal < std::numeric_limits<double>::min()) ? dtMin_ :
                3.0*std::min(vol, nvol)/(dtLocal*h);
              if(dtLocal < dtMin_) dtMin_ = dtLocal;
              //std::cout << "dtLocal: " << dtLocal << std::endl;
        
              if (adaptive_) {
                // Gradient estimator for density (first component)
                //FieldVector<Real, dim> dist = en.geometry().global(centerPoint_);
                //dist -= nit.outside()->geometry().global(centerPoint_);
                // * zeroth element is the density
                /* Original formulation of Gessner
                   Real eta = std::fabs(valEl_[0] - valNeigh_[0]) / dist.two_norm();
                */
                // * Formulation by Dedner
                Real eta = 2.0 * std::fabs(valEl_[0] - valNeigh_[0]) / (valEl_[0] + valNeigh_[0]);
                //std::cout << "Eta = " << eta << std::endl;
                LocalFuncType neighError = errorFunc_->localFunction(neighbour);
                if (enError[0] < eta) enError[0] = eta;
                if (neighError[0] < eta) neighError[0] = eta;
              }

              // add local flux to update function 
              for(int l=0; l< oldLf.numDofs(); l++) {
                // * bugfix? multiply with h as well
                up[l]      += vol_1 * flux_[l] * h;
                upNeigh[l] -= nvol_1  * flux_[l] * h;
              }
        
            }
          }
          // if intersection is with boundary 
          if(nit.boundary()) {
            //std::cout << "At boundary\n";
            // evaluate boundary function
        
            typedef typename BoundaryManagerType::BoundaryType BoundaryType;
            const BoundaryType& bc = 
              bcManager_.getBoundaryCondition(nit.boundaryId());
  
            midPoint = nit.intersectionGlobal().global(mid);

            // evaluate the local function on the middle point of the actual face
            //oldLf_->evaluate   ( en , quad_ , nit.numberInSelf() , valEl_ );
            oldLf.evaluate   ( en , quad_ , 0 , valEl_ );

            if (bc.boundaryType() == BoundaryType::Dirichlet) {
           
              // evaluate boundary on the middle point of the actual face
              ValueType bval;
              bc.evaluate(valEl_, midPoint, normal, bval);
           
              // calculate numerical flux 
              dtLocal = 
                fluxFcn_(valEl_, bval, normal, normal, flux_);
              dtLocal =   
                (dtLocal < std::numeric_limits<double>::min()) ? dtMin_ :
                3.0*vol/(dtLocal*h);
              if (dtLocal < dtMin_) dtMin_ = dtLocal;
              //std::cout << "dtLocal: " << dtLocal << std::endl;

            } else if (bc.boundaryType() == BoundaryType::Neumann) {
              //    std::cout << "n: " << nit.integrationOuterNormal(mid) << ",x: " << midPoint <<std::endl;
              bc.evaluate(valEl_, midPoint, normal, flux_);
            }
            // add local flux to update function 
            // * bugfix? multiply with area as well...
            flux_ *= vol_1 * h;
            for(int l=0; l<oldLf.numDofs(); l++)
              {
                up[ l ] += flux_[l];
              }
        
          } // end Boundary 
      
        } // end Neighbor 

      // calculate local time step size 
      //dtSum = vol/dtSum;
      //std::cout << "dtSum: " << dtSum << std::endl;
      //if( dtSum < dtMin_) dtMin_ = dtSum;

    } // end onCell 

  }; // end class ScalarFV 

} // end namespace Dune

#endif

