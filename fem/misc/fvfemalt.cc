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
      typedef typename GridType::LeafIterator LeafIterator;
      typedef typename DiscreteFunctionType::LocalFunctionType 
        LocalFunctionType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 

      const GridType & grid = functionSpace_.getGrid();

      LeafIterator endit = grid.leafend   ( level );
      LeafIterator it    = grid.leafbegin ( level ); 

      FieldVector<RangeFieldType,GridType::dimension> tmp;
    
      FixedOrderQuad < RangeFieldType,
        typename FunctionSpaceType::DomainType , polOrd > quad(*it);
      double volRefelem = 0.0;
      for(int i=0; i<quad.nop(); i++)
        volRefelem += quad.weight(i);
      double l1norm = 0.0;

      LocalFunctionType lf = discFunc.newLocalFunction ();

      for(it; it != endit ; ++it)
        {
          discFunc.localFunction( *it , lf  );
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



  //********************************************************************
  // 
  // Project from constant space to linear space
  //
  //********************************************************************
  class ProjectSpaces
  {
    // remember the distances from barycenter to point 
    Array<double> distances_;
  

  public:

    // project constant values to linear values 
    template <class GridIteratorType, class QuadratureType , 
              class ArgFuncType, class DestFuncType, class FuncSpace>
    void projectLocal (GridIteratorType &it, QuadratureType & baryquad, 
                       ArgFuncType &arglf, DestFuncType &destlf, const FuncSpace & functionSpace)
    {
      enum { dimrange = FuncSpace::DimRange };
      enum { dim = GridIteratorType::dimension };
    
      typedef typename FuncSpace::RangeFieldType RangeFieldType; 

      FieldVector<RangeFieldType,dim> bary = it->geometry().global(baryquad.point(0));
      FieldVector<RangeFieldType,dim> dist; 
   
      // here we know that we have 
      // constant and linear Spaces 
      int locNum = arglf.numberOfDofs();
      for(int i=0; i<it->geometry().corners(); i++)
        {
          dist = it->geometry()[i] - bary;
          double r_1 = (1.0/dist.two_norm());

          for(int l=0; l<locNum; l++)
            {
              int k = functionSpace.mapToGlobal(*it,i*dimrange + l);
              distances_[k] += r_1;
              destlf[i*dimrange + l] += (r_1 * arglf[l]);
            }
        }
    }
  
    // project constant values to linear values 
    template <class ArgFuncType, class DestFuncType>
    void project (const ArgFuncType &arg, DestFuncType &dest, int level)
    {
         
      typedef typename DestFuncType::DofIteratorType DestItType;
      typedef typename DestFuncType::FunctionSpace FuSpace;
      typedef typename FuSpace::GridType GridType;
      enum { dim = GridType::dimension };
      enum { dimrange = FuSpace::DimRange };
      const typename DestFuncType::FunctionSpace
        & functionSpace= dest.getFunctionSpace();
    
      GridType &grid = functionSpace.getGrid(); 
   
      int length = functionSpace.size();
      if( distances_.size() < length ) distances_.resize( length );

      for(int i=0; i<length; i++) distances_[i] = 0.0;
    
      typedef typename GridType::template Traits<0>::LevelIterator LevelIterator;
    
      dest.clear();
    
      LevelIterator it = grid. template lbegin<0>( level);
      LevelIterator endit = grid.template lend<0>(level);
      FixedOrderQuad < typename FuSpace::RangeFieldType, 
        typename FuSpace::DomainType , 1 > baryquad(*it);

      typedef typename ArgFuncType::LocalFunctionType   ArgLFType; 
      typedef typename DestFuncType::LocalFunctionType  DestLFType; 
    
      ArgLFType  arglf   =  const_cast<ArgFuncType &> (arg).newLocalFunction ();
      DestLFType destlf  = dest.newLocalFunction ();
    
      // walktrough grid 
      for( ; it != endit; ++it)
        {
          const_cast<ArgFuncType &> (arg).localFunction  ( *it , arglf );
          dest.localFunction ( *it , destlf);
          projectLocal(it,baryquad,arglf,destlf,functionSpace);
        }
    
      typedef typename DestFuncType::DofIteratorType DofIteratorType;
      DofIteratorType dofit = dest.dbegin( level );
      DofIteratorType enddof = dest.dend ( level ); 
      for( ; dofit != enddof; ++dofit )
        {
          (*dofit) /= distances_[dofit.index()];
        }
      // now we have a piecwise linear function global continous 
    }
  };

  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************
  //************************************************************************

  template <class LinFuncType>
  class ProjectDiscontinous 
  {
    typedef typename LinFuncType::FunctionSpaceType LinFuncSpaceType;
    typedef typename LinFuncType::FunctionSpaceType::GridType GridType;
    typedef typename LinFuncSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename LinFuncSpaceType::DomainType DomainType;

    enum { dim = GridType::dimension };
    enum { dimrange = LinFuncType::FunctionSpaceType::DimRange };

    typedef typename LinFuncType::FunctionSpaceType::RangeFieldType  RangeFieldType;
  
    //! the project to piecewise linear operator
    ProjectSpaces pro_;
    
    Array<double> beta_;

    typedef BaryCenterQuad < typename LinFuncSpaceType::RangeFieldType , 
                             typename LinFuncSpaceType::DomainType, 0 > QuadratureType;
    QuadratureType *quad_;

    LinFuncSpaceType & linSpace_;
    FieldMatrix<double,dim,dimrange> grad_; 

    LinFuncType reCon_;
    int level_;
  
    typedef typename LinFuncType::LocalFunctionType LocalFuncType; 
    LocalFuncType recLf_;
  
    JacobianRangeType tmpGrad_;

    FieldVector<RangeFieldType,dim> bary_;
    FieldVector<RangeFieldType,dim> othBary_;
    FieldVector<RangeFieldType,dim> point_;

    DomainType tmpVec_;

    typedef typename GridType::Traits::IntersectionIterator 
    IntersectionIterator; 

    IntersectionIterator nit;
    IntersectionIterator endnit;

  public:

    //! create linear function for reconstruction 
    ProjectDiscontinous ( LinFuncSpaceType & l ) : 
      linSpace_ ( l ) 
                                                 , reCon_ ( l , l.getGrid().maxlevel() , 0 , true ) 
                                                 , recLf_ ( reCon_.newLocalFunction () )
                                                 , nit  ( l.getGrid().template lbegin<0>(0)->ibegin() ) 
                                                 , endnit ( l.getGrid().template lbegin<0>(0)->iend() )   
    {
      GridType &grid = l.getGrid(); 
      quad_ = new QuadratureType ( *(grid.template lbegin<0>(0)) );
    }
  
    // project constant values to linear values 
    template <class ConstFuncType, class DGFuncType>
    void project (const ConstFuncType &Arg, DGFuncType &dest, int level)
    {
      typedef typename ConstFuncType::FunctionSpaceType ConstFuncSpaceType;
      ConstFuncType & arg = const_cast<ConstFuncType &> (Arg); 
      ConstFuncSpaceType & constSpace = arg.getFunctionSpace();
    
      typedef typename GridType::template Traits<0>::LevelIterator LevelIterator; 
      GridType & grid = constSpace.getGrid();
   
      int length = constSpace.size( );
      if( beta_.size () < length ) beta_.resize ( length );
   
    
      // the level we walk on 
      level_ = level;     
    
      pro_.project( arg , reCon_ , level_ );
   
      LevelIterator it = grid.template lbegin<0>  ( level_ );
      LevelIterator endit = grid.template lend<0> ( level_ );
   

      typedef typename ConstFuncType::LocalFunctionType LocalFuncType;
      LocalFuncType lf  = arg.newLocalFunction ();
      LocalFuncType nlf = arg.newLocalFunction ();

      typedef typename DGFuncType::LocalFunctionType DGLocalFuncType;
      DGLocalFuncType destlf = dest.newLocalFunction ();
    
      for( ; it != endit; ++it)
        { 
          arg.localFunction ( *it , lf );
          dest.localFunction ( *it , destlf );
          doLocalPreGlobal( *it, arg ,lf, nlf, destlf ,constSpace);
        }
    }

  private:

    //! do the limitation an stuff 
    template <class EntityType, class ConstFuncType , class LocalFuncType, 
              class DGLocalFuncType, class ConstFuncSpaceType> 
    void doLocalPreGlobal ( EntityType & en, ConstFuncType & arg, LocalFuncType &lf, LocalFuncType &nlf, DGLocalFuncType &dglf, ConstFuncSpaceType &constSpace )
    {
      bary_ = en.geometry().global(quad_->point(0));
    
      int el;
      for(int l=0; l<lf.numberOfDofs(); l++)
        {
          el = constSpace.mapToGlobal(en,l); 
          beta_[el] = 1.0;
        }
    
      { 
        nit = en.ibegin();
        endnit = en.iend();

        for( ; nit != endnit; ++nit)
          {
    
            if(nit.neighbor())
              {
                othBary_ = nit.outside()->geometry().global(quad_->point(0));
         
                arg.localFunction ( *nit.outside() , nlf );

                calcGrad( en , lf);

                for(int l=0; l<lf.numberOfDofs(); l++)
                  {
                    // note here *it, because we need the actual element, not the
                    // neighbors 
                    el = constSpace.mapToGlobal( en ,l); 
                    beta_[el] = std::min( beta_[el] , limitReCon(bary_,othBary_,lf[l],nlf[l],l));
                  }
              }
            if(nit.boundary())
              {
                for(int l=0; l<lf.numberOfDofs(); l++)
                  {
                    // note here *it, because we need the actual element, not the
                    // neighbors 
                    el = constSpace.mapToGlobal(en,l); 
                    beta_[el] = 0.0;
                  }
              }
          } // end neighbor  
      }
  
      for(int i=0; i<en.geometry().corners(); i++)
        {
          point_ = en.geometry()[i] - bary_;
      
          // stimmt so noch nicht 
          for(int l=0; l<lf.numberOfDofs(); l++)
            {
              el = constSpace.mapToGlobal(en,l); 
              dglf[i*dimrange + l] = evalFunc(point_,beta_[el],lf[l],l); 
            }
        }
    }

    template <class EntityType , class LFType>
    void calcGrad (EntityType &en, LFType &lf )  
    {
      typedef typename LinFuncSpaceType::BaseFunctionSetType BaseFuncSetType;
      typedef typename LinFuncSpaceType::RangeType  DRangeType;
      typedef typename LinFuncSpaceType::DomainType DomainType;

      reCon_.localFunction( en , recLf_ );
    
      FieldMatrix<double,GridType::dimension,GridType::dimension>& inv =
        en.geometry().Jacobian_inverse((*quad_).point(0));

      const BaseFuncSetType &baseSet = linSpace_.getBaseFunctionSet( en );

      // gradient of the function
      for(int l=0; l<lf.numberOfDofs(); l++)
        {
          tmpVec_ = 0.0;
          for(int i=l; i<baseSet.getNumberOfBaseFunctions(); i+=dimrange)
            {
              // calc grad on point 0 because grad is constant for triangles 
              baseSet.jacobian(i,(*quad_),0,tmpGrad_);
    
              for(int j=0; j<dim; j++)
                tmpVec_(j) += recLf_[i] * tmpGrad_(j,0); 
            }
    
          // grad with respect to actual element is gradient on reference
          // element multiplied with transpose of jacobian inverse 
          grad_(l) = inv.mult_t(tmpVec_);
        }
    }
  
    template <class ctype, int dim>
    double limitReCon (FieldVector<ctype,dim> &bary, 
                       FieldVector<ctype,dim> &othBary , double  u , double u_jl , int l)  
    {
      typedef typename LinFuncSpaceType::BaseFunctionSetType BaseFuncSetType;
      typedef typename LinFuncSpaceType::RangeType  DRangeType;

      // note that this only wotks for structured Grid 
      tmpVec_ = othBary - bary;
  
      double g = grad_(l) * tmpVec_;
    
      double d = u_jl - u;

      // product < 0 means gradient point in wrong direction 
    
      if(std::abs(g) > 1.0E-25 )
        {
          d /= g;
      
          if( d <= 0.0)
            return 0.0; 
   
          return std::min(1.0, d );
        }
      return 0.0;
    }
  
    template <class ctype, int dim>
    double evalFunc (FieldVector<ctype,dim> &point, 
                     double beta, double  u , int l )  
    {
      if(std::abs(beta) > 0.0)
        {
          // gradient is calced from before 
          double add = grad_(l) * point;
          return u + beta*add;
        }
      return u;
    }
  }; // end ProjectDiscontinous 


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

    mutable LocalFuncType *oldLf_;
    mutable LocalFuncType *neighLf_;
    mutable LocalFuncType *up_;
    mutable LocalFuncType *upNeigh_;

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
    LocalFuncType enError_; 
    LocalFuncType neighError_; 
  
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
      oldLf_ ( NULL ), neighLf_ (NULL) , up_ (NULL) , upNeigh_ (NULL),
      quad_ ( *(f.begin()) ),
      errorFunc_ ( error ),
      enError_ ( errorFunc_->newLocalFunction () ) ,
      neighError_ ( errorFunc_->newLocalFunction () )
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

      if(!oldLf_)
        oldLf_   = new LocalFuncType ( argTemp.newLocalFunction() );     

      if(!neighLf_)
        neighLf_ = new LocalFuncType ( argTemp.newLocalFunction() );     
  
      if(!up_)
        up_      = new LocalFuncType ( dest.newLocalFunction() );     

      if(!upNeigh_)
        upNeigh_ = new LocalFuncType ( dest.newLocalFunction() );     
    }
 
    // return calculated time step size 
    void finalizeGlobal ()
    {
      dtOld_ = dtMin_;

      if( oldLf_ ) delete oldLf_;     oldLf_ = 0;
      if( neighLf_ ) delete neighLf_; neighLf_ = 0;
      if( up_ ) delete up_;           up_ = 0;
      if( upNeigh_ ) delete upNeigh_; upNeigh_ = 0;

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
    
      oldSol.localFunction ( en , *oldLf_ ); 
      update.localFunction ( en , *up_    ); 

      errorFunc_->localFunction( en, enError_ );

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
              
          // if intersection is with neighbor and neighbors number is bigger
          // then ours ==> calculate flux 
          if( nit.neighbor() ) {

            if( nit.outside()->globalIndex() > en.globalIndex()
                || nit.outside()->partitionType() == GhostEntity) {

              // edge or face volume
                         
              // volume of neighbor 
              nvol = 
                volRefelem_ * nit.outside()->geometry().integrationElement(centerPoint_);
              double nvol_1 = 1/nvol;
          
              // get access to local functions  
              oldSol.localFunction ( *nit.outside() , *neighLf_ );
              update.localFunction ( *nit.outside() , *upNeigh_ );
        
              // evaluate the local function on the middle point of the actual face
              // * Why numberInSelf, numberInNeighbor here
              //oldLf_->evaluate(en, quad_, nit.numberInSelf(), valEl_);
              //neighLf_->evaluate(*nit.outside(), quad_, nit.numberInNeighbor(), valNeigh_);

              oldLf_->evaluate(en, quad_, 0, valEl_);
              neighLf_->evaluate(*nit.outside(), quad_, 0, valNeigh_);
        
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
                errorFunc_->localFunction(*nit.outside(), neighError_);
                if (enError_[0] < eta) enError_[0] = eta;
                if (neighError_[0] < eta) neighError_[0] = eta;
              }

              // add local flux to update function 
              for(int l=0; l< oldLf_->numDofs(); l++) {
                // * bugfix? multiply with h as well
                (*up_)[l]      += vol_1 * flux_[l] * h;
                (*upNeigh_)[l] -= nvol_1  * flux_[l] * h;
              }
        
            }
          }
          // if intersection is with boundary 
          if(nit.boundary()) {
            //std::cout << "At boundary\n";
            // evaluate boundary function
        
            typedef typename BoundaryManagerType::BoundaryType BoundaryType;
            const BoundaryType& bc = 
              bcManager_.getBoundaryCondition(nit.boundaryEntity().id());
  
            midPoint = nit.intersectionGlobal().global(mid);

            // evaluate the local function on the middle point of the actual face
            //oldLf_->evaluate   ( en , quad_ , nit.numberInSelf() , valEl_ );
            oldLf_->evaluate   ( en , quad_ , 0 , valEl_ );

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
            for(int l=0; l<oldLf_->numDofs(); l++)
              {
                (*up_)[ l ] += flux_[l];
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

