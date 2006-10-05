#ifndef APOSTERIORI_HH 
#define APOSTERIORI_HH 

#include <fem/quadrature/cachequad.hh>
#include <fem/quadrature/gausspoints.hh>
#include <fem/pass/tuples.hh>

namespace Dune {
template <class TupType>
struct TupleToPair {
  typedef Pair<typename TupType::Type1,
	       typename TupleToPair<typename TupType::Type2>::ReturnType> ReturnType;
  static ReturnType convert(TupType tup) {
    return ReturnType(tup.first(),
		      TupleToPair<typename TupType::Type2>::convert(tup.second()));
  }
};
template <class T>
struct TupleToPair<Pair<T,Nil> > {
  typedef Pair<T,Nil> ReturnType;
  static ReturnType convert(Pair<T,Nil> tup) {
    return ReturnType(tup.first(),Nil());
  }
};
  template <class DiscModelType,class DiscreteFunctionType> 
  class LocalResiduum {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionType::RangeType RangeType;
    typedef typename DiscreteFunctionType::DomainType DomainType;
    typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    enum { dimR = RangeType :: dimension };
    enum { dimD = DomainType :: dimension };
    enum {nopTime = 2}; 
  public:  
    template <class EntityType>
    RangeType element (const DiscModelType& model,DiscreteFunctionType &discFunc,
		       const EntityType &en,
		       double time,
		       double dt,
		       int maxp,
		       RangeType& infProj) {
      typedef typename FunctionSpaceType::GridType GridType;
      const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
      int quadOrd = 2 * space.polynomOrder() + 2;
      const GaussPts& timeQuad = GaussPts::instance();
      discFunc.setEntity(en);
      CachingQuadrature <GridType , 0 > quad(en,quadOrd);
      RangeType u    (0.0);
      RangeType dtu  (0.0);
      RangeType divfu(0.0);
      JacobianRangeType dxu  (0.0);
      RangeType error(0.0);
      std::vector<RangeType> U(quad.nop()*nopTime);
      FieldVector<RangeType,nopTime> average(0);
      double vol = en.geometry().volume();
      for(int qP = 0; qP < quad.nop(); qP++) {
	double det = quad.weight(qP) 
	  * en.geometry().integrationElement(quad.point(qP));
	for (int qT=0;qT<nopTime;++qT) {
	  double s = timeQuad.point(nopTime,qT);
	  double ldet = det*dt*timeQuad.weight(nopTime,qT);
	  u   = discFunc.uval(en,quad,qP,s,maxp);
	  dxu = discFunc.dxuval(en,quad,qP,s,maxp);
	  dtu = discFunc.dtuval(en,quad,qP,s,maxp);
	  model.model().jacobian(en,time+s*dt,quad.point(qP),u,dxu,divfu);
	  for(int k=0; k<dimR; ++k) {
	    error[k] += ldet * fabs(dtu[k] + divfu[k]);
	    average[qT][k] += det/vol * u[k];
	  }
	  U[qP*nopTime+qT]=u;
	}
      }
      infProj = 0.0;
      for(int qP = 0; qP < quad.nop(); qP++) {
	for (int qT=0; qT<nopTime;qT++) {
	  for(int k=0; k<dimR; ++k) {
	    double linfProj = fabs(U[qP*nopTime+qT][k]-average[qT][k]);
	    infProj[k] = (infProj[k]<linfProj)?linfProj:infProj[k];
	  }
	}
      }
      return error;
    }
    template <class IIteratorType>
    RangeType jump (const DiscModelType& model,DiscreteFunctionType &discFunc,
		    IIteratorType &nit,
		    double time,
		    double dt,
		    int maxpen,int maxpnb) {
      typedef typename IIteratorType::Entity EntityType;
      typedef typename EntityType :: EntityPointer EntityPointerType;
      typedef typename FunctionSpaceType::GridType GridType;
      const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
      const GaussPts& timeQuad = GaussPts::instance();
      TwistUtility<GridType> twistUtil(space.grid());
      int quadOrd = 2 * space.polynomOrder() + 2;
      RangeType error(0.0);
      EntityPointerType np = nit.outside();
      EntityType & nb = *np;
      EntityPointerType ep = nit.inside();
      EntityType & en = *ep;
      int twistSelf = twistUtil.twistInSelf(nit); 
      CachingQuadrature<GridType,1> 
	faceQuadInner(nit, quadOrd, twistSelf, 
		      CachingQuadrature<GridType,1>::INSIDE);
      int faceQuadInner_nop = faceQuadInner.nop();
      int twistNeighbor = twistUtil.twistInNeighbor(nit);
      CachingQuadrature<GridType,1> 
	faceQuadOuter(nit, quadOrd, twistNeighbor,
		      CachingQuadrature<GridType,1>::OUTSIDE);
      RangeType uLeft(0.0),uRight(0.0);
      RangeType gLeft(0.0),gRight(0.0);
      RangeType normalFLeft(0.0),normalFRight(0.0);
      JacobianRangeType fLeft(0.0),fRight(0.0);
      for (int l = 0; l < faceQuadInner_nop; ++l) {
	// double det=nit.intersectionGlobal().
	//	integrationElement(faceQuadInner.localPoint(l));
	DomainType normal = 
	  nit.integrationOuterNormal(faceQuadInner.localPoint(l));
	for (int qT=0;qT<nopTime;++qT) {
	  double s = timeQuad.point(nopTime,qT);	
	  double ldet = dt*timeQuad.weight(nopTime,qT);
	  discFunc.setEntity(en);
	  uLeft = discFunc.uval(en,faceQuadInner,l,s,maxpen);
	  discFunc.setEntity(nb);
	  uRight = discFunc.uval(nb,faceQuadOuter,l,s,maxpnb);	
	  model.numericalFlux(nit,time+dt*s,faceQuadInner.localPoint(l),
			      uLeft,uRight,gLeft,gRight);
	  model.model().analyticalFlux(en,time+dt*s,faceQuadInner.point(l),
				       uLeft,fLeft);
	  model.model().analyticalFlux(nb,time+dt*s,faceQuadOuter.point(l),
				       uRight,fRight);
	  normalFLeft = 0.0;
	  normalFRight = 0.0;
	  fLeft.umv(normal,normalFLeft);
	  fRight.umv(normal,normalFRight);
	  for(int k=0; k<dimR; ++k) {
	    error[k] += ldet * 
	      fabs(gLeft[k]+gRight[k]-normalFLeft[k]-normalFRight[k]);
	  }
	}
      }
      return error;
    }
    template <class EntityType>
    int computePolDeg(const DiscModelType& model,
		      DiscreteFunctionType &discFunc,
		      const EntityType& en,
		      int maxp,
		      double lambda,
		      double& projErr,
		      bool doPAdapt,
		      bool start) {
      discFunc.setEntity(en);
      typedef typename FunctionSpaceType::GridType GridType;
      const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
      int polOrd = space.polynomOrder();
      int quadOrd = 2 * space.polynomOrder() + 2;
      projErr = 0.;
      if (!doPAdapt) {
	return polOrd;
      }
      CachingQuadrature <GridType , 0 > quad(en,quadOrd);
      RangeType maxOrd(polOrd); 
      for (int r=0;r<dimR;r++) {
	for(int qP = 0; qP < quad.nop(); qP++) {
	  RangeType average = discFunc.uval(en,quad,qP,1.,0);
	  int p = int(maxOrd[r]);
	  for (; p>0;--p) { 
	    RangeType u = discFunc.uval(en,quad,qP,1.,p);
	    u -= average;
	    if (fabs(u[r]) < lambda) {
	      break;
	    } 
	  }
	  maxOrd[r] = p;
	}
      }
      for(int qP = 0; qP < quad.nop(); qP++) {
	double det = quad.weight(qP) 
	  * en.geometry().integrationElement(quad.point(qP));
	RangeType proj(0);
	if (start) 
	  model.model().problem().
	    evaluate(en.geometry().global(quad.point(qP)),proj);
	else
	  proj = discFunc.uval(en,quad,qP,1.0,maxp);
	if (start)
	  proj -= discFunc.uval(en,quad,qP,0.0,int(maxOrd[0]));
	else
	  proj -= discFunc.uval(en,quad,qP,1.0,int(maxOrd[0]));
	projErr += proj.one_norm()*det; 
      }
      return int(maxOrd[0]);
    }
  };
  template <class GridPartType,class DiscModelType,class DiscreteFunctionType> 
  class Residuum {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionType::RangeType RangeType;
    typedef typename DiscreteFunctionType::DomainType DomainType;
    typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionType::DestinationType DiscFuncType;
    typedef typename GridPartType::GridType GridType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    enum { dimR = RangeType :: dimension };
    enum { dimD = DomainType :: dimension };
    // Result Functions
    // typedef LeafGridPart<GridType> GridPartType;
    //typedef DGAdaptiveLeafIndexSet<GridType> DGIndexSetType;
    // typedef typename LeafGridPart<GridType>::IndexSetType DGIndexSetType;
    //typedef DefaultGridPart<GridType,DGIndexSetType> GridPartType;

    typedef FunctionSpace < double , double, dimD , 1 > ScalarFSType;
    typedef DiscontinuousGalerkinSpace<ScalarFSType, GridPartType, 0> 
    ConstDiscSType;
    typedef DFAdapt<ConstDiscSType> 
    ConstDiscFSType;
    typedef typename ConstDiscFSType::LocalFunctionType LConstDiscFSType;
    bool doPAdapt_;
  public:  
    typedef ConstDiscFSType DestinationType;
    // DGIndexSetType iset_;
    FieldVector<int,10> padapt_num;
    GridPartType& part_;
    ConstDiscSType space_;
    DestinationType ind_;
    DestinationType RT_,RS_,RP_,rho_,lambda_,maxPol_,maxPolNew_;
    typedef Tuple<DestinationType*,DestinationType*,DestinationType*,
		  DestinationType*,DestinationType*> OutputTType;
    typedef typename TupleToPair<OutputTType>::ReturnType OutputType;
    OutputType output() {
      return TupleToPair<OutputTType>::convert
	(OutputTType(&RT_,&RS_,&rho_,&lambda_,&maxPol_));
    }
    Residuum(GridPartType& part,int maxOrd,
	     bool doPAdapt=true) : 
      doPAdapt_(doPAdapt),
      // iset_(grid),
      // part_(grid,iset_),
      part_(part),
      space_(part_),
      ind_("Indicator",space_),
      RT_("Element Res.",space_),
      RS_("Jump Res.",space_),
      RP_("Proj Res.",space_),
      rho_("rho",space_),
      lambda_("lambda",space_),
      maxPol_("PolDeg",space_),
      maxPolNew_("PolDeg",space_),
      padapt_num(0) {
	typedef typename FunctionSpaceType::IteratorType IteratorType;
	IteratorType endit = space_.end();
	for(IteratorType it = space_.begin(); 
	    it != endit ; ++it) {
	  LConstDiscFSType lmaxPol = maxPol_.localFunction(*it);
	  LConstDiscFSType lmaxPolNew = maxPolNew_.localFunction(*it);
	  lmaxPol[0] = maxOrd;
	  lmaxPolNew[0] = maxOrd;
	}
      }
    template <class EntityType>
    int usePolDeg(const EntityType& en) {
      LConstDiscFSType lmaxPol = maxPol_.localFunction(en);
      return int(lmaxPol[0]);
    }
    void reset() {
      maxPol_.assign(maxPolNew_);
    }
    int numPAdapt() {
      return padapt_num[0]+padapt_num[1];
    }
    template <class Adapt>
    RangeType calc (const DiscModelType& model,
		    Adapt& adapt,
		    DiscreteFunctionType &discFunc,
		    double time,double dt,
		    bool start = false) {
      typedef typename FunctionSpaceType::IteratorType IteratorType;
      const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
      int polOrd = space.polynomOrder();
      if (adapt)
	adapt->clearIndicator();
      RangeType ret(0.);
      LocalResiduum<DiscModelType,DiscreteFunctionType> localRes;
      IteratorType endit = space.end();
      RS_.clear();
      RT_.clear();
      padapt_num = 0;
      for(IteratorType it = space.begin(); 
	  it != endit ; ++it) {
	LConstDiscFSType lmaxPol = maxPol_.localFunction(*it);
	LConstDiscFSType lRT = RT_.localFunction(*it);
	LConstDiscFSType lrho = rho_.localFunction(*it);
	double h = sqrt(it->geometry().volume());
	RangeType infProj;
	RangeType resid = 
	  localRes.element(model,discFunc,*it,time,dt,int(lmaxPol[0]),infProj);
	lRT[0] = resid.one_norm();
	lrho[0] = h + infProj.one_norm();
	assert(lrho[0] == lrho[0]);
	assert(lrho[0] < 1e10);
	resid *= lrho[0];
	ret += resid;

	IntersectionIterator endnit = it->ileafend();
	IntersectionIterator nit = it->ileafbegin();
      
	double maxh = lrho[0]; 
	for (; nit != endnit; ++nit) {
	  if (nit.neighbor()) {
	    if (part_.indexSet().index(*(nit.outside())) > 
		part_.indexSet().index(*it)) {
	      LConstDiscFSType lmaxPolnb = 
		maxPol_.localFunction(*(nit.outside()));
	      RangeType resjmp = 
		localRes.jump(model,discFunc,nit,time,dt,
			      int(lmaxPol[0]),int(lmaxPolnb[0]));
	      resjmp *= maxh;
	      double jump = resjmp.one_norm();
	      {
		LConstDiscFSType lRS = RS_.localFunction(*it);
		lRS[0] += jump*0.5;
	      }
	      {
		LConstDiscFSType lRS = RS_.localFunction(*(nit.outside()));
		lRS[0] += jump*0.5;
	      }
	      ret += resjmp;
	    } // end if ...
	  } // end if neighbor
	  if (nit.boundary()) {
	  } // end if boundary
	}
  }
      ret = 0;
      double pot = 2.*(polOrd+1); // double(polOrd+2)/double(polOrd+1);
      for(IteratorType it = space.begin(); 
	  it != endit ; ++it) {
	double vol = it->geometry().volume();
	double h = sqrt(vol); // 0.5*sqrt(it->geometry().integrationElement(DomainType(0))); // 0.001
	LConstDiscFSType lRT = RT_.localFunction(*it);
	LConstDiscFSType lRS = RS_.localFunction(*it);
	LConstDiscFSType lRP = RP_.localFunction(*it);
	LConstDiscFSType lrho = rho_.localFunction(*it);
	LConstDiscFSType llam = lambda_.localFunction(*it);
	LConstDiscFSType lmaxPol = maxPol_.localFunction(*it);
	LConstDiscFSType lmaxPolNew = maxPolNew_.localFunction(*it);
	llam[0] = h / pow(h + h*(lRT[0]+lRS[0])/(vol*dt),pot);
	lmaxPolNew[0] = 
	  localRes.computePolDeg(model,
				 discFunc,*it,int(lmaxPol[0]),
				 llam[0],lRP[0],doPAdapt_,start);
	if (lmaxPol[0]<lmaxPolNew[0]) ;
	   lmaxPol[0] = lmaxPolNew[0];
	if (polOrd>lmaxPolNew[0] || polOrd>lmaxPol[0])
	  padapt_num[int(lmaxPolNew[0])] += 1;
	LConstDiscFSType lind = ind_.localFunction(*it);
	if (start)
	  lind[0] = lRP[0];
	else
	  lind[0] = 2.*(lRT[0]+lRS[0]+lRP[0])*lrho[0];
	ret += lind[0];
	if (adapt)
	  adapt->addToLocalIndicator(*it,lind[0]);
      }
      return ret;
    }
  };
} // end namespace 

#endif

