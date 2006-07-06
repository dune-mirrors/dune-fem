int numPol[3]={1,3,6};
template <class LocalFuncType>
struct LocalFuncHelper {
  //enum { dimRange = LocalFuncType::dimRange };
  //enum { dimDomain = LocalFuncType::dimDomain };
  enum { dimRange = 1 };
  enum { dimDomain = 2 };
  typedef typename LocalFuncType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename LocalFuncType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename LocalFuncType::RangeFieldType RangeFieldType;
  typedef typename LocalFuncType::DomainType DomainType;
  typedef typename LocalFuncType::RangeType RangeType;
  typedef typename LocalFuncType::JacobianRangeType JacobianRangeType;
  typedef FieldMatrix<double, dimDomain, dimDomain> JacobianInverseType;
  const DiscreteFunctionSpaceType& spc_;
  const BaseFunctionSetType* bSet_;
  JacobianRangeType tmpJRT_,tmpGrad_;
  RangeType tmpRT_;
  int numDofs_;
  vector<double> dat_;
  /****************/
  LocalFuncHelper(const DiscreteFunctionSpaceType& spc) :
    spc_(spc),
    bSet_(0), tmpJRT_(0), tmpGrad_(0), tmpRT_(0), numDofs_(-1), dat_(0) {
  }
  template <class EntityType>
  void init(EntityType& en) {
    bSet_ = &(spc_.getBaseFunctionSet(en));
    numDofs_ = bSet_->numBaseFunctions();
    dat_.resize(numDofs_);
  }
  void assign(const LocalFuncType& lf) {
    for (int i=0;i<numDofs_;++i) 
      dat_[i] = lf[i];
  }
  void assign(const double a) {
    for (int i=0;i<numDofs_;++i) 
      dat_[i] = a;
  }
  void addscaled(LocalFuncType& lf,double l) {
    for (int i=0;i<numDofs_;++i) 
      dat_[i] += l*lf[i];
  }
  /****************/
  RangeFieldType& operator [] (int num) {
    return dat_[num];
  }
  const RangeFieldType& operator [] (int num) const {
    return dat_[num];
  }
  int numDofs () const {
    return numDofs_;
  }
  template <class EntityType,class QuadratureType>
  void ucomponents(const EntityType& en,
		   const QuadratureType& quad,int l,
		   std::vector<RangeType>& comp) const {
    for (int i = 0; i < numDofs_; ++i) {
      {
        bSet_->eval(i, quad,l, tmpRT_);
        for (int l = 0; l < dimRange; ++l) {
  	      comp[i][l] = dat_[i] * tmpRT_[l];
        }
      }
    }
  }
  template <class EntityType,class QuadratureType>
  void evaluateLocal(const EntityType& en,
		     const QuadratureType& quad,int p,
		     int maxp,
		     RangeType & ret) {
    ret    = 0.0;
    tmpRT_ = 0.0;
    for (int i = 0; i < numDofs_; ++i) {
      if (i<numPol[maxp]) {
        bSet_->eval(i, quad,p, tmpRT_);
       for (int l = 0; l < dimRange; ++l) {
	  ret[l] += dat_[i] * tmpRT_[l];
        }
      }
    }
  }
  template <class EntityType,class QuadratureType>
  void jacobianLocal(const EntityType& en,
		     const QuadratureType& quad,int p,
         int maxp,
         JacobianRangeType& ret) {
    tmpGrad_ = 0;
    tmpJRT_ = 0;
    ret = 0.0;
    const JacobianInverseType& jti =
      en.geometry().jacobianInverseTransposed(quad.point(p));
    for (int i = 0; i < numDofs_; ++i) {
      // tmpGrad_ = 0;
      if (i<numPol[maxp]) { 
        bSet_->jacobian(i, quad,p, tmpGrad_);
        tmpGrad_ *= dat_[i];
        tmpJRT_  += tmpGrad_;
      }
    }
    for (int l = 0; l < dimRange; ++l)
      jti.umv(tmpJRT_[l],ret[l]);
  }
};
