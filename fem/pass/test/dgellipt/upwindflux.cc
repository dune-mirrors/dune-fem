#ifndef UPWINDFLUX_CC
#define UPWINDFLUX_CC

class ElliptFluxBase
{
  // see 
  // Unified Analysis for DG Methods for Eliiptic problems 
  // Arnold, Brezzi, Cockburn, Marini (2001)
  
protected:
  // implements {sigma}, which is the just the average 
  // and for simplicity here also multiplied by the normal 
  template <class DomainType, class RangeType, class FluxType > 
  void averageFlux(const DomainType & normal,
        const RangeType & wLeft, const RangeType & wRight,
        FluxType & result) const 
  {
    result  = normal * ( wLeft + wRight );
    result *= 0.5;
  }

  // implements -alpha_j(u_h)
  template <class DomainType, class RangeFieldType> 
  void alphaJot(const DomainType & normal, const RangeFieldType & edgeInv, 
                const RangeFieldType & uLeft, 
                const RangeFieldType & uRight,
                DomainType & result) const
  {
    result  = normal; 
    result *= -(uLeft-uRight)*edgeInv;
  }

  template <class DomainType, class RangeFieldType> 
  void jump(const DomainType & normal, const RangeFieldType & edgeInv, 
                const RangeFieldType & uLeft, 
                const RangeFieldType & uRight,
                DomainType & result) const
  {
    result  = normal; 
    result *= -(uLeft-uRight)*edgeInv;
  }
};


// Numerical Upwind-Flux
template <class ModelImp>
class UpwindFlux {
public:
  typedef ModelImp Model;
  typedef typename Model::Traits Traits;
  typedef typename Traits::FluxRangeType FluxRangeType; 
  //enum { dimRange = Model::dimDomain Range };
public:
  UpwindFlux(const Model& mod) : model_(mod) , velocity_(1.0) {}
  const Model& model() const {return model_;}
  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline double numericalFlux(IntersectionIterator& it,
             double time, 
             const FaceDomainType& x,
             const ArgType& uLeft, 
             const ArgType& uRight,
             RanType& gLeft,
             RanType& gRight) const 
  {
    typedef typename Traits::DomainType DomainType;
    const DomainType normal = it.integrationOuterNormal(x);    

    enum { dim = DomainType :: dimension };

    FluxRangeType ret;

    gLeft  = 0.0;
    gRight = 0.0;
    
    model_.analyticalFlux(*it.inside(),time,x,uLeft,ret);

    ret.umv(normal,gLeft);  
  
    if(it.neighbor())
    {
      FluxRangeType ret1;
      model_.analyticalFlux(*it.outside(),time,x,uRight,ret1);
      ret1.umv(normal,gRight);  
      gLeft += gRight; 
      gLeft *= 0.5;
    }

    double edge = 1.0/normal.two_norm();
    DomainType tmp = normal;
    tmp *= -(uLeft[dim]-uRight[dim])*edge*edge; 
   
    gLeft[dim]  += (normal * tmp);
    
    gRight = gLeft;
    return std::abs(normal.two_norm()); 
  }
private:
  const Model& model_;
  const typename Traits::DomainType velocity_;
};

// Numerical Upwind-Flux
template <class ModelImp>
class IPUpwindFlux {
public:
  typedef ModelImp Model;
  typedef typename Model::Traits Traits;
  typedef typename Traits::FluxRangeType FluxRangeType; 
  enum { dimRange = Model::dimRange };
public:
  IPUpwindFlux(const Model& mod) : model_(mod) , velocity_(1.0) {}
  const Model& model() const {return model_;}
  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline double numericalFlux(IntersectionIterator& it,
             double time, 
             const FaceDomainType& x,
             const ArgType& uLeft, 
             const ArgType& uRight,
             RanType& gLeft,
             RanType& gRight) const 
  {
    typedef typename Traits::DomainType DomainType;
    const DomainType normal = it.integrationOuterNormal(x);    

    enum { dim = DomainType :: dimension };

    FluxRangeType ret;
    gLeft  *= 0.0;
    gRight *= 0.0;
    model_.analyticalFlux(*it.inside(),time,x,uLeft,ret);
    ret.umv(normal,gLeft);  
    
    if(it.neighbor())
    {
      FluxRangeType ret1;
      model_.analyticalFlux(*it.outside(),time,x,uRight,ret1);
      ret1.umv(normal,gRight);  
      gLeft += gRight; 
      gLeft *= 0.5;
    }

    double edge = 1.0/normal.two_norm();

    DomainType w(normal); 
    w *= uLeft[dim];  
    
    DomainType v(normal); 
    // minus because normal from other side 
    v *= -uRight[dim];

    v += w; 
 
    double beta = 100.0;
    gLeft[dim] -= beta*(normal * v);
    double jump = normal * v;  

    for(int i=0; i<dim; ++i) gLeft[i] += jump; 

    DomainType tmp = normal;
    tmp *= -( uLeft[dim] - uRight[dim])*beta*edge; 

    gLeft[dim] += (normal * tmp);
    
    gRight = gLeft;
    return std::abs(normal.two_norm()); 
  }
private:
  const Model& model_;
  const typename Traits::DomainType velocity_;
};

// Numerical Upwind-Flux
template <class ModelImp>
class LDGUpwindFlux : public ElliptFluxBase {
public:
  typedef ModelImp Model;
  typedef typename Model::Traits Traits;
  typedef typename Traits::DomainType DomainImp;
  typedef typename Traits::FluxRangeType FluxRangeType; 
  enum { dimRange = Model::dimRange };
public:
  LDGUpwindFlux(const Model& mod) 
    : model_(mod)
    , upwindDir_(-1.0)
  {
    /*
    upwindDir_[0] = M_PI;
    if(DomainImp::dimension > 1) 
      upwindDir_[1] = M_LN2;
    if(DomainImp::dimension > 2) 
      upwindDir_[2] = std::sqrt(2.0);
      */
  }
  const Model& model() const {return model_;}
  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline double numericalFlux(IntersectionIterator& it,
             double time, 
             const FaceDomainType& x,
             const ArgType& uLeft, 
             const ArgType& uRight,
             RanType& gLeft,
             RanType& gRight) const 
  {
    typedef typename Traits::DomainType DomainType;
    DomainType normal = it.integrationOuterNormal(x);    

    double upwind = normal * upwindDir_;
    //if(upwind > 0.0) normal *= -1.0;

    enum { dim = DomainType :: dimension };
    
    FluxRangeType ret;
    gLeft  = 0.0;
    gRight = 0.0;
    model_.analyticalFlux(*it.inside(),time,x,uLeft,ret);
    ret.umv(normal,gLeft);  
   
    if(it.neighbor())
    {
      model_.analyticalFlux(*it.outside(),time,x,uRight,ret);
      ret.umv(normal,gRight);  
   
      gLeft += gRight; 
      gLeft *= 0.5;
    }

    double edge = normal.two_norm();
    //****************************************************
    // u_h 
    {
      DomainType normalDir(normal); 
      //if(upwind < 0.0) normalDir *= -1.0;
      
      DomainType jump(normal); 
      jump *= uLeft[dim]*edge; 
      DomainType j2 (normal);
      j2 *= -uRight[dim]*edge;
      jump += j2; 

      DomainType beta(normalDir);
      beta *= -0.5*edge;

      double val = jump * beta; 
      for(int i=0; i<dim; ++i) gLeft[i] -=val;
    }
    //****************************************************

    //****************************************************
    // sigma 
    {
      DomainType unitnormal(normal);
      //if(upwind < 0.0) unitnormal *= -1.0;
      
      DomainType w; 
      for(int i=0; i<dim; ++i) w[i] = uLeft[i];
      DomainType v; 
      for(int i=0; i<dim; ++i) v[i] = uRight[i];
      unitnormal *= edge;

      double jump = w * unitnormal; 
      jump -= unitnormal * v;  

      // choice of beta =-0.5*normal
      DomainType beta(normal);
      //if(upwind <0.0) beta *= -1.0;
      beta *= -edge*edge*0.5*jump;

      DomainType tmp;
      
      this->alphaJot(normal,edge,uLeft[dim],uRight[dim],tmp);
      beta += tmp;

      gLeft[dim] += normal * beta;    
    }
    //***************************************************

    gRight = gLeft;
    return std::abs(normal.two_norm()); 
  }
private:
  const Model& model_;
  DomainImp upwindDir_;
};

#if 0
// Numerical Upwind-Flux
template <class ModelImp>
class LDGFlux {
public:
  typedef ModelImp Model;
  typedef typename Model::Traits Traits;
  typedef typename Traits::FluxRangeType FluxRangeType; 
public:
  LDGFlux(const Model& mod) : model_(mod) , velocity_(1.0) {}
  const Model& model() const { return model_; }

  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline void numFluxHalf(IntersectionIterator& it,
             double time, 
             const DomainType & normal,
             const double scaling,  
             const FaceDomainType& x,
             const ArgType& u, 
             const double factor,
             RanType& g) const
  {
    RanType tmp = u * normal;
    g = tmp;
    g *= 0.5;   

    tmp *= factor*scaling; 
    g += tmp; 
  }
  
  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline double numericalFlux(IntersectionIterator& it,
             double time, 
             const FaceDomainType& x,
             const ArgType& uLeft, 
             const ArgWType &wLeft,
             const ArgType& uRight,
             const ArgWType& wRight, 
             RanType& gLeft,
             RanType& gRight) const 
  {
    typedef typename Traits::DomainType DomainType;
    const DomainType normal = it.integrationOuterNormal(x);    

    enum { dim = DomainType :: dimension };

    gLeft  = uLeft*normal;
    gLeft *= 0.5;   

    gRight = uRight*normal;

    double scaling = normal.two_norm();
    scaling *=beta_;
    
    argType tmp;  
    argType tmp2; 
    DomainType v(1.0);
    //tmp=tauleft+tauright;
    tmp *= 0.5;   
                  
    double scaling = normal.two_norm();
    scaling *=beta_;
    tmp2=tauleft-tauright;
    tmp2*=scaling;
    tmp+=tmp2;

    ret=tmp*normal;
    ret +=(phileft-phiright);

    
    /*
    FluxRangeType ret;
    
    model_.analyticalFlux(*it.inside(),time,x,uLeft,ret);

    ret.umv(normal,gLeft);  
  
    if(it.neighbor())
    {
      FluxRangeType ret1;
      model_.analyticalFlux(*it.outside(),time,x,uRight,ret1);
      ret1.umv(normal,gRight);  
      gLeft += gRight; 
      gLeft *= 0.5;
    }

    double edge = 1.0/normal.two_norm();
    DomainType tmp = normal;
    tmp *= -(uLeft[dim]-uRight[dim])*edge*edge; 
   
    gLeft[dim]  += (normal * tmp);
    
    gRight = gLeft;
    */
    return std::abs(normal.two_norm()); 
  }
  
  template <class IntersectionIterator, class FaceDomainType, 
            class ArgType, class RanType> 
  inline double gradientFlux(IntersectionIterator& it,
             double time, 
             const FaceDomainType& x,
             const ArgType& uLeft, 
             const ArgType& uRight,
             RanType& gLeft,
             RanType& gRight) const 
  {
    typedef typename Traits::DomainType DomainType;
    const DomainType normal = it.integrationOuterNormal(x);    

    enum { dim = DomainType :: dimension };
  }

private:
  const Model& model_;
  const typename Traits::DomainType velocity_;
};
#endif


#endif
