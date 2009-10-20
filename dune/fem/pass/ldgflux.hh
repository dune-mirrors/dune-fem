#ifndef DUNE_LDGFLUX_CC
#define DUNE_LDGFLUX_CC

#include <cassert>
#include <cmath>

// Numerical Upwind-Flux
template <class ModelImp>
class LDGFlux
{
public:
  typedef ModelImp ModelType;
  typedef typename ModelType::Traits::DomainType DomainType;
public:
  LDGFlux(const ModelType& mod, 
          const double beta, 
          const double power,
          const double eta) 
    : model_(mod)
    , beta_(beta)
    , power_(power)
    , eta_(eta)
    , betaNotZero_( std::abs( beta_ ) > 0.0 )
  {
  }

  const ModelType& model() const {return model_;}

  //! evaluates { sigma } = 0.5 ( sigmaLeft + sigmaRight )
  template <class URangeType>
  inline double sigmaFluxBetaZero(
                       const DomainType & unitNormal,
                       const double faceVol, 
                       const URangeType & uLeft,
                       const URangeType & uRight, 
                       URangeType & sigmaLeft,
                       URangeType & sigmaRight) const
  {
    sigmaLeft  = uLeft;
    sigmaRight = uRight;
    
    sigmaLeft  *= 0.5;
    sigmaRight *= 0.5;

    return 0.0;
  }
  
  //! evaluates sigmaBetaZero + stabilization 
  template <class URangeType>
  inline double sigmaFluxStability(const DomainType & unitNormal,
                            const double faceVol, 
                            const URangeType & uLeft,
                            const URangeType & uRight, 
                            URangeType & gLeft,
                            URangeType & gRight) const
  {
    // stabilization term 
    const double factor = eta_ / faceVol;
    gLeft   =  uLeft;
    gRight  = -uRight;

    gLeft  *= factor;
    gRight *= factor;

    return 0.0;
  }

  //! evaluates sigmaBetaZero + stabilization 
  template <class URangeType>
  inline double sigmaFlux(const DomainType & unitNormal,
                   const double faceVol, 
                   const URangeType & uLeft,
                   const URangeType & uRight, 
                   URangeType & sigmaLeft,
                   URangeType & sigmaRight,
                   URangeType & gLeft,
                   URangeType & gRight) const
  {
    // call part of flux for beta == 0
    sigmaFluxBetaZero(unitNormal,faceVol,uLeft,uRight,sigmaLeft,sigmaRight);

    // only calculate this part if |beta| > 0
    if( betaNotZero_ )
    {
      const double scaling = beta_ * std::pow(faceVol,power_);

      URangeType jumpLeft(uLeft);
      URangeType jumpRight(uRight);

      jumpLeft  *= scaling;
      jumpRight *= scaling;

      sigmaLeft  -= jumpLeft;
      sigmaRight += jumpRight;
    }

    sigmaFluxStability(unitNormal,faceVol,uLeft,uRight,gLeft,gRight);
    return 0.0;
  }

  //! evaluates { u } = 0.5 * ( uLeft + uRight )
  template <class SigmaRangeType>
  inline double uFluxBetaZero(const DomainType & unitNormal,
               const double faceVol, 
               const SigmaRangeType & uLeft,
               const SigmaRangeType & uRight, 
               SigmaRangeType & sigmaLeft,
               SigmaRangeType & sigmaRight,
               SigmaRangeType & gLeft,
               SigmaRangeType & gRight) const
  {
    // this flux has not sigma parts 
    sigmaLeft  = 0.0;
    sigmaRight = 0.0;

    gLeft  = uLeft; 
    gRight = uRight; 
    
    gLeft  *= 0.5;
    gRight *= 0.5;
    
    return 0.0;
  }


  //! evaluates uFLuxBetaZero + stabilization 
  template <class SigmaRangeType>
  inline double uFlux(const DomainType & unitNormal,
               const double faceVol, 
               const SigmaRangeType & uLeft,
               const SigmaRangeType & uRight, 
               SigmaRangeType & sigmaLeft,
               SigmaRangeType & sigmaRight,
               SigmaRangeType & gLeft,
               SigmaRangeType & gRight) const
  {
    // call part of flux for beta == 0
    uFluxBetaZero(unitNormal,faceVol,uLeft,uRight,
                  sigmaLeft,sigmaRight,gLeft,gRight);

    // only calculate this part if |beta| > 0
    if( betaNotZero_ )
    {
      SigmaRangeType left(uLeft);
      SigmaRangeType right(uRight); 
      
      const double scaling = beta_ * std::pow(faceVol,power_);
      //const double scaling = beta_ * faceVol;

      left  *= scaling;
      right *= scaling;

      gLeft  -= left;
      gRight += right;
    }
    return 0.0;
  }

private:
  const ModelType& model_;
  const double beta_;
  const double power_;
  const double eta_;
  const bool betaNotZero_;
};


//! Flux for Gradient calculation 
class GradientFlux
{
public:
  //! constructor taking beta and power 
  GradientFlux(const double beta, 
               const double power)
    : beta_(beta)
    , power_(power)
  {
    assert( std::abs( beta_ ) > 0.0 );
  }
  
  //! copy constructor 
  GradientFlux(const GradientFlux& org) 
    : beta_(org.beta_), power_(org.power_) {} 

  //! evaluates: result = { u } + n * [ u ] 
  //! ( see Brezzi et al )
  template <class URangeType>
  inline double uFlux(const double faceVol, 
               const URangeType & uLeft,
               const URangeType & uRight, 
               URangeType & result) const 
  {
    const double scaling = beta_ * std::pow(faceVol,power_);
   
    // the following is done:
    // flux = uLeft + 0.5 * (uLeft - uRight) * scaling
    
    result  = uLeft; 
    result -= uRight;
  
    result *= 0.5 * scaling;
    
    result += uLeft;
    
    return 0.0;
  }
  
private:
  const double beta_;
  const double power_;
};

//! Flux for Gradient calculation 
class AverageFlux
{
public:
  //! constructor taking beta and power 
  AverageFlux(const double beta, 
               const double power)
  {
  }
  
  //! evaluates: result = { u } + n * [ u ] 
  //! ( see Brezzi et al )
  template <class URangeType>
  inline double uFlux(const double faceVol, 
               const URangeType & uLeft,
               const URangeType & uRight, 
               URangeType & result) const 
  {
    // the following is done:
    // flux = uLeft + 0.5 * (uLeft - uRight)
    
    result  = uLeft; 
    result -= uRight;
  
    result *= 0.5;
    
    result += uLeft;
    
    return 0.0;
  }
};

#endif
