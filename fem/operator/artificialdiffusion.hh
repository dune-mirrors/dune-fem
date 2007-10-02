#ifndef DUNE_ARTIFICIALDIFFUSION_HH 
#define DUNE_ARTIFICIALDIFFUSION_HH 

//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

//- Dune-fem includes 
#include <dune/fem/function/common/temporarylocalfunction.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/space/fvspace.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/misc/gridwidth.hh>

namespace Dune {

/////////////////////////////////////////////////////////////////////////////
//
//  See: P.O. Persson and J. Peraire. Sub-Cell Shock Capturing 
//       for Discontinuous Galerkin Methods.
//
/////////////////////////////////////////////////////////////////////////////
template <class DiscreteFunctionImp>
class ArtificialDiffusion
{
  typedef DiscreteFunctionImp DiscreteFunctionType;
public:
  typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;

  typedef typename GridType :: template Codim<0> :: Entity EntityType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  typedef FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0> DiffusionSpaceType; 
  typedef AdaptiveDiscreteFunction< DiffusionSpaceType > DiffusionFunctionType;

  const DiscreteFunctionType& df_;
  DiffusionSpaceType diffusionSpace_;
  DiffusionFunctionType diffusion_;
public:
  //! constructor 
  ArtificialDiffusion(const DiscreteFunctionType& df) 
    : df_(df), diffusionSpace_( const_cast<GridPartType&> (df.space().gridPart()))
    , diffusion_(df_.name() + "-art-diff", diffusionSpace_ )
  {}

  //! update artificial diffusion 
  void update()
  {
    calculate(df_, diffusion_ );
  }

  //! return const reference to artificial diffusion 
  const DiffusionFunctionType& diffusion() const { return diffusion_; }

  //! return const reference to artificial diffusion 
  void evaluateLocal (const EntityType& en, RangeType& diff) const 
  { 
    typedef typename DiffusionFunctionType :: LocalFunctionType LocalFunctionType;   
    LocalFunctionType lf = diffusion_.localFunction(en); 
    assert( lf.numDofs() == RangeType :: dimension );
    for(int i=0; i<RangeType :: dimension; ++i) 
    {
      diff[i] = lf[i];
    }
  }

  static double artificialDiffusion(const double s_0,
                                    const double kappa,
                                    const double epsilon_0,
                                    const double s_e)
  {
    const double sKappaLower = s_0 - kappa;
    const double sKappaUpper = s_0 + kappa;
    
    if( s_e < sKappaLower) return 0.0;
    if( (sKappaLower <= s_e) && (s_e <= sKappaUpper) ) 
    {
      return 0.5 * epsilon_0 * (1.0 + sin(M_PI * (s_e - s_0) * 0.5 / kappa));
    }
    if(s_e > sKappaUpper) return epsilon_0;

    DUNE_THROW(InvalidStateException,"Aborted in ArtificialDiffusion kappa");
  }
   
  static double  calculate(DiscreteFunctionType &discFunc, 
                           DiffusionFunctionType &diffusion) 
  {
    return Calc<0,DiscreteFunctionSpaceType::polynomialOrder>::doCalc(discFunc,diffusion);
  }
  
  // calculation of artificial diffusion 
  template <int dummy, int polOrd> 
  struct Calc 
  {
    template <class DiscreteFunctionType, class DiffusionFunctionType> 
    static double doCalc(DiscreteFunctionType &discFunc, 
                         DiffusionFunctionType &diffusion) 
    {
      typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;
      typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
      typedef typename GridType :: template Codim<0> :: Entity EntityType;
      typedef typename GridType :: template Codim<0> :: Geometry Geometry;
      typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
      typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 

      enum { polOrd = DiscreteFunctionSpaceType :: polynomialOrder };

      typedef  DiscontinuousGalerkinSpace<FunctionSpaceType,GridPartType,polOrd-1,
        CachingStorage> LowerSpaceType; 
     
      enum { dim = GridType::dimension };
      enum { dimRange = DiscreteFunctionSpaceType:: DimRange };
      
      // get space 
      const DiscreteFunctionSpaceType& space =  discFunc.space();
      // get grid part  
      GridPartType& gridPart = const_cast<GridPartType&> (space.gridPart());  
      
      // clear destination 
      diffusion.clear();

      // create space with polOrd - 1 
      LowerSpaceType lowerSpace(gridPart);

      // type of temporary local function belonging to lower space 
      typedef TemporaryLocalFunction< LowerSpaceType > TemporaryLocalFunctionType; 
      TemporaryLocalFunctionType uHat( lowerSpace );

      typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      typedef typename DiffusionFunctionType::LocalFunctionType LocalDiffType;
      
      // Get quadraturae rule
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
       
      RangeType ret (0.0);
      RangeType tmp (0.0);

      // maximal diffusion factor 
      double maxDiff = 0.0;

      // constant s_0 
      const double s_0 = 1./(polOrd*polOrd*polOrd*polOrd);

      // constant kappa, sufficiently large 
      const double kappa = 1e2 * s_0;
        
      // grid width / polOrd 
      const double h = GridWidth ::calcGridWidth( gridPart ); 
      const double epsilon_0 = h/polOrd;  

      const int quadOrd = (2 * space.order());
      
      Iterator endit = space.end();
      for(Iterator it = space.begin(); it != endit ; ++it) 
      {
        // get entity 
        const EntityType& en = *it;
        // get geometry 
        const Geometry& geo = en.geometry();

        // get local functions 
        LocalFuncType lf = discFunc.localFunction(en);

        // initialize uHat and set to zero 
        uHat.init( en );
        const int numDofs = uHat.numDofs();
        for(int i=0; i<numDofs; ++i)
        {
          uHat[i] = 0;
        }
       
        // get quadrature 
        typedef CachingQuadrature <GridPartType , 0> QuadratureType; 
        QuadratureType quad(en,quadOrd);
        const int quadNop = quad.nop();

        // L2 Projection 
        {
          typedef typename LowerSpaceType :: BaseFunctionSetType  BaseFunctionSetType; 
          //! Note: BaseFunctions must be ortho-normal!!!!
          const BaseFunctionSetType& baseset = uHat.baseFunctionSet();
          for(int qP = 0; qP < quadNop ; ++qP)
          {
            lf.evaluate(quad,qP, ret);
            for(int i=0; i<numDofs; ++i)
            {
              baseset.evaluate(i,quad,qP, tmp);
              uHat[i] += quad.weight(qP) * (ret * tmp) ;
            }
          }
        }
        
        RangeType uPhi_1(0.0);
        RangeType uPhi(0.0);

        // now evaluate L2 scalar product 
        for(int qp=0; qp<quadNop; ++qp)
        {
          // calculate integration factor 
          const double intel = geo.integrationElement(quad.point(qp)) * quad.weight(qp);

          // eval function 
          lf.evaluate(quad,qp,ret);
          
          for(int i=0; i<dimRange; ++i)
          {
            uPhi[i] += intel * (ret[i] * ret[i]);
          }
          
          // eval p-1 function 
          uHat.evaluate(quad,qp,tmp);
          ret -= tmp;

          for(int i=0; i<dimRange; ++i)
          {
            uPhi_1[i] += intel * (ret[i] * ret[i]);
          }
        }

        LocalDiffType diff = diffusion.localFunction(en);
        for(int i=0; i<dimRange; ++i)
        {
          const RangeFieldType uphi = uPhi[i];
          // if uphi is zero the hole function is zero and 
          // we don't need to do any artificial diffusion 
          if( std::abs(uphi) > 0.0 )
          {
            // if S_e == 0 , then function already constant 
            const double S_e = uPhi_1[i] / uphi;
            if( S_e > 1e-15 )
            {
              const double s_e = log10(S_e);
              const double artDiff = artificialDiffusion(s_0,kappa,epsilon_0,s_e);
              diff[i] = artDiff; 
              maxDiff = std::max( artDiff, maxDiff );
            }
          }
        }
      } // end iteration over elements 
      
      //diffusion.print( std::cout );
      return maxDiff;
    }
  };
  
  template <int dummy> 
  struct Calc<dummy,0> 
  {
    template <class DiscreteFunctionType, class DiffusionFunctionType> 
    static double doCalc(DiscreteFunctionType &discFunc, 
                         DiffusionFunctionType &diffusion) 
    {
      diffusion.clear();
      std::cerr << "PolOrd = 0, doing nothing! \n";
      return 0.0;
    }
  };
};  

} // end namespace 
#endif
