// Gitterdefinition?
#ifndef DUNE_DGOPERATORS_HH
#define DUNE_DGOPERATORS_HH

//- system includes 
#include <string>

// Dune includes
// #include <dune/common/utility.hh>
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include "discretemodels.hh"

//#include "limitpass.hh"
#include <dune/fem/pass/limitpass.hh>

//*************************************************************
namespace Dune {  
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGAdvectionDiffusionOperator : 
  public SpaceOperatorInterface<
	   //typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType
     typename PassTraits< Model , Model::Traits::dimRange , polOrd >::DestinationType
  >
#if 0
  public Operator<double,double,
	   typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType,
	   typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType> 
#endif
  {
  public:

    enum PassIdType{ u , pass1 , pass2 };  
    
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef AdvDiffDModel1< Model , NumFluxType , polOrd
              , (int)u > DiscreteModel1Type;
    
    typedef AdvDiffDModel2< Model , NumFluxType , polOrd
              , (int)u , (int)pass1 > DiscreteModel2Type;
    
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass< DiscreteFunction2Type , (int)u > Pass0Type;
    typedef LocalDGPass< DiscreteModel1Type , Pass0Type , (int)pass1 > Pass1Type;
    typedef LocalDGPass< DiscreteModel2Type , Pass1Type , (int)pass2 > Pass2Type;

    typedef typename Traits1::DiscreteFunctionSpaceType 
      Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType  
      Space2Type;
    typedef typename Traits1::DestinationType Destination1Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits1::GridPartType GridPartType;

  public:
    DGAdvectionDiffusionOperator(GridType& grid,
				 const NumFluxType& numf,
				 DomainType& upwind) :
      upwind_(upwind),
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      gridPart_(grid_),
      space1_(gridPart_),
      space2_(gridPart_),
      problem1_(upwind_,model_,numflux_),
      problem2_(upwind_,model_,numflux_),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_) 
    {}
    void setTime(double time) {
      pass2_.setTime(time);
    }
    double timeStepEstimate() const {
      return pass2_.timeStepEstimate();
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
      // upwind_*=(-1.);
    }
    void limit(const DestinationType& arg,DestinationType& dest) const {
    }
    void switchupwind() {upwind_*=(-1.);}
    const SpaceType& space() const {
      return space2_;
    }
    
    void printmyInfo(std::string filename) const {
      std::ostringstream filestream;
      filestream << filename;
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Advection-Diffusion Op., polynomial order: " << polOrd << "\\\\\n\n";
      ofs.close();
    }
  private:
    mutable DomainType upwind_;
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType gridPart_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
/**************************************************************/
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGAdvectionOperator : 
  public SpaceOperatorInterface<
	   //typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType
     typename PassTraits< Model , Model::Traits::dimRange , polOrd >::DestinationType
  >
#if 0
    public Operator<double,double,
		    typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,false,true>::Traits::DiscreteFunctionType,
		    typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,false,true>::Traits::DiscreteFunctionType> 
#endif
  {
  public:

    enum PassIdType{ u = 25 , pass1 = 5 , pass2 = 24 };  

    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef AdvDiffDModel3< Model , NumFluxType , polOrd , u > DiscreteModel2Type;
     typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass< DiscreteFunction2Type , u > Pass0Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass0Type, pass2 > Pass2Type;

    // typedef typename Traits2::SingleDiscreteFunctionSpaceType SDFSType;
    typedef typename Traits2::DiscreteFunctionSpaceType Space2Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits2::GridPartType GridPartType;
    typedef Space2Type DiscreteFunctionSpaceType;

  public:
    DGAdvectionOperator(GridType& grid,
			const NumFluxType& numf,
			DomainType& upwind) :
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      gridPart_(grid_),
      // sspace2_(gridPart_),
      // space2_(sspace2_),
      space2_(gridPart_),
      problem2_(model_,numflux_),
      pass2_(problem2_, pass0_, space2_) 
    {}
    void setTime(double time) {
      pass2_.setTime(time);
    }
    double timeStepEstimate() const {
      return pass2_.timeStepEstimate();
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
    }
    const SpaceType& space() const {
      return space2_;
    }
    void switchupwind() {}
    void limit(const DestinationType& arg,DestinationType& dest) const {
    }
    void printmyInfo(std::string filename) const {
      std::ostringstream filestream;
      filestream << filename;
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Advection Op., polynomial order: " << polOrd << "\\\\\n\n";
      ofs.close();
    }
  private:
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType gridPart_;
    // SDFSType sspace2_;
    Space2Type space2_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass2Type pass2_;
  };
/**************************************************************/
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGDiffusionOperator : 
  public SpaceOperatorInterface<
	   //typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType
     typename PassTraits< Model , Model::Traits::dimRange , polOrd >::DestinationType
  >
#if 0
    public Operator<double,double,
		    typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,false>::Traits::DiscreteFunctionType,
		    typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,false>::Traits::DiscreteFunctionType> 
#endif 
  {
  public:

    enum PassIdType{ u = 25 , pass1 = 5 , pass2 = 24 };  
    
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef AdvDiffDModel1<Model,NumFluxType,polOrd , u > DiscreteModel1Type;
    
    typedef AdvDiffDModel4< Model , NumFluxType , polOrd
              , u , pass1 > DiscreteModel2Type;
    
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass< DiscreteFunction2Type , u > Pass0Type;
    typedef LocalDGPass<DiscreteModel1Type, Pass0Type , pass1 > Pass1Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass1Type , pass2 > Pass2Type;

    typedef typename Traits1::DiscreteFunctionSpaceType 
      Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType  
      Space2Type;
    typedef typename Traits1::DestinationType Destination1Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits1::GridPartType GridPartType;

  public:
    DGDiffusionOperator(GridType& grid,
			const NumFluxType& numf,
			DomainType& upwind) :
      upwind_(upwind),
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      gridPart_(grid_),
      space1_(gridPart_),
      space2_(gridPart_),
      problem1_(upwind_,model_,numflux_),
      problem2_(upwind_,model_,numflux_),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_) 
    {}
    void setTime(double time) {
      pass2_.setTime(time);
    }
    double timeStepEstimate() const {
      return pass2_.timeStepEstimate();
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
      // upwind_*=(-1.);
    }
    void limit(const DestinationType& arg,DestinationType& dest) const {
		}
    void switchupwind() {upwind_*=(-1.);}
    const SpaceType& space() const {
      return space2_;
    }
    void printmyInfo(std::string filename) const {
      std::ostringstream filestream;
      filestream << filename;
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Diffusion Op., polynomial order: " << polOrd << "\\\\\n\n";
      ofs.close();
    }
  private:
    mutable DomainType upwind_;
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType gridPart_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
  // *********************************************************
  // *********************************************************
  // *********************************************************
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGLimitedAdvectionOperator : 
    public Operator<double,double,
	   //typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType,
	   //typename AdvDiffDModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType
     typename PassTraits< Model , Model::Traits::dimRange , polOrd >::DestinationType ,
     typename PassTraits< Model , Model::Traits::dimRange , polOrd >::DestinationType
       > 
  {
  public:
    
    enum PassIdType{ u = 25 , pass1 = 5 , pass2 = 24 };  
    
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    // define limiter discrete model 
    typedef LimiterDefaultDiscreteModel<PassTraits<Model,dimRange,polOrd>, Model> DiscreteModel1Type;

    typedef AdvDiffDModel4< Model , NumFluxType , polOrd
              , u > DiscreteModel2Type;
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass< DiscreteFunction2Type , u > Pass0Type;
    typedef LimitDGPass<DiscreteModel1Type, Pass0Type , pass1 > Pass1Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass1Type , pass2 > Pass2Type;

    typedef typename Traits1::DiscreteFunctionSpaceType 
      Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType  
      Space2Type;
    typedef typename Traits1::DestinationType Destination1Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits1::GridPartType GridPartType;

  public:
    DGLimitedAdvectionOperator(GridType& grid,
			       const NumFluxType& numf,
			       DomainType& upwind) :
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      gridPart_(grid_),
      space1_(gridPart_),
      space2_(gridPart_),
      problem1_(model_),
      problem2_(model_,numf),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_) 
    {}
    void setTime(double time) {
      pass2_.setTime(time);
    }
    double timeStepEstimate() const {
      return pass2_.timeStepEstimate();
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
    }
    void limit(const DestinationType& arg,DestinationType& dest) const {
      pass1_(arg,dest);
    }
    void switchupwind() {}
    const SpaceType& space() const {
      return space2_;
    }
    void printmyInfo(std::string filename) const {
      std::ostringstream filestream;
      filestream << filename;
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Limit-Advection Op., polynomial order: " << polOrd << "\\\\\n\n";
      ofs.close();
    }
  private:
    mutable DomainType upwind_;
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType gridPart_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
}
#endif
