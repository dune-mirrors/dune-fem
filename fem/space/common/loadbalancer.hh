#ifndef DUNE_LOADBALANCER_HH
#define DUNE_LOADBALANCER_HH

//- system includes 
#include <vector>

//- local includes 
#include <dune/fem/space/common/datacollector.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune{

/** @defgroup Adaptation LoadBalancer 
 * Concept for calling the grids load balance method 
 * @{
 **/


  
/** @defgroup LoadBalancer 
    @{
**/

class LoadBalancerInterface 
{
public:
  //! default constructor 
  LoadBalancerInterface () : lb_ (0) {}
  
  //! copy constructor 
  LoadBalancerInterface(const LoadBalancerInterface & org) 
    : lb_(org.lb_) 
  {}

  //! destructor 
  virtual ~LoadBalancerInterface () {}
  
  virtual bool loadBalance () const  
  {
    if(lb_) return lb_->loadBalance();  
    else 
    {
      std::cerr << "WARNING: loadBalance with empty lb_ object! \n";
      return false; 
    }
  }

  //! return current balance counter 
  virtual int balanceCounter () const 
  {
    return (lb_) ? lb_->balanceCounter() : 0;
  }

  //! Assignement operator
  LoadBalancerInterface & operator = (const LoadBalancerInterface & lb)
  {
    /** \todo This const-casting seems strange to me! */
    lb_ = &lb;
    return (*this);
  }

private: 
  const LoadBalancerInterface *lb_;
};

/*! \brief This class manages the adaptation process. 
 If the method adapt is called, then the grid is adapted and also 
 all the data belonging to the given dof manager will be rearranged 
 for data set where it is necessary to keep the data.
 */
template <class GridType>
class LoadBalancer 
: public LoadBalancerInterface 
{  
  template <class LBManager, class GridImp, bool isGoodGrid> 
  struct CallLoadBalance
  {
    template <class DofManagerImp>
    static bool balance(const LBManager& lb, 
                        GridImp & grid, 
                        DofManagerImp& dm)
    {
      return grid.loadBalance( dm );
    }
  };
  
  template <class LBManager, class GridImp> 
  struct CallLoadBalance<LBManager,GridImp,false>
  {
    template <class DofManagerImp>
    static bool balance(const LBManager& lb, 
                        GridImp & grid, 
                        DofManagerImp& dm)
    {
      return false;
    }
  };
  
  typedef LoadBalancer<GridType> ThisType;
  typedef DofManager<GridType> DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  // type of data collector during load balance 
  typedef typename DofManagerType :: DataCollectorType DataCollectorType; 

  // type of local data collector interface 
  typedef typename DataCollectorType :: LocalInterfaceType
    LocalDataCollectorInterfaceType;
public:
  /** \brief constructor of LoadBalancer  
     \param grid Grid that load balancing is done for 
     \param paramFile optional parameter file which contains 
        the following two lines:
     # BalanceStep, balancing is done every x-th step, 0 means no balancing    
     BalanceStep: 1 # (do balancing every step)
     \param balanceCounter actual counter, default is zero 
  */   
  LoadBalancer(GridType & grid,
               std::string paramFile = "",
               int balanceCounter = 0) 
    : grid_(grid) 
    , dm_ ( DofManagerFactoryType::getDofManager(grid_) )
    , balanceStep_(0)
    , balanceCounter_(balanceCounter)
    , localList_()
    , collList_()
  {
    if( paramFile != "")
    {
      readParameter(paramFile,"BalanceStep",balanceStep_);
    }

    std::cout << "Created LoadBalancer: balanceStep = " << balanceStep_ << std::endl;
  }

  //! destructor 
  virtual ~LoadBalancer () 
  {
    // remove data collectors 
    for(size_t i=0; i<collList_.size(); ++i)
    {
      delete collList_[i];
    }
     
    // remove local data handler 
    for(size_t i=0; i<localList_.size(); ++i)
    {
      delete localList_[i];
    }
  }

  //! returns actual balanceCounter for checkpointing 
  int balanceCounter () const { return balanceCounter_; }
  
  //! do load balance every balanceStep_ step 
  bool loadBalance () const 
  {
    bool changed = false;
    // if balance counter has readed balanceStep do load balance
    if( (balanceCounter_ >= balanceStep_) && (balanceStep_ > 0) )
    {
      // call load balance 
      changed = CallLoadBalance<ThisType,GridType,
                  Conversion<GridType,HasHierarchicIndexSet>::exists>::
                    balance(*this,grid_,dm_);

      // reset balance counter 
      balanceCounter_ = 0;
    }

    // increase balanceCounter if balancing is enabled 
    if( balanceStep_ > 0 ) ++balanceCounter_;

    return changed;
  }

  //! add discrete function to data inliner/xtractor list 
  template <class DiscreteFunctionType> 
  void addDiscreteFunction(DiscreteFunctionType& df) 
  {
    // to bew revised
    const bool leaf = true; 
    
    ////////////////////////////
    // data inliners 
    ////////////////////////////
    {
      const bool readData = false; // readData is described by false 
      typedef DataInliner<DiscreteFunctionType> LocalInlinerType; 

      LocalInlinerType * di = new LocalInlinerType(df);

      // for later removal 
      localList_.push_back( di );
    
      typedef DataCollector<GridType,LocalInlinerType> DataCollectorImp;

      DataCollectorImp* gdi = 
        new DataCollectorImp( grid_, dm_ , *di , readData , leaf );
      
      // for later removal 
      collList_.push_back(gdi);

      dm_.addDataInliner( *gdi );
    }
   
    ////////////////////////////
    // data xtractors 
    ////////////////////////////
    {
      typedef DataXtractor<DiscreteFunctionType> LocalXtractorType; 
      LocalXtractorType * dx = new LocalXtractorType(df);

      // for later removal 
      localList_.push_back( dx );

      const bool writeData = true; // writedata is described by true 
      
      typedef DataCollector<GridType,LocalXtractorType> DataCollectorImp;
      
      DataCollectorImp* gdx = 
        new DataCollectorImp( grid_, dm_ , *dx , writeData , leaf );
      
      // for later removal 
      collList_.push_back(gdx);

      dm_.addDataXtractor( *gdx );
    }
    
  }
 
protected: 
  //! corresponding grid 
  mutable GridType & grid_;

  //! DofManager corresponding to grid
  mutable DofManagerType & dm_;

  // call loadBalance ervery balanceStep_ step
  int balanceStep_ ;
  // count actual balance call
  mutable int balanceCounter_;

  // list of created local data collectors 
  std::vector<LocalDataCollectorInterfaceType*> localList_;
  std::vector<DataCollectorType*> collList_;
};

/** @} end documentation group */
} // end namespace Dune 
#endif
