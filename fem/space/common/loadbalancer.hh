#ifndef DUNE_LOADBALANCER_HH
#define DUNE_LOADBALANCER_HH

//- system includes 
#include <cassert>
#include <vector>
#include <set>
#include <iostream>

//- local includes 
#include <dune/common/timer.hh>
#include <dune/fem/space/common/datacollector.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/io/file/asciiparser.hh>


namespace Dune{

/** @addtogroup LoadBalancer 
    In this module a concept for calling the grids load balance method
    is described and implemented.

    \remarks The interface for a LoadBalancer object is described by the
    class LoadBalancerInterface.
    @{
 **/

/** \brief Interface class for load balancing. 
*/
class LoadBalancerInterface 
{
protected:  
  //! default constructor 
  LoadBalancerInterface () {}
  
public:
  //! destructor 
  virtual ~LoadBalancerInterface () {}
  
  /** \brief call load balance, returns true if grid was changed 
    \return \b true if grid was changed, \b false otherwise 
  */
  virtual bool loadBalance () = 0; 

  /** \brief return number of cycles since last application of load balance 
    \return number of cycles since last application of load balance 
  */
  virtual int balanceCounter () const = 0;

  /** \brief time that last load balance cycle took */
  virtual double loadBalanceTime () const 
  {
    return 0.0;
  }
};

/*! \brief This class manages the adaptation process. 
 If the method adapt is called, then the grid is adapted and also 
 all the data belonging to the given dof manager will be rearranged 
 for data set where it is necessary to keep the data.
 */
template <class GridType>
class LoadBalancer 
: virtual public LoadBalancerInterface 
{  
  // type of this 
  typedef LoadBalancer<GridType> ThisType;
  // dof manager 
  typedef DofManager<GridType> DofManagerType; 
  // factory 
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
    const bool output = (grid_.comm().rank() == 0);
    if( paramFile != "")
    {
      readParameter(paramFile,"BalanceStep",balanceStep_, output);
    }

    if( output )
    {
      std::cout << "Created LoadBalancer: balanceStep = " << balanceStep_ << std::endl;
    }
  }

  /** \brief constructor of LoadBalancer  
     \param grid Grid that load balancing is done for 
     \param rpOp restrict prolong tpye 
     \param paramFile optional parameter file which contains 
        the following two lines:
     # BalanceStep, balancing is done every x-th step, 0 means no balancing    
     BalanceStep: 1 # (do balancing every step)
     \param balanceCounter actual counter, default is zero 
  */   
  template <class RestrictProlongTpye> 
  LoadBalancer(GridType & grid,
               RestrictProlongTpye& rpOp,
               std::string paramFile = "",
               int balanceCounter = 0) 
    : grid_(grid) 
    , dm_ ( DofManagerFactoryType::getDofManager(grid_) )
    , balanceStep_(0)
    , balanceCounter_(balanceCounter)
    , localList_()
    , collList_()
    , balanceTime_(0.0)
  {
    const bool output = (grid_.comm().rank() == 0);
    if( paramFile != "")
    {
      readParameter(paramFile,"BalanceStep",balanceStep_, output);
    }

    rpOp.addToList(*this);
    if( output ) 
    {
      std::cout << "Created LoadBalancer: balanceStep = " << balanceStep_ << std::endl;
    }
  }

  //! destructor 
  virtual ~LoadBalancer () 
  {
    // clear objects from dof managers list 
    dm_.clearDataInliners();
    dm_.clearDataXtractors();

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
  bool loadBalance () 
  {
    // get stopwatch 
    Timer timer ; 
    
    bool changed = false;
    // if balance counter has readed balanceStep do load balance
    if( (balanceCounter_ >= balanceStep_) && (balanceStep_ > 0) )
    {
      // call grids load balance, only implemented in ALUGrid right now
      changed = grid_.loadBalance( dm_ ); 

      // reset balance counter 
      balanceCounter_ = 0;
    }

    // increase balanceCounter if balancing is enabled 
    if( balanceStep_ > 0 ) ++balanceCounter_;

    // get time 
    balanceTime_ = timer.elapsed();

    return changed;
  }

  /** @copydoc LoadBalancerInterface::loadBalanceTime */
  virtual double loadBalanceTime() const 
  {
    return balanceTime_;
  }

  //! add discrete function to data inliner/xtractor list 
  template <class DiscreteFunctionType>
  void addToList(DiscreteFunctionType& df)
  {
    addDiscreteFunction(df);    
  }

  //! add discrete function to data inliner/xtractor list 
  template <class DiscreteFunctionType> 
  void addDiscreteFunction(DiscreteFunctionType& df) 
  {
    CompileTimeChecker< Conversion<DiscreteFunctionType,IsDiscreteFunction>::exists > 
      only_valid_for_discretefunctions();

    const IsDiscreteFunction * fct = &df;

    // if discrete functions is not in list already 
    if( listOfFcts_.find(fct) == listOfFcts_.end() )
    {
      // insert into set 
      listOfFcts_.insert( fct );

      // to be revised
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

      // enable this discrete function for dof compression 
      df.enableDofCompression();
    }
  }
 
protected: 
  //! corresponding grid 
  GridType & grid_;

  //! DofManager corresponding to grid
  DofManagerType & dm_;

  // call loadBalance ervery balanceStep_ step
  int balanceStep_ ;
  // count actual balance call
  int balanceCounter_;

  // list of created local data collectors 
  std::vector<LocalDataCollectorInterfaceType*> localList_;
  std::vector<DataCollectorType*> collList_;

  // list of already added discrete functions 
  std::set< const IsDiscreteFunction * > listOfFcts_;

  // time for last load balance call
  double balanceTime_;
};

/** @} end documentation group */
} // end namespace Dune 
#endif
