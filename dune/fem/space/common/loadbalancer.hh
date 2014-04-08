#ifndef DUNE_FEM_LOADBALANCER_HH
#define DUNE_FEM_LOADBALANCER_HH

#include <cassert>
#include <iostream>
#include <set>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/tuples.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem/space/common/datacollector.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem 
  {

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
    : virtual public LoadBalancerInterface ,
      public AutoPersistentObject
    {  
      // type of this 
      typedef LoadBalancer<GridType> ThisType;
      // dof manager 
      typedef DofManager<GridType> DofManagerType; 

      // type of data inlining during load balance 
      typedef typename DofManagerType :: DataInlinerType DataInlinerType; 

      // type of data extraction during load balance 
      typedef typename DofManagerType :: DataXtractorType DataXtractorType; 

      // type of local data collector interface 
      typedef typename DataInlinerType :: LocalInterfaceType   LocalDataInlinerInterfaceType;
      // type of local data collector interface 
      typedef typename DataXtractorType :: LocalInterfaceType  LocalDataXtractorInterfaceType;

      typedef std::pair< LocalDataInlinerInterfaceType*, LocalDataXtractorInterfaceType*  >  LocalDataCollectorPairType;
      typedef std::pair< DataInlinerType* , DataXtractorType* > DataCollectorPairType;
    protected:
      /** \brief constructor of LoadBalancer  
         The following optional parameter is used from the Parameter class:
         # BalanceStep, balancing is done every x-th step, 0 means no balancing    
         BalanceStep: 1 # (do balancing every step)
         \param grid Grid that load balancing is done for 
         \param rpOp restrict prolong tpye 
         \param balanceCounter actual counter, default is zero 
      */   
      template< class RestrictProlongOperator > 
      LoadBalancer ( GridType & grid,
                     RestrictProlongOperator &rpOp,
                     int balanceCounter = 0 )
      : grid_( grid ),
        dm_ ( DofManagerType::instance( grid_ ) ),
        balanceStep_( getBalanceStep( balanceCounter ) ),
        balanceCounter_( balanceCounter ),
        localList_(),
        collList_(),
        balanceTime_( 0.0 )
      {
        rpOp.addToLoadBalancer( *this );
        if( Parameter::verbose() )
          std::cout << "Created LoadBalancer: balanceStep = " << balanceStep_ << std::endl;
      }

      /** \brief constructor of LoadBalancer  
         The following optional parameter is used from the Parameter class:
         # BalanceStep, balancing is done every x-th step, 0 means no balancing    
         BalanceStep: 1 # (do balancing every step)
         \param grid Grid that load balancing is done for 
         BalanceStep: 1 # (do balancing every step)
         \param balanceCounter actual counter, default is zero 
      */   
      explicit LoadBalancer ( GridType &grid,
                              int balanceCounter = 0 )
      : grid_( grid ),
        dm_ ( DofManagerType::instance( grid_ ) ),
        balanceStep_( getBalanceStep( balanceCounter ) ),
        balanceCounter_( balanceCounter ),
        localList_(),
        collList_(),
        balanceTime_( 0.0 )
      {
        if( Parameter::verbose() )
          std::cout << "Created LoadBalancer: balanceStep = " << balanceStep_ << std::endl;
      }

      int getBalanceStep( int balanceCounter ) const 
      {
        int step = balanceCounter;
        if( Parameter :: exists("BalanceStep") )
        {
          std::cout << "Warning: deprecated parameter 'BalanceStep', please use 'fem.loadbalancing.step' instead!" << std::endl;
          step = Parameter::getValue< int >( "BalanceStep", balanceCounter );
        }
        else 
          step = Parameter::getValue< int >( "fem.loadbalancing.step", balanceCounter );

        return step;
      }

    public:  
      //! destructor 
      virtual ~LoadBalancer () 
      {
        // clear objects from dof managers list 
        dm_.clearDataInliners();
        dm_.clearDataXtractors();

        // remove data collectors 
        for(size_t i=0; i<collList_.size(); ++i)
        {
          delete collList_[ i ].first  ;
          delete collList_[ i ].second ;
        }
         
        // remove local data handler 
        for(size_t i=0; i<localList_.size(); ++i)
        {
          delete localList_[ i ].first  ;
          delete localList_[ i ].second ;
        }
      }

      //! returns actual balanceCounter for checkpointing 
      int balanceCounter () const { return balanceCounter_; }
      
      //! do load balance every balanceStep_ step 
      bool loadBalance () 
      {
        // make sure this is only called in single thread mode 
        assert( Fem :: ThreadManager :: singleThreadMode() );

        // get stopwatch 
        Dune::Timer timer ; 
        
        bool changed = false;

        // if balance counter has readed balanceStep do load balance
        const bool callBalance = ( (balanceCounter_ >= balanceStep_) && (balanceStep_ > 0) );

#ifndef NDEBUG
        // make sure load balance is called on every process 
        int willCall = (callBalance) ? 1 : 0;
        const int iCall = willCall;

        // send info from rank 0 to all other 
        grid_.comm().broadcast(&willCall, 1 , 0);

        assert( willCall == iCall );
#endif
        
        // if balance counter has reached balanceStep do load balance
        if( callBalance )
        {
          try {
            // call grids load balance, only implemented in ALUGrid right now
            changed = grid_.loadBalance( dm_ ); 
          }
          catch (...) 
          {
            std::cout << "P[" << grid_.comm().rank() << "] : Cought an exepction during load balance" << std::endl;
            abort();
          }
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

      //! backup internal data 
      void backup() const 
      { 
        Dune::tuple<const int& > value( balanceCounter_ );
        PersistenceManager::backupValue("loadbalancer",value);
      }

      //! retore internal data 
      void restore() 
      {
        Dune::tuple< int& > value( balanceCounter_ );
        PersistenceManager::restoreValue("loadbalancer",value);
      }

      //! add discrete function to data inliner/xtractor list 
      template <class DiscreteFunctionType>
      void addToLoadBalancer(DiscreteFunctionType& df)
      {
        addDiscreteFunction(df);    
      }

      //! add discrete function to data inliner/xtractor list 
      template <class DiscreteFunctionType> 
      void addDiscreteFunction( DiscreteFunctionType& df ) 
      {
        addDiscreteFunction( df, df.defaultLoadBalanceContainsCheck() );
      }

      //! add discrete function to data inliner/xtractor list 
      template <class DiscreteFunctionType, class ContainsCheck > 
      void addDiscreteFunction(DiscreteFunctionType& df, const ContainsCheck& containsCheck ) 
      {
        static_assert( (Conversion< DiscreteFunctionType, IsDiscreteFunction >::exists),
                            "Only valid for discrete functions" );

        const IsDiscreteFunction * fct = &df;

        // if discrete functions is not in list already 
        if( listOfFcts_.find(fct) == listOfFcts_.end() )
        {
          // insert into set 
          listOfFcts_.insert( fct );

          ////////////////////////////
          // data inliners 
          ////////////////////////////
          LocalDataCollectorPairType localPair; 
          DataCollectorPairType collPair; 
          {
            typedef LocalDataInliner<DiscreteFunctionType, ContainsCheck > LocalInlinerType; 
            LocalInlinerType * di = new LocalInlinerType(df, containsCheck );
            localPair.first = di ;
          
            typedef DataCollector<GridType, LocalInlinerType > DataCollectorImp;
            DataCollectorImp* gdi = new DataCollectorImp( grid_, dm_ , *di, di->readWriteInfo() );
            collPair.first = gdi ;
            
            dm_.addDataInliner( *gdi );
          }
         
          ////////////////////////////
          // data xtractors 
          ////////////////////////////
          {
            typedef LocalDataXtractor< DiscreteFunctionType, ContainsCheck > LocalXtractorType; 
            LocalXtractorType * dx = new LocalXtractorType(df, containsCheck );
            localPair.second = dx ;

            typedef DataCollector<GridType,LocalXtractorType> DataCollectorImp;
            DataCollectorImp* gdx = new DataCollectorImp( grid_, dm_ , *dx, dx->readWriteInfo() );
            collPair.second = gdx ;
            
            dm_.addDataXtractor( *gdx );
          }

          // for later removal 
          localList_.push_back( localPair );
          collList_.push_back( collPair );

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
      const int balanceStep_ ;
      // count actual balance call
      int balanceCounter_;

      // list of created local data collectors 
      std::vector< LocalDataCollectorPairType > localList_;
      std::vector< DataCollectorPairType > collList_;

      // list of already added discrete functions 
      std::set< const IsDiscreteFunction * > listOfFcts_;

      // time for last load balance call
      double balanceTime_;
    };

    /** @} end documentation group */

  } // namespace Fem 

} // namespace Dune 
#endif // #ifndef DUNE_FEM_LOADBALANCER_HH
