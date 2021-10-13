#ifndef DUNE_FEM_LOADBALANCER_HH
#define DUNE_FEM_LOADBALANCER_HH

#include <cassert>
#include <iostream>
#include <set>
#include <type_traits>
#include <vector>

#include <dune/common/timer.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/datacollector.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/commoperations.hh>

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
      /** \brief constructor of LoadBalancer **/
      template< class RestrictProlongOperator >
      LoadBalancer ( GridType &grid, RestrictProlongOperator &rpOp )
      : grid_( grid ),
        dm_ ( DofManagerType::instance( grid_ ) ),
        localList_(),
        collList_(),
        commList_(rpOp),
        balanceTime_( 0.0 )
      {
        rpOp.addToLoadBalancer( *this );
      }

      explicit LoadBalancer ( GridType &grid )
      : grid_( grid ),
        dm_ ( DofManagerType::instance( grid_ ) ),
        localList_(),
        collList_(),
        commList_(),
        balanceTime_( 0.0 )
      {}

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

      void communicate () const
      {
        // if overlap or ghost elements are available
        // these need to be synchronized here
        const auto gv = grid_.leafGridView();
        if( (gv.overlapSize( 0 ) > 0) || (gv.ghostSize( 0 ) > 0) )
        {
          // exchange all modified data
          // this also rebuilds the dependecy cache of the
          // cached communication manager if used
          commList_.exchange();
        }
#ifndef NDEBUG
        // make sure every process is on the same page
        gv.comm().barrier();
#endif // #ifndef NDEBUG
      }

      //! do load balance
      bool loadBalance ()
      {
        bool changed = false;

        // if only one core don't do anything
        if( grid_.comm().size() <= 1 )
          return changed;

        // make sure this is only called in single thread mode
        if( ! Fem :: MPIManager :: singleThreadMode() )
        {
          assert( Fem :: MPIManager :: singleThreadMode() );
          DUNE_THROW(InvalidStateException,"LoadBalancer::loadBalance::adapt: only call in single thread mode!");
        }

        // get stopwatch
        Dune::Timer timer ;


        try {
          // call grids load balance, only implemented in ALUGrid right now
          changed = grid_.loadBalance( dm_ );
        }
        catch (...)
        {
          std::cout << "P[" << grid_.comm().rank() << "] : Caught an exepction during load balance" << std::endl;
          abort();
        }

        // get time
        balanceTime_ = timer.elapsed();

        // restore data consistency
        communicate();

        return changed;
      }

      /** @copydoc LoadBalancerInterface::loadBalanceTime */
      virtual double loadBalanceTime() const
      {
        return balanceTime_;
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
        static_assert( std::is_convertible< DiscreteFunctionType, IsDiscreteFunction >::value,
                       "Only valid for discrete functions" );

        //////////////////////////////////////////////////////////
        //
        //  Note: DiscreteFunctionType here can also be
        //        FemPy::DiscreteFunctionList for
        //        python adaptation and load balance
        //
        //////////////////////////////////////////////////////////

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

      // list of created local data collectors
      std::vector< LocalDataCollectorPairType > localList_;
      std::vector< DataCollectorPairType > collList_;

      // list of already added discrete functions
      std::set< const IsDiscreteFunction * > listOfFcts_;

      mutable CommunicationManagerList commList_;

      // time for last load balance call
      double balanceTime_;
    };

    /** @} end documentation group */

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_LOADBALANCER_HH
