#ifndef DUNE_COMMUNICATION_MANAGER_HH
#define DUNE_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/space/common/commoperations.hh>

// include ALUGrid to check whether the 
// parallel version is avaiable 
#if defined ENABLE_ALUGRID 
#include <dune/grid/alugrid.hh>
#endif

// default is: enabled 
#ifndef WANT_CACHED_COMM_MANAGER 
#define WANT_CACHED_COMM_MANAGER 1
#endif

#if ALU3DGRID_PARALLEL && WANT_CACHED_COMM_MANAGER
#define USE_CACHED_COMM_MANAGER 
#else 
#ifndef NDEBUG 
#if HAVE_MPI == 0
  #warning "HAVE_MPI == 0, therefore default CommunicationManager is used!"
#elif !ALU3DGRID_PARALLEL 
  #warning "No Parallel ALUGrid found, using default CommunicationManager!"
#elif ! WANT_CACHED_COMM_MANAGER 
  #warning "CachedCommunication Manager disabled by WANT_CACHED_COMM_MANAGER=0!"
#endif 
#endif
#endif

#undef WANT_CACHED_COMM_MANAGER 

#ifdef USE_CACHED_COMM_MANAGER
#include "cachedcommmanager.hh"
#endif

//- Dune-fem includes 
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune
{
  
/** @addtogroup Communication Communication 
    @{
**/

  /** \class DefaultCommunicationManager
   *  \ingroup Communication
   *  \brief default communication manager using just the grids communicate
   *         method
   */
  template< class Space >
  class DefaultCommunicationManager 
  {
  public:
    typedef Space SpaceType; 

  protected:
    const SpaceType& space_;

    const InterfaceType interface_;
    const CommunicationDirection dir_;

    mutable double exchangeTime_;
    
  public:
    //! constructor taking space and communication interface/direction
    explicit DefaultCommunicationManager
      ( const SpaceType &space,
        const InterfaceType interface,
        const CommunicationDirection dir)
    : space_( space ),
      interface_( interface ),
      dir_ ( dir ),
      exchangeTime_(0.0)
    {}

  private:
    // prohibit copying
    DefaultCommunicationManager ( const DefaultCommunicationManager & );
    
  public:
    /** \brief return communication interface */
    InterfaceType communicationInterface() const {
      return interface_;
    }

    /** \brief return communication direction */
    CommunicationDirection communicationDirection() const
    {
      return dir_;
    }

    /** \brief return time needed for last build 

        \return time needed for last build of caches (if needed)
    */
    double buildTime() const { return 0.0; }

    /** \brief return time needed for last exchange of data  

        \return time needed for last exchange of data 
    */
    double exchangeTime() const { return exchangeTime_; }

    /** \brief exchange data for a discrete function using the copy operation
     *  
     *  \param  discreteFunction  discrete function to communicate
     */
    template< class DiscreteFunction >
    inline void exchange ( DiscreteFunction &discreteFunction ) const
    {
      // get type of default operation 
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: 
        template CommDataHandle< DiscreteFunction > :: OperationType DefaultOperationType;
      
      exchange( discreteFunction, (DefaultOperationType *) 0 );
    }

    /** \brief exchange data for a discrete function using the given operation
     *
     *  The used operation is derived from the type of the op-pointer. The
     *  actual pointer is not used.
     *  
     *  \param      discreteFunction  discrete function to communicate
     *  \param[in]  operation         a (phony) pointer to an operation
     */
    template< class DiscreteFunction, class Operation >
    inline void exchange ( DiscreteFunction &discreteFunction,
                           const Operation *operation ) const
    {
      // get type of data handle from the discrete function space
      typedef typename DiscreteFunction
        :: template CommDataHandle< Operation > :: Type
        DataHandleType;
      
      // on serial runs: do nothing
      if( space_.grid().comm().size() <= 1 )
        return;

      // get stopwatch 
      Timer exchangeT;
    
      // communicate data
      DataHandleType dataHandle = discreteFunction.dataHandle( operation );
      space_.gridPart().communicate( dataHandle, interface_ , dir_ );

      // store time 
      exchangeTime_ = exchangeT.elapsed();
    }
  };


  
#ifndef USE_CACHED_COMM_MANAGER
  // if no ALUGrid found, supply default implementation 
  //! \brief use Default CommunicationManager as Communication Manager 
  template <class SpaceImp> 
  class CommunicationManager 
  : public DefaultCommunicationManager<SpaceImp> 
  {
    typedef DefaultCommunicationManager<SpaceImp> BaseType;
    CommunicationManager(const CommunicationManager &);
  public:
    //! constructor taking space and communication interface/direction
    CommunicationManager(const SpaceImp & space,  
        const InterfaceType interface,
        const CommunicationDirection dir)
      : BaseType(space,interface,dir) 
    {}
    //! constructor taking space
    CommunicationManager(const SpaceImp & space)
      : BaseType(space,
                 space.communicationInterface(),
                 space.communicationDirection() )
    {}
  };



  //! Proxy class to DependencyCache which is singleton per space 
  class CommunicationManagerList  
  {
    //! communicated object interface 
    class DiscreteFunctionCommunicatorInterface  
    {
    protected:
      DiscreteFunctionCommunicatorInterface () {}
    public:
      virtual ~DiscreteFunctionCommunicatorInterface () {}
      virtual void exchange () const = 0;
    };
    
    //! communicated object implementation  
    template <class DiscreteFunctionImp>
    class DiscreteFunctionCommunicator 
    : public DiscreteFunctionCommunicatorInterface 
    {
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      
      typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
    
      DiscreteFunctionType& df_;
      CommunicationManagerType comm_;
    public:  
      //! constructor taking disctete function 
      DiscreteFunctionCommunicator(DiscreteFunctionType& df) 
        : df_(df), comm_(df_.space())
      {
      }

      // exchange discrete function 
      void exchange () const
      {
        comm_.exchange(df_);
      }
    };

    typedef DiscreteFunctionCommunicatorInterface CommObjIFType;
    typedef std::list < DiscreteFunctionCommunicatorInterface * > CommObjListType;
    CommObjListType objList_;

    CommunicationManagerList(const CommunicationManagerList&); 
  public: 
    //! constructor
    template <class CombinedObjectType>
    CommunicationManagerList(CombinedObjectType& cObj) 
    {
      cObj.addToList(*this);
    }

    //! remove object comm
    ~CommunicationManagerList() 
    {
      // delete all entries 
      while( ! objList_.size() == 0 )
      {
        CommObjIFType * obj = objList_.back();
        objList_.pop_back();
        delete obj;
      }
    }

    //! add discrete function to communication list 
    template <class DiscreteFunctionImp> 
    void addToList(DiscreteFunctionImp &df)
    {
      typedef DiscreteFunctionCommunicator<DiscreteFunctionImp> CommObjType;
      CommObjType* obj = new CommObjType(df);
      objList_.push_back(obj);
    }

    //! exchange discrete function to all procs we share data with 
    //! by using given OperationImp when receiving data from other procs 
    void exchange() const
    {
      typedef CommObjListType :: const_iterator iterator; 
      {
        iterator end = objList_.end();
        for(iterator it = objList_.begin(); it != end; ++it) 
        {
          (*it)->exchange();
        }
      }
    }
  };

  // end toggle AULGrid yes/no
#endif 
  //@}
  
} // end namespace Dune

#endif
