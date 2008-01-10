#ifndef DUNE_COMMUNICATION_MANAGER_HH
#define DUNE_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_ALUGRID && ALU3DGRID_PARALLEL 
#define USE_CACHED_COMM_MANAGER 
#else 
#ifndef NDEBUG 
#if HAVE_MPI == 0
  #warning "HAVE_MPI == 0, therefore default CommunicationManager is used!"
#elif !ALU3DGRID_PARALLEL 
  #warning "No Parallel ALUGrid found, using default CommunicationManager!"
#endif 
#endif
#endif

#ifdef USE_CACHED_COMM_MANAGER
#include "cachedcommmanager.hh"
#endif

//- Dune-fem includes 
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune { 
  
/** @addtogroup Communication Communication 
    @{
**/
  
  //! \brief Default CommunicationManager class just using the grids communicate
  //! method 
  template <class SpaceImp> 
    //, 
    //        class OperationImp = DFCommunicationOperation::Copy >
  class DefaultCommunicationManager 
  {
    typedef SpaceImp SpaceType; 
    typedef typename SpaceType :: GridPartType GridPartType; 

    // gridPart for communication 
    const GridPartType & gridPart_; 

    const InterfaceType interFace_;
    const CommunicationDirection dir_;
    
    DefaultCommunicationManager(const DefaultCommunicationManager &);
  public:
    //! constructor taking space, but here only storing gridPart for
    //! communication
    DefaultCommunicationManager(const SpaceType & space,
                                const InterfaceType interFace = InteriorBorder_All_Interface ,
                                const CommunicationDirection dir = ForwardCommunication)
      : gridPart_(space.gridPart()) 
      , interFace_( interFace )
      , dir_ ( dir )
    {}

    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) 
    {
      // if serial run, just return   
      if(gridPart_.grid().comm().size() <= 1) return;
     
      // get data handler type from space  
      typedef typename SpaceType :: template CommDataHandle<DiscreteFunctionType> :: Type DataHandleType;
      DataHandleType dataHandle = df.dataHandle();

      // communicate data 
      gridPart_.communicate( dataHandle, interFace_ , dir_ );
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
    //! constructor taking space, but here only storing gridPart for
    //! communication
    CommunicationManager(const SpaceImp & space) 
      : BaseType(space) 
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
      virtual void exchange () = 0;
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
      void exchange () 
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
    void exchange() 
    {
      typedef CommObjListType :: iterator iterator; 
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
