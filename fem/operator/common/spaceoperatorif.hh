#ifndef DUNE_SPACEOPERATORIF_HH
#define DUNE_SPACEOPERATORIF_HH

//- system includes 
#include <utility>

//-Dune fem includes 
#include <dune/fem/misc/timeutility.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/objpointer.hh>

namespace Dune {

/** \brief SpaceOperatorInterface for Operators of the type 
  L: X --> X where X is a discrete function space
*/
template <class DestinationType>
class SpaceOperatorInterface 
: public Operator< typename DestinationType :: RangeFieldType,
                   typename DestinationType :: RangeFieldType,
                   DestinationType,
                   DestinationType>
{
protected:
  SpaceOperatorInterface() {}

public:
  // destructor 
  virtual ~SpaceOperatorInterface() {}

  // apply operator 
  virtual void operator () (const DestinationType& arg, DestinationType& dest) const = 0;

  // pass time provider to underlying operator 
  virtual void timeProvider(TimeProvider* tp) {}
};

//! only for kepping the pointer 
template <class OperatorType>
class SpaceOperatorPtr
: public SpaceOperatorInterface<typename OperatorType::DestinationType>,
  public ObjPointerStorage                  
{
  typedef typename OperatorType::DestinationType DestinationType;
  OperatorType* op_;
  
  //! copying not allowed 
  SpaceOperatorPtr(const SpaceOperatorPtr& org);
  SpaceOperatorPtr& operator = (const SpaceOperatorPtr& org);

public:
  //! constructor storing pointer 
  SpaceOperatorPtr(OperatorType * op)
    : op_(op)
  {}

  //! destructor deletes operator 
  ~SpaceOperatorPtr()
  {
    // delete operator before destructor of base class is called
    delete op_; op_ = 0;
  }

  //! application operator does nothing here
  void operator () (const DestinationType& arg, DestinationType& dest) const
  {
    // this method should not be called 
    assert(false);
    abort();
  }

  //! return reference to pass 
  OperatorType& pass() 
  { 
    assert( op_ );
    return (*op_); 
  }
};

//! apply wrapper 
template <class OperatorType>
class SpaceOperatorWrapper
: public SpaceOperatorInterface<typename OperatorType::DestinationType>,
  public ObjPointerStorage                  
{
  // operator pointer 
  OperatorType* op_;

  //! copying not allowed
  SpaceOperatorWrapper(const SpaceOperatorWrapper& org);
  SpaceOperatorWrapper& operator = (const SpaceOperatorWrapper& org);
public:
  //! type of Argument and Destination 
  typedef typename OperatorType::DestinationType DestinationType;
  
  //! constructor storing pointer 
  SpaceOperatorWrapper(OperatorType * op)
    : op_(op)
  {}
    
  //! destructor deletes operator 
  ~SpaceOperatorWrapper()
  {
    // delete operator before destructor of base class is called
    delete op_; op_ = 0;
  }

  //! call application operator of internal operator  
  void operator () (const DestinationType& arg, DestinationType& dest) const
  {
    assert( op_ );
    (*op_)(arg,dest);
  }
  
  //! pass timeprovider to internal operator 
  void timeProvider(TimeProvider* tp) 
  { 
    assert( op_ );
    (*op_).timeProvider(tp); 
  }
};

/** \brief CreatePass takes a discrete model and a PassType (like LocalDGPass)
 and creates with the parameter PreviousPass in the method create the
 desired pass. The advantage here is, that no typedefs have to be done.
*/
template <class Model, template <class,class> class PassType>
class CreatePass
{
public:  
  typedef typename Model :: Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename Model :: Traits :: DiscreteFunctionType DestinationType;

protected:
  Model& model_;
  const DiscreteFunctionSpaceType& space_;

public:
  //! constructor 
  //! \param model DiscreteModel 
  //! \param space DiscreteFunctionSpace
  CreatePass(Model& model, const DiscreteFunctionSpaceType& space)
    : model_(model) , space_(space)
  {
  }

  //! copy constructor
  CreatePass(const CreatePass& org)
    : model_(org.model_) , space_(org.space_)
  {
  }

  //! creation method
  template <class PreviousPass>
  SpaceOperatorPtr<PassType<Model,PreviousPass> >*
  create(SpaceOperatorPtr<PreviousPass>* prevObj)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevObj->pass(),space_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    // remember previous object for delete 
    obj->saveObjPointer(prevObj);
    return obj;
  }
  
  //! creation method
  template <class PreviousPass>
  SpaceOperatorPtr<PassType<Model,PreviousPass> >*
  createFirst(PreviousPass& prevPass)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevPass,space_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    return obj;
  }
  
  //! creation method
  template <class PreviousPass>
  SpaceOperatorWrapper<PassType<Model,PreviousPass> >*
  createLast(PreviousPass& prevPass)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevPass,space_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    return obj;
  }
  
  //! last creation method 
  template <class PreviousPass>
  SpaceOperatorWrapper<PassType<Model,PreviousPass> >*
  createLast(SpaceOperatorPtr<PreviousPass>* prevObj)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevObj->pass(),space_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    // remember previous object for delete 
    obj->saveObjPointer(prevObj);
    return obj;
  }
};

template <class Model, template <class,class> class PassType, 
          class SpaceType = typename Model :: Traits :: DiscreteFunctionSpaceType >
class CreateFeaturedPass
{
public:  
  // could also be discrete function 
  typedef SpaceType DiscreteFunctionSpaceType;
  typedef typename Model :: Traits :: DiscreteFunctionType DestinationType;

protected:
  Model& model_;
  DiscreteFunctionSpaceType& space_;
  const std::string paramFile_;

public:
  //! constructor 
  //! \param model DiscreteModel 
  //! \param space DiscreteFunctionSpace
  //! \param paramfile parameter file passes through 
  CreateFeaturedPass(Model& model, DiscreteFunctionSpaceType& space, std::string paramfile = "")
    : model_(model) , space_(space) , paramFile_(paramfile)
  {
  }

  //! copy constructor
  CreateFeaturedPass(const CreateFeaturedPass& org)
    : model_(org.model_) , space_(org.space_) , paramFile_(org.paramFile_)
  {
  }

  //! creation method
  template <class PreviousPass>
  SpaceOperatorPtr<PassType<Model,PreviousPass> >*
  create(SpaceOperatorPtr<PreviousPass>* prevObj)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevObj->pass(),space_,paramFile_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    // remember previous object for delete 
    obj->saveObjPointer(prevObj);
    return obj;
  }
  
  //! creation method
  template <class PreviousPass>
  SpaceOperatorPtr<PassType<Model,PreviousPass> >*
  createFirst(PreviousPass& prevPass)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevPass,space_,paramFile_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    return obj;
  }
  
  //! creation method
  template <class PreviousPass>
  SpaceOperatorWrapper<PassType<Model,PreviousPass> >*
  createLast(PreviousPass& prevPass)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevPass,space_,paramFile_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    return obj;
  }
  
  //! last creation method 
  template <class PreviousPass>
  SpaceOperatorWrapper<PassType<Model,PreviousPass> >*
  createLast(SpaceOperatorPtr<PreviousPass>* prevObj)
  {
    typedef PassType<Model,PreviousPass> RealPassType;
    typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
    // create pass 
    RealPassType* pass = new RealPassType(model_,prevObj->pass(),space_,paramFile_);
    // create pass storage 
    ObjPtrType* obj = new ObjPtrType(pass);
    // remember previous object for delete 
    obj->saveObjPointer(prevObj);
    return obj;
  }
};

struct CreatePassTree
{
  //! create 2 passes 
  template <class FirstModel,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*  
  create(FirstModel& mf,
               LastModel& ml)
  {
    return ml.createLast( mf );
  }

  //! create 3 passes 
  template <class StartModel,
            class Mod1,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               LastModel& mlast)
  {
    return mlast.createLast( m1.createFirst( sm ) );
  }

  //! create 4 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               LastModel& mlast)
  {
    return mlast.createLast( m2.create( m1.createFirst( sm ) ) );
  }

  //! create 5 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class Mod3,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               Mod3& m3,
               LastModel& mlast)
  {
    return
      mlast.createLast( 
        m3.create(
          m2.create(
            m1.createFirst( sm ))));
  }
  
  //! create 6 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class Mod3,
            class Mod4,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               Mod3& m3,
               Mod4& m4,
               LastModel& mlast)
  {
    return
      mlast.createLast( 
              m4.create(
                m3.create(
                  m2.create(
                    m1.createFirst( sm )))) );
  }
  //! create 7 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class Mod3,
            class Mod4,
            class Mod5,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               Mod3& m3,
               Mod4& m4,
               Mod5& m5,
               LastModel& mlast)
  {
    return
      mlast.createLast( 
            m5.create(
              m4.create(
                m3.create(
                  m2.create(
                    m1.createFirst( sm )))) ));
  }
  
  //! create 8 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class Mod3,
            class Mod4,
            class Mod5,
            class Mod6,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               Mod3& m3,
               Mod4& m4,
               Mod5& m5,
               Mod6& m6,
               LastModel& mlast)
  {
    return
      mlast.createLast( 
          m6.create(
            m5.create(
              m4.create(
                m3.create(
                  m2.create(
                    m1.createFirst( sm )))) )));
  }
  //! create 9 passes 
  template <class StartModel,
            class Mod1,
            class Mod2,
            class Mod3,
            class Mod4,
            class Mod5,
            class Mod6,
            class Mod7,
            class LastModel>
  inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
  create(StartModel& sm,
               Mod1& m1,
               Mod2& m2,
               Mod3& m3,
               Mod4& m4,
               Mod5& m5,
               Mod6& m6,
               Mod7& m7,
               LastModel& mlast)
  {
    return
      mlast.createLast( 
        m7.create(
          m6.create(
            m5.create(
              m4.create(
                m3.create(
                  m2.create(
                    m1.createFirst( sm )))) ))));
  }
};

} // end namespace Dune 
#endif
