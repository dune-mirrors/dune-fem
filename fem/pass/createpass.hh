#ifndef DUNE_CREATEPASS_HH
#define DUNE_CREATEPASS_HH

#include <dune/fem/operator/common/spaceoperatorif.hh>

namespace Dune {

/** @ingroup Pass 
  \brief CreatePass takes a discrete model and a PassType (like LocalDGPass)
  and creates with the parameter PreviousPass in the method create the
  desired pass. The advantage here is, that no typedefs have to be done.

  To generate a pass tree one  only has to use the following example code:
  @code 
  // diffusion pass
  typedef CreatePass<DiscreteModel1Type,LocalDGPass> Pass1Type;
  Pass1Type pass1 ( problem1_, space1_ );

  // advection pass
  typedef CreatePass<DiscreteModel2Type,LocalDGPass> Pass2Type;
  Pass2Type pass2 ( problem2_, space2_ );

  // create pass tree and return pointer to resulting 
  // operator satisfying the SpaceOperatorInterface.
  SpaceOperatorInterface<DestinationType>* passTree 
    = CreatePassTree::create( pass0_ , pass1 , pass2 );
  @endcode
*/
template <class Model, template <class,class> class PassType>
class CreatePass
{
public:  
  //! type of discrete function space
  typedef typename Model :: Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  //! destination type 
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

/** 
 \brief CreateFeaturedPass takes a discrete model and a PassType (like LocalDGEllliptPass)
 and creates with the parameter PreviousPass in the method create the
 desired pass. The advantage here is, that no typedefs have to be done.
*/
template <class Model, template <class,class> class PassType, 
          class SpaceType = typename Model :: Traits :: DiscreteFunctionSpaceType >
class CreateFeaturedPass
{
public:  
  //! type of discrete functions space 
  typedef SpaceType DiscreteFunctionSpaceType;
  //! destination type 
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

/** \brief create pass tree from given list of discrete models 
 */
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
