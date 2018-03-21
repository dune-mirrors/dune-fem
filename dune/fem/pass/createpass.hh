#ifndef DUNE_FEM_CREATEPASS_HH
#define DUNE_FEM_CREATEPASS_HH

#include <memory>

#include <dune/fem/common/memory.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/pass/common/pass.hh>

namespace Dune
{

  namespace Fem
  {

    /** @ingroup Pass
      \brief CreatePass takes a discrete model and a PassType (like LocalDGPass)
      and creates with the parameter PreviousPass in the method create the
      desired pass. The advantage here is, that no typedefs have to be done.

      To generate a pass tree one only has to use the following example code:
      @code
      // diffusion pass
      typedef CreatePass<DiscreteModel1Type,LocalDGPass> Pass1Type;
      Pass1Type pass1 ( discreteModel1_, space1_ );

      // advection pass
      typedef CreatePass<DiscreteModel2Type,LocalDGPass> Pass2Type;
      Pass2Type pass2 ( discreteModel2_, space2_ );

      // create pass tree and return pointer to resulting
      // operator satisfying the SpaceOperatorInterface.
      SpaceOperatorInterface<DestinationType>* passTree
        = CreatePassTree::create( pass1 , pass2 );
      @endcode
    */
    template< class Model , template <class,class,int> class PassType , int pId  = -1 >
    class CreatePass
    {
    public:
      //! forward pass id
      enum { passId = pId };

      //! type of discrete function space
      typedef typename Model :: Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! destination type
      typedef typename Model :: Traits :: DiscreteFunctionType DestinationType;

      //! type of space operator
      typedef SpaceOperatorInterface<DestinationType> SpaceOperatorIFType;
    protected:
      Model& model_;
      std::shared_ptr< const DiscreteFunctionSpaceType > space_;
      SpaceOperatorIFType* passPointer_ = nullptr;

    public:
      //! constructor
      //! \param model DiscreteModel
      //! \param space DiscreteFunctionSpace
      CreatePass ( Model &model, const DiscreteFunctionSpaceType &space )
        : model_( model ), space_( referenceToSharedPtr( space ) )
      {}

      //! constructor
      //! \param model DiscreteModel (or discrete function)
      //! \param space DiscreteFunctionSpace
      CreatePass ( const Model &model, const DiscreteFunctionSpaceType &space )
        : CreatePass( const_cast< Model & >( model ), space )
      {}

      //! copy constructor
      CreatePass ( const CreatePass & ) = default;

      //! creation method
      template <class PreviousPass>
      SpaceOperatorPtr< PassType< Model , PreviousPass , passId > >*
      create(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType< Model , PreviousPass , passId > RealPassType;
        typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(model_,prevObj->pass(), *space_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass);

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! last creation method
      template <class PreviousPass>
      SpaceOperatorWrapper< PassType< Model , PreviousPass , passId > >*
      createLast(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType< Model , PreviousPass , passId > RealPassType;
        typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(model_,prevObj->pass(), *space_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass);

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! return pointer to space operator if
      SpaceOperatorIFType* pass()
      {
        assert( passPointer_ );
        return passPointer_;
      }

      //! return pointer to destination
      const DestinationType* destination() const
      {
        assert( passPointer_ );
        return passPointer_->destination();
      }
    };


    //! DiscreteModelWrapper to combine DiscreteModel and Selector
    template <class DiscreteModelImp, class SelectorImp>
    class DiscreteModelWrapper
    : public ObjPointerStorage,
      public DiscreteModelImp
    {
      // type of base class
      typedef DiscreteModelImp BaseType;
    public:
      //! exporting given type of selector
      typedef SelectorImp SelectorType;
      //! constructor calling the copy constructor of the base type
      DiscreteModelWrapper(const BaseType& base)
       : BaseType(base) // calling copy constructor of base
      {}

      //! copy constructor
      DiscreteModelWrapper(const DiscreteModelWrapper& other)
       : BaseType(other) // calling copy constructor of base
      {}
    };

    //! create pass with previous unknown selector
    template <class Model, class SelectorImp, template <class,class> class PassType>
    class CreateSelectedPass
    {
    public:
      //! type of discrete function space
      typedef typename Model :: Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! destination type
      typedef typename Model :: Traits :: DiscreteFunctionType DestinationType;

      //! type of space operator
      typedef SpaceOperatorInterface<DestinationType> SpaceOperatorIFType;

      //! type of discrete model
      typedef DiscreteModelWrapper<Model,SelectorImp> DiscreteModelType;
    protected:
      DiscreteModelType* model_;
      std::shared_ptr< const DiscreteFunctionSpaceType > space_;
      SpaceOperatorIFType *passPointer_ = nullptr;
      bool owner_ = true;

    public:
      //! constructor
      //! \param model DiscreteModel
      //! \param space DiscreteFunctionSpace
      CreateSelectedPass(Model& model, const DiscreteFunctionSpaceType& space)
        : model_(new DiscreteModelType(model)), space_( referenceToSharedPtr( space ) )
      {}

      //! copy constructor
      CreateSelectedPass(const CreateSelectedPass& org)
        : model_(new DiscreteModelType(org.model_)),
          space_(org.space_),
          passPointer_( org.passPointer_ )
      {}

      //! destructor deleting model if still owner
      ~CreateSelectedPass()
      {
        if( owner_ ) delete model_;
        model_ = nullptr;
      }

      //! creation method
      template <class PreviousPass>
      SpaceOperatorPtr<PassType<DiscreteModelType,PreviousPass> >*
      create(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType<DiscreteModelType,PreviousPass> RealPassType;
        typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(*model_,prevObj->pass(), *space_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass, model_ );

        // don't delete model anymore
        owner_ = false;

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! last creation method
      template <class PreviousPass>
      SpaceOperatorWrapper<PassType<DiscreteModelType,PreviousPass> >*
      createLast(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType<DiscreteModelType,PreviousPass> RealPassType;
        typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(*model_,prevObj->pass(), *space_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass, model_);

        // don't delete model anymore
        owner_ = false;

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! return pointer to space operator if
      SpaceOperatorIFType* pass()
      {
        assert( passPointer_ );
        return passPointer_;
      }

      //! return pointer to destination
      const DestinationType* destination() const
      {
        assert( passPointer_ );
        return passPointer_->destination();
      }
    };

    /**
     \brief CreateFeaturedPass takes a discrete model and a PassType (like LocalDGEllliptPass)
     and creates with the parameter PreviousPass in the method create the
     desired pass. The advantage here is, that no typedefs have to be done.
    */
    template <class Model, template <class,class,int> class PassType,
              class SpaceType, // = typename Model :: Traits :: DiscreteFunctionSpaceType ,
              int pId
              >
    class CreateFeaturedPass
    {
    public:
      //! forward pass id
      enum { passId = pId };

      //! type of discrete functions space
      typedef SpaceType DiscreteFunctionSpaceType;
      //! destination type
      typedef typename Model :: Traits :: DiscreteFunctionType DestinationType;

      //! type of space operator
      typedef SpaceOperatorInterface<DestinationType> SpaceOperatorIFType;
    protected:
      Model& model_;
      std::shared_ptr< DiscreteFunctionSpaceType > space_;
      const std::string paramFile_;
      SpaceOperatorIFType* passPointer_ = nullptr;

    public:
      //! constructor
      //! \param model DiscreteModel
      //! \param space DiscreteFunctionSpace
      //! \param paramfile parameter file passes through
      CreateFeaturedPass(Model& model, DiscreteFunctionSpaceType& space, std::string paramfile = "")
        : model_(model),
          space_( referenceToSharedPtr( space) ),
          paramFile_(paramfile)
      {}

      //! copy constructor
      CreateFeaturedPass ( const CreateFeaturedPass & ) = default;

      //! creation method
      template <class PreviousPass>
      SpaceOperatorPtr< PassType<Model,PreviousPass,passId> >*
      create(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType<Model,PreviousPass,passId> RealPassType;
        typedef SpaceOperatorPtr<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(model_,prevObj->pass(), *space_,paramFile_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass);

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! last creation method
      template <class PreviousPass>
      SpaceOperatorWrapper< PassType<Model,PreviousPass,passId> >*
      createLast(SpaceOperatorStorage<PreviousPass>* prevObj)
      {
        typedef PassType<Model,PreviousPass,passId> RealPassType;
        typedef SpaceOperatorWrapper<RealPassType> ObjPtrType;
        // create pass
        RealPassType* pass = new RealPassType(model_,prevObj->pass(), *space_,paramFile_);

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(pass);

        // remember pass obj
        passPointer_ = obj;

        // remember previous object for delete
        obj->saveObjPointer(prevObj);
        return obj;
      }

      //! return pointer to destination
      const DestinationType* destination() const
      {
        assert( passPointer_ );
        return passPointer_->destination();
      }
    };

    /** \brief create pass tree from given list of discrete models
        the passId is deliviered to the start pass and stands for the global variable
        calculated by this pass
     */
    template <int startPassId = -1>
    class CreatePassTree
    {
    protected:
      //! method that creates first pass
      template <class DestinationType>
      inline static SpaceOperatorStorage< StartPass<DestinationType,startPassId> >* createStartPass()
      {
        typedef StartPass<DestinationType, startPassId> StartPassType;
        typedef SpaceOperatorStorage<StartPassType> ObjPtrType;

        // create start pass
        StartPassType* startPass = new StartPassType ();

        // create pass storage
        ObjPtrType* obj = new ObjPtrType(startPass);

        // return obj
        return obj;
      }

    public:
      //! create 1 pass
      template <class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(LastModel& ml)
      {
        return ml.createLast( createStartPass<typename LastModel :: DestinationType> () );
      }

      //! create 2 passes
      template <class FirstModel,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(FirstModel& mf,
             LastModel& ml)
      {
        return ml.createLast( mf.create( createStartPass<typename LastModel :: DestinationType>() ) );
      }

      //! create 3 passes
      template <class Mod0,
                class Mod1,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
             Mod1& m1,
             LastModel& mlast)
      {
        return mlast.createLast(
              m1.create(
                m0.create( createStartPass<typename LastModel :: DestinationType> () )
                       )
                               );
      }

      //! create 4 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
             Mod1& m1,
             Mod2& m2,
             LastModel& mlast)
      {
        return mlast.createLast( m2.create( m1.create(
                m0.create( createStartPass<typename LastModel :: DestinationType> () )
                  ) ) );
      }

      //! create 5 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class Mod3,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
             Mod1& m1,
             Mod2& m2,
             Mod3& m3,
             LastModel& mlast)
      {
        return
          mlast.createLast(
            m3.create(
              m2.create(
                m1.create(
                  m0.create( createStartPass<typename LastModel :: DestinationType> () )
                  ))));
      }

      //! create 6 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class Mod3,
                class Mod4,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
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
                        m1.create(
                          m0.create( createStartPass<typename LastModel :: DestinationType> () )
                          )))) );
      }
      //! create 7 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class Mod3,
                class Mod4,
                class Mod5,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
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
                        m1.create(
                          m0.create( createStartPass<typename LastModel :: DestinationType> () )
                          )))) ));
      }

      //! create 8 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class Mod3,
                class Mod4,
                class Mod5,
                class Mod6,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
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
                        m1.create(
                          m0.create( createStartPass<typename LastModel :: DestinationType> () )
                          )))) )));
      }
      //! create 9 passes
      template <class Mod0,
                class Mod1,
                class Mod2,
                class Mod3,
                class Mod4,
                class Mod5,
                class Mod6,
                class Mod7,
                class LastModel>
      inline static SpaceOperatorInterface<typename LastModel :: DestinationType>*
      create(Mod0& m0,
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
                        m1.create(
                          m0.create( createStartPass<typename LastModel :: DestinationType> () )
                          )))) ))));
      }
    };

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_CREATEPASS_HH
