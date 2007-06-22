#ifndef DUNE_DISCRETEFUNCTION_HH
#define DUNE_DISCRETEFUNCTION_HH

//- system includes 
#include <string>

//- Dune inlcudes 
#include <dune/grid/common/grid.hh>

//- local includes 
#include "function.hh"
#include <dune/fem/space/common/discretefunctionspace.hh>
#include "dofiterator.hh"
#include "localfunctionwrapper.hh"

namespace Dune{


  /** @defgroup DiscreteFunction DiscreteFunction
      @ingroup FunctionCommon
      The DiscreteFunction is responsible for the dof storage. This can be
      done in various ways an is left to the user. The user has to derive his
      own implementation from the DiscreteFunctionDefault class. If some of
      the implementations in the default class are for ineffecient for the
      dof storage in the derived class these functions can be overloaded.
  
      @{
  */

  //! empty object for Conversion Utility 
  class IsDiscreteFunction
  {
  };
  
  //! empty object for Conversion Utility 
  class HasLocalFunction
  {
  };
  
  //************************************************************************
  //
  //  --DiscreteFunctionInterface
  //
  //! This is the minimal interface of a discrete function which has to be
  //! implemented. It contains a local function and a dof iterator which can 
  //! iterate over all dofs of one level. Via the method access the local
  //! dofs and basis functions can be accessed for a given entity.
  //! The DOF-Iterators are STL-like Iterators, i.e. they can be dereferenced
  //! giving the corresponding DOF.
  //! 
  //************************************************************************
  template<class DiscreteFunctionTraits>
  class DiscreteFunctionInterface : 
    public IsDiscreteFunction , 
    public HasLocalFunction , 
    public Function<typename DiscreteFunctionTraits::DiscreteFunctionSpaceType,
                    DiscreteFunctionInterface<DiscreteFunctionTraits> > 
  {
  public:
    //- Typedefs and enums

    //! types that we sometimes need outside 
    typedef Function<
      typename DiscreteFunctionTraits::DiscreteFunctionSpaceType,
      DiscreteFunctionInterface<DiscreteFunctionTraits> 
    > FunctionType;

    typedef typename DiscreteFunctionTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    //! Domain vector
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! Range vector
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  
    //! Domain field type (usually a float type)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    //! Range field type (usually a float type)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  
    //! Type of the underlying grid
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    //! Type of the discrete function implementation
    typedef typename DiscreteFunctionTraits::DiscreteFunctionType DiscreteFunctionType;

    //! Type of the local function implementation
    typedef typename DiscreteFunctionTraits::LocalFunctionType LocalFunctionType;

    //! Type of the dof iterator used in the discrete function implementation.
    typedef typename DiscreteFunctionTraits::DofIteratorType DofIteratorType;

    //! Type of the constantdof iterator used in the discrete function implementation
    typedef typename DiscreteFunctionTraits::ConstDofIteratorType ConstDofIteratorType;

  public:
    //- Public Methods

    //! Constructor
    //! Needs to be called in derived classes
    DiscreteFunctionInterface (const DiscreteFunctionSpaceType& f) 
      : FunctionType ( f ) {}

    //! Name of the discrete function
    std::string name() const {
      return asImp().name();
    }

    //! The size of the discrete function
    int size() const {
      return asImp().size();
    }

    //! the implementation of an iterator to iterate efficient 
    //! over all dofs of a discrete function 
    DofIteratorType dbegin () 
    {
      return asImp().dbegin ();
    }

    //! the implementation of an iterator to iterate efficient 
    //! over all dofs of a discrete function 
    DofIteratorType dend () 
    {
      return asImp().dend ();
    }


    //! const version of dbegin 
    ConstDofIteratorType dbegin () const  
    {
      return asImp().dbegin ();
    };

    //! const version of dend 
    ConstDofIteratorType dend () const  
    {
      return asImp().dend ();
    };

    //! return local function for given entity
    template <class EntityType>
    LocalFunctionType localFunction(const EntityType& en) const {
      return asImp().localFunction(en);
    }

  private:
    // Barton-Nackman trick 
    DiscreteFunctionType& asImp() 
    { 
      return static_cast<DiscreteFunctionType&>(*this); 
    }

    //! const version of asImp 
    const DiscreteFunctionType &asImp() const 
    { 
      return static_cast<const DiscreteFunctionType&>(*this); 
    }
  };
  //*************************************************************************
  //
  //  --DiscreteFunctionDefault
  //  
  //! Default implementation of the discrete function. This class provides 
  //! is responsible for the dof storage. Different implementations of the
  //! discrete function use different dof storage. 
  //! The default implementation provides +=, -= and so on operators and 
  //! a DofIterator access, which can run over all dofs in an efficient way. 
  //! Furthermore with an entity you can access a local function to evaluate
  //! the discrete function by multiplying the dofs and the basefunctions. 
  //! 
  //*************************************************************************
  template <class DiscreteFunctionTraits>
  class DiscreteFunctionDefault : 
    public DiscreteFunctionInterface<DiscreteFunctionTraits> 
  { 

    typedef DiscreteFunctionInterface< 
      DiscreteFunctionTraits
    > DiscreteFunctionInterfaceType;

    typedef DiscreteFunctionDefault<
      DiscreteFunctionTraits
    > DiscreteFunctionDefaultType;

    enum { myId_ = 0 };  
  
  public:
    typedef typename DiscreteFunctionTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Type of domain vector
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! Type of range vector
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  
    //! Type of domain field (usually a float type)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    //! Type of range field (usually a float type)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

    //! Type of the discrete function (Barton-Nackman parameter)
    typedef typename DiscreteFunctionTraits::DiscreteFunctionType DiscreteFunctionType;
  
    //! Type of the local function
    typedef typename DiscreteFunctionTraits::LocalFunctionType LocalFunctionType;

    //! Type of the local function implementation 
    typedef typename DiscreteFunctionTraits::LocalFunctionImp LocalFunctionImp;
    //! Type of the dof iterator
    typedef typename DiscreteFunctionTraits::DofIteratorType DofIteratorType;
  
    //! Type of the const dof iterator
    typedef typename DiscreteFunctionTraits::ConstDofIteratorType ConstDofIteratorType;

    typedef LocalFunctionStorage < DiscreteFunctionDefaultType > LocalFunctionStorageType;
    friend class LocalFunctionStorage < DiscreteFunctionDefaultType >;
    friend class LocalFunctionWrapper < DiscreteFunctionType >;
  public:
    //- Methods
    //! pass the function space to the interface class
    DiscreteFunctionDefault (const DiscreteFunctionSpaceType & f ) :
      DiscreteFunctionInterfaceType ( f ) , lfStorage_ (*this) {}

    //! Continuous data
    bool continuous() const DUNE_DEPRECATED {
      return this->functionSpace_.continuous();
    }

    //! print all dof values of this function to stream 
    void print(std::ostream & ) const ;
    
    //! Set all elements to zero
    void clear();

    //! daxpy operation
    void addScaled(const DiscreteFunctionType& g, const RangeFieldType& c);

    //! Evaluate a scalar product of the dofs of two DiscreteFunctions
    //! on the top level of the underlying grid
    RangeFieldType scalarProductDofs(const DiscreteFunctionType& g) const;

    //! Assignment
    virtual DiscreteFunctionDefaultType& assign(const MappingType& g);

    //! Addition of g to discrete function 
    virtual DiscreteFunctionDefaultType& operator += (const MappingType& g);

    //! subtract g from discrete function 
    virtual DiscreteFunctionDefaultType& operator -= (const MappingType &g);
 
    //! multiply with scalar 
    virtual DiscreteFunctionDefaultType& operator *=(const RangeFieldType &scalar);

    //! Division by a scalar
    virtual DiscreteFunctionDefaultType& operator /= (const RangeFieldType &scalar);

    //! evaluate Function f  
    //! \param arg: global coordinate
    //! \param dest: f(arg)
    void evaluate (const DomainType & arg, RangeType & dest) const {
      // Die a horrible death! Never call that one...
      assert(false); abort();
    }

    //! evaluate function and derivatives (just dies)
    //! \param arg: global coordinate
    //! \param dest: f(arg)
    template <int derivation>
    void evaluate  ( const FieldVector<deriType, derivation> &diffVariable, 
                     const DomainType& arg, RangeType & dest) const { 
      // Die a horrible death! Never call that one...
      assert(false); abort();
    }

protected: 
  //this methods are used by the LocalFunctionStorage class 

  //! return pointer to local function implementation 
  LocalFunctionImp* newLocalFunctionObject() const {
    return asImp().newLocalFunctionObject();
  }

  //! return reference for local function storage  
  LocalFunctionStorageType& localFunctionStorage() const { 
    return lfStorage_; 
  }

  // the local function storage stack 
  mutable LocalFunctionStorageType lfStorage_;

private:
    // Barton-Nackman trick 
    DiscreteFunctionType &asImp() 
    { 
      return static_cast<DiscreteFunctionType&>(*this); 
    }
    const DiscreteFunctionType &asImp() const 
    { 
      return static_cast<const DiscreteFunctionType&>(*this); 
    }
  }; // end class DiscreteFunctionDefault 

  template <class FunctionImp, class GridPartImp>
  class DiscreteFunctionAdapter;

  template <class FunctionImp, class GridPartImp> 
  struct DiscreteFunctionAdapterTraits 
  {
    typedef typename FunctionImp :: FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscreteFunctionSpaceAdapter<FunctionSpaceType,GridPartImp> DiscreteFunctionSpaceType;

    typedef GridPartImp GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    //! type of iterator 
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 
    //! type of IndexSet 
    typedef typename GridPartType :: IndexSetType IndexSetType; 

    typedef DiscreteFunctionAdapter<FunctionImp,GridPartImp> DiscreteFunctionType;
  };

  //! adapter to provide local function for a function
  template <class FunctionImp, class GridPartImp>
  class DiscreteFunctionAdapter 
  : public HasLocalFunction , 
    public Function<DiscreteFunctionSpaceAdapter<typename FunctionImp :: FunctionSpaceType,GridPartImp>, 
                    DiscreteFunctionAdapter<FunctionImp,GridPartImp> >
  {
    typedef Function<DiscreteFunctionSpaceAdapter<typename FunctionImp :: FunctionSpaceType,GridPartImp>, 
                     DiscreteFunctionAdapter<FunctionImp,GridPartImp> > BaseType;
  public:  
    typedef GridPartImp GridPartType;
                       
    //! type of traits 
    typedef DiscreteFunctionAdapterTraits<FunctionImp,GridPartImp> Traits;
    //! type of this pointer 
    typedef DiscreteFunctionAdapter<FunctionImp,GridPartImp> ThisType;

    // type of function 
    typedef FunctionImp FunctionType;

    // type of discrete function space
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType; 

    //! type of discrete function space 
    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    //! type of grid 
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
    //! domain type (from function space)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
    //! range type (from function space)
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
    //! jacobian type (from function space)
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    //! type of codim 0 entity
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 

    private:
    class LocalFunction
    {
      //! type of geometry 
      typedef typename EntityType :: Geometry GeometryImp;
    public:  
      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
      //! jacobian type (from function space)
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

      //! constructor initializing local function 
      LocalFunction(const EntityType& en, const ThisType& a)
        : function_(a.function_) 
        , geometry_(&(en.geometry())) 
        , global_(0)
        , result_(0)
        , quadId_(-1)
        , qp_(-1)
      {}

      LocalFunction(const ThisType& a)
        : function_(a.function_) 
        , geometry_(0) 
        , global_(0)
        , result_(0)
        , quadId_(-1)
        , qp_(-1)
      {}

      //! copy constructor 
      LocalFunction(const LocalFunction& org) 
        : function_(org.function_) 
        , geometry_(org.geometry_)  
        , global_(org.global_)
        , result_(org.result_)
        , quadId_(org.quadId_)
        , qp_(org.qp_)
      {}

      //! evaluate local function 
      void evaluate(const DomainType& local, RangeType& result) const
      {
        global_ = geometry_->global(local);
        function_.evaluate(global_,result);
      }

      //! evaluate local function 
      template <class QuadratureType>
      void evaluate(const QuadratureType& quad,
                    const int quadPoint, 
                    RangeType& result) const 
      {
        if(quadId_ != (int) quad.id() && qp_ != quadPoint) 
        {
          // evaluate function and cache value 
          global_ = geometry_->global(quad.point(quadPoint));
          function_.evaluate(global_,result_);
          // remember cached point 
          quadId_ = quad.id();
          qp_ = quadPoint;
        }
        result = result_;
      }

      //! jacobian of local function 
      void jacobian(const DomainType& local, JacobianRangeType& result) const
      {
        assert(false);
        abort();
      }

      //! jacobian of local function 
      template <class QuadratureType>
      void jacobian(const QuadratureType& quad,
                    const int quadPoint, 
                    JacobianRangeType& result) const 
      {
        assert(false);
        abort();
      }

      //! init local function
      void init(const EntityType& en) 
      {
        geometry_ = &(en.geometry());
        quadId_ = -1;
        qp_ = -1;
      } 

    private:
      const FunctionType& function_;
      const GeometryImp* geometry_;
      mutable DomainType global_;
      mutable RangeType result_;
      mutable int quadId_; 
      mutable int qp_;
    };

    public:
    //! type of local function to export 
    typedef LocalFunction LocalFunctionType; 

    // reference to function this local belongs to
    DiscreteFunctionAdapter(const std::string name, const FunctionType& f, const GridPartType& gridPart) 
      : BaseType(space_)
      , space_(gridPart)
      , function_(f)
      , name_(name)
    {}

    // reference to function this local belongs to
    DiscreteFunctionAdapter(const DiscreteFunctionAdapter& org) 
      : BaseType(org)
      , space_(org.space_)
      , function_(org.function_)
      , name_(org.name_)
    {}

    //! evaluate function on local coordinate local 
    void evaluate(const DomainType& global, RangeType& result) const 
    {
      function_.evaluate(global,result);  
    }

    //! return const local function object 
    const LocalFunctionType localFunction(const EntityType& en) const 
    {
      return LocalFunctionType(en,*this);
    }

    //! return local function object 
    LocalFunctionType localFunction(const EntityType& en) 
    {
      return LocalFunctionType(en,*this);
    }

    //! return name of function
    const std::string& name() const 
    {
      return name_;
    }

    bool write_xdr(std::string filename) const { return true; }
    bool write_ascii(std::string filename) const { return true; }
    bool write_pgm(std::string filename) const { return true; }

    bool read_xdr(std::string filename) const { return true; }
    bool read_ascii(std::string filename) const { return true; }
    bool read_pgm(std::string filename) const { return true; }

  private:    
    DiscreteFunctionSpaceType space_;
    //! reference to function 
    const FunctionType& function_; 
    
    const std::string name_;
  };



  /** @} end documentation group */

} // end namespace Dune

#include "discretefunction.cc"

#endif
