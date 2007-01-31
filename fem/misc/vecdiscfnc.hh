#ifndef ADI_VECDISCFUNC_HH
#define ADI_VECDISCFUNC_HH

// Dune Includes
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/localfunction.hh>
#include <dune/common/array.hh>

using namespace Dune;

namespace Adi {

  //- Forward declarations
  template <class CompositeFunctionType>
  class LocalFunctionComposite;
  template <class CompositeFunctionType>
  class DofIteratorComposite;

  //- class CompositeDiscreteFunction

  //! \brief Composite of the composite pattern for discrete function
  template <class CompositeFunctionSpaceImp, template <class> class ContainedFunctionImp>
  class CompositeDiscreteFunction :
    public DiscreteFunctionDefault<
    CompositeFunctionSpaceImp,
    DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                   ContainedFunctionImp> >,
    LocalFunctionComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                     ContainedFunctionImp> >,
    CompositeDiscreteFunction<CompositeFunctionSpaceImp, ContainedFunctionImp> >
  {
    //- Local typedefs
    typedef DiscreteFunctionDefault<
      CompositeFunctionSpaceImp,
      DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                     ContainedFunctionImp> >,
      LocalFunctionComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                       ContainedFunctionImp> >,
      CompositeDiscreteFunction<CompositeFunctionSpaceImp, 
                                ContainedFunctionImp>
    > DiscreteFunctionDefaultType;

    typedef CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                      ContainedFunctionImp> MyType;
    
    enum { DimRange = CompositeFunctionSpaceImp::DimRange };
    
    //- Friends
    friend class DofIteratorComposite<MyType>;
    friend class LocalFunctionComposite<MyType>;

  public:
    //- Public typedefs
    // Range definition
    typedef typename CompositeFunctionSpaceImp::RangeField RangeFieldType;

    // Function space definition
    typedef CompositeFunctionSpaceImp CompositeFunctionSpaceType;
    typedef typename CompositeFunctionSpaceType::ContainedFunctionSpaceType ContainedFunctionSpaceType;
    // * Export the composite function space as official function space
    typedef CompositeFunctionSpaceType FunctionSpaceType;
    
    // Function definition (global and local)
    typedef ContainedFunctionImp<ContainedFunctionSpaceType> ContainedFunctionType;
    typedef CompositeDiscreteFunction<CompositeFunctionSpaceType, 
                                               ContainedFunctionImp> CompositeFunctionType;
    typedef typename ContainedFunctionType::LocalFunctionType ContainedLocalFunctionType;
    typedef typename CompositeFunctionType::LocalFunctionType CompositeLocalFunctionType;
    // * Export the composite local and global functions
    typedef CompositeLocalFunctionType LocalFunctionType;
    typedef MyType DiscreteFunctionType;
    
    // Dof iterator types
    typedef DofIteratorComposite<MyType> DofIteratorType;
    typedef ConstDofIteratorDefault<DofIteratorType> ConstDofIteratorType;

    //- Constructors & destructors
    //! Constructor
    CompositeDiscreteFunction(const CompositeFunctionSpaceType& f);

    //! Destructor
    ~CompositeDiscreteFunction();

    //! \brief Copy constructor
    //! The copy constructor implements a deep copy
    CompositeDiscreteFunction(const CompositeFunctionType& other);

    //- Methods
    //virtual MyType& operator=(const MyType& other);

    //! Get a reference to a contained function
    ContainedFunctionType& getDiscreteFunction(int position);

    //! return object of type LocalFunctionType
    LocalFunctionComposite<CompositeFunctionType> newLocalFunction();

    //! return reference to this 
    //! this methods is only to fullfill the interface as parameter classes 
    CompositeFunctionType & argument    () { return *this; }

    //! return reference to this 
    //! this methods is only to fullfill the interface as parameter classes 
    const CompositeFunctionType & argument () const { return *this; }
    
    //! return reference to this 
    //! this methods is only to fullfill the interface as parameter classes 
    CompositeFunctionType & destination () { return *this; }

    //! update LocalFunction to given Entity en  
    template <class EntityType> 
    void localFunction(EntityType& en, 
                       LocalFunctionComposite<CompositeFunctionType>& lf);

        //! we use the default implementation 
    DofIteratorType dbegin ();
    
    //! points behind the last dof of type cc
    DofIteratorType dend   ();

    //! the const versions 
    //! we use the default implementation 
    ConstDofIteratorType dbegin () const;
    
    //! points behind the last dof of type cc
    ConstDofIteratorType dend   () const;

    //! set all dofs to zero  
    void clear( );

    //! set all dof to value val 
    void set( RangeFieldType val );

  private:
    // * Helper function
    ContainedLocalFunctionType& getChild(int i, LocalFunctionType& lf) {
      return *(lf.children_[i]);
    }
    

    //- Private typedefs
    typedef FieldVector<ContainedFunctionType*,
                        DimRange> ContainerType;
    
    //- Data members
    ContainerType children_;

  }; // end class CompositeDiscreteFunction

  //- class LocalFunctionComposite
  //! \brief Local function for CompositeDiscreteFunction
  template <class CompositeFunctionImp>
  class LocalFunctionComposite :
    public LocalFunctionDefault<typename CompositeFunctionImp::CompositeFunctionSpaceType,
                                LocalFunctionComposite<CompositeFunctionImp> >
  {
  public:
    //- General typedefs
    // Function spaces
    typedef typename CompositeFunctionImp::CompositeFunctionSpaceType CompositeFunctionSpaceType;
    typedef typename CompositeFunctionImp::ContainedFunctionSpaceType ContainedFunctionSpaceType;
   typedef typename CompositeFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
 
   
    // Discrete functions (global and local)
    typedef CompositeFunctionImp CompositeFunctionType;
    typedef typename CompositeFunctionImp::ContainedFunctionType ContainedFunctionType;
    typedef typename CompositeFunctionImp::LocalFunctionType CompositeLocalFunctionType;
    typedef typename ContainedFunctionType::LocalFunctionType ContainedLocalFunctionType;
    typedef CompositeLocalFunctionType LocalFunctionType;

    // This type
    typedef LocalFunctionComposite<CompositeFunctionImp> ThisType;

    // (Global) range types
    typedef typename FunctionSpaceType::RangeField RangeFieldType;
    typedef typename FunctionSpaceType::Range RangeType;
    typedef typename FunctionSpaceType::Domain DomainType;

  private:
    //- Local typedefs
 
    // * The range vector of the scalar function space
    typedef typename ContainedFunctionSpaceType::Range ContainedRangeType;
    
    enum { DimRange = CompositeFunctionSpaceType::DimRange };

    //- Friends
    //friend class CompositeFunctionImp;
    friend 
    ContainedLocalFunctionType& 
    CompositeFunctionImp::getChild(int i, ThisType& lf);

  public:
    //- Constructors & destructors
    //! Constructor
    LocalFunctionComposite(const CompositeFunctionSpaceType& spc,
                           CompositeFunctionType& f);

    //! Copy constructor
    LocalFunctionComposite(const ThisType& other);

    //! Destructor
    ~LocalFunctionComposite();

    //- Methods
    void setNDof() { nDofChildren_ = children_[0]->numberOfDofs(); }

    //! access to dof number num, all dofs of the dof entity
    RangeFieldType& operator [] (int num);
  
    //! access to dof number num, all dofs of the dof entity
    const RangeFieldType& operator [] (int num) const;

    //! return number of degrees of freedom 
    int numberOfDofs () const;

    //! sum over all local base functions 
    template <class EntityType> 
    void evaluate (EntityType& en,
                   const DomainType& x, 
                   RangeType& ret) const;
  
    //! sum over all local base functions evaluated on given quadrature point
    template <class EntityType, class QuadratureType> 
    void evaluate (EntityType& en,
                   QuadratureType &quad,
                   int quadPoint,
                   RangeType& ret) const;

    //! update local function for given Entity  
    template <class EntityType> 
    bool init (CompositeFunctionImp& df, EntityType &en ) const;

  private:
    //- Local typedefs
    typedef FieldVector<ContainedLocalFunctionType*,
                        DimRange> ContainerType;

    //- Local methods

    
    //- Data members
    //! the corresponding function space which provides the base function set
    const CompositeFunctionSpaceType& fSpace_;
    
    //! The collection of local functions
    ContainerType children_;

    //! The dof number of the contained function (cached)
    mutable int nDofChildren_;

  }; // end class LocalFunctionComposite

  //- class DofIteratorComposite
  //! \brief DofIterator for CompositeDiscreteFunction
  template <class CompositeFunctionImp>
  class DofIteratorComposite :
    public DofIteratorDefault<typename CompositeFunctionImp::RangeFieldType,
                              DofIteratorComposite<CompositeFunctionImp> >
  {
    //- Local typedefs
    typedef typename CompositeFunctionImp::ContainedFunctionType ContainedFunctionType;
    typedef typename ContainedFunctionType::DofIteratorType ContainedDofIteratorType;
    typedef DofIteratorComposite<CompositeFunctionImp> ThisType;
    typedef typename CompositeFunctionImp::CompositeFunctionSpaceType CompositeFunctionSpaceType;

  public:
    //- Public typedefs
    typedef typename CompositeFunctionImp::RangeFieldType DofType;
    typedef CompositeFunctionImp CompositeFunctionType;

    //- Constructors & destructors
    //! Default constructor
    DofIteratorComposite() :
      children_(0) {}

    //! Constructor
    //! Produces an end iterator
    DofIteratorComposite(CompositeFunctionType& df, int dummy) :
      children_(&df.children_),
      funcIter_(df.children_.begin()),
      endFuncIter_(df.children_.end()),
      dofIter_((*funcIter_)->dbegin()),
      endDofIter_((*funcIter_)->dend()) {
      while(funcIter_ != endFuncIter_) {
        dofIter_ = (*funcIter_)->dend();
        ++funcIter_;
      }
    }

    //! Constructor
    //! Generates a begin iterator
    DofIteratorComposite(CompositeFunctionType& df) : 
      children_(&df.children_),
      funcIter_(df.children_.begin()),
      endFuncIter_(df.children_.end()),
      dofIter_((*funcIter_)->dbegin()),
      endDofIter_((*funcIter_)->dend()) {}

    //! Copy constructor
    DofIteratorComposite(const ThisType& other);

    //- Methods
    //! Assignment operator
    ThisType& operator= (const ThisType& other);
    
    //! return dof
    DofType & operator* ();

    //! return dof read only 
    const DofType & operator* () const;

    //! go next dof
    DofIteratorComposite<CompositeFunctionType> & operator++ ();
  
    //! compare
    bool operator == (const ThisType& I ) const;

    //! compare 
    bool operator != (const ThisType& I ) const; 

    //! set dof iterator back to begin , for const and not const Iterators
    void reset ();
  
   private:
    //- Local typedefs
    typedef typename CompositeFunctionType::ContainerType ContainerType;

    //- Data members
    //! Reference to the container
    ContainerType* children_;

    //! Contained discrete functions
    typename ContainerType::Iterator funcIter_;

    //! End iterator of the discrete function container
    typename ContainerType::Iterator endFuncIter_;

    //! The actual dof iterator
    ContainedDofIteratorType dofIter_;

    //! The end dof iterator of the actual function
    ContainedDofIteratorType endDofIter_;
  }; // end class DofIteratorComposite

} // end namespace Adi

#include "vecdiscfnc.cc"

#endif
