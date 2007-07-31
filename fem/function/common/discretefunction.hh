#ifndef DUNE_DISCRETEFUNCTION_HH
#define DUNE_DISCRETEFUNCTION_HH

//- system includes 
#include <string>

//- Dune inlcudes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/grid.hh>


//- local includes 
#include "function.hh"
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/objectstack.hh>
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

  /** Base class for determing whether a class is a discrete function or not. 
  */
  class IsDiscreteFunction
  {
  };
  
  /** Base class for determing whether a function has local functions or not.
  */
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

    //! type of discrete function space for discrete function 
    typedef typename DiscreteFunctionTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    //! type of domain field, i.e. type of coordinate component
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! type of range field, i.e. dof type 
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    //! type of domain, i.e. type of coordinates 
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! type of range, i.e. result of evaluation 
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! type of jacobian, i.e. type of evaluated gradient 
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
 
    //! Type of the underlying grid
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    //! Type of the discrete function implementation
    typedef typename DiscreteFunctionTraits::DiscreteFunctionType DiscreteFunctionType;

    //! Type of exported local function
    typedef typename DiscreteFunctionTraits::LocalFunctionType LocalFunctionType;

    //! Type of the local function implementation
    typedef typename DiscreteFunctionTraits::LocalFunctionImp LocalFunctionImp;

    //! Type of the dof iterator used in the discrete function implementation.
    typedef typename DiscreteFunctionTraits::DofIteratorType DofIteratorType;

    //! Type of the constantdof iterator used in the discrete function implementation
    typedef typename DiscreteFunctionTraits::ConstDofIteratorType ConstDofIteratorType;

  public:
    //- Public Methods

    /** \brief Constructor storing discrete function space 
        \param[in] f discrete function space 
    */
    DiscreteFunctionInterface (const DiscreteFunctionSpaceType& f) 
      : FunctionType ( f ) {}

    /** \brief returns name of discrete function 
        \return string holding name of discrete function 
    */
    std::string name() const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().name()); 
      return asImp().name();
    }

    /** \brief returns total number of degrees of freedom, i.e. size of discrete function space 
        \return total number of dofs 
    */
    int size() const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size()); 
      return asImp().size();
    }

    /** \brief returns dof iterator pointing to the first degree of freedom of this discrete function 
        \return dof iterator pointing to first dof 
    */
    DofIteratorType dbegin () 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().dbegin ());
      return asImp().dbegin ();
    }

    /** \brief returns dof iterator pointing behind the last degree of freedom of this discrete function 
        \return dof iterator pointing behind the last dof 
    */
    DofIteratorType dend () 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().dend ());
      return asImp().dend ();
    }

    /** \brief returns dof iterator pointing to the first degree of freedom of this discrete function 
        \return dof iterator pointing to first dof 
    */
    ConstDofIteratorType dbegin () const  
    {
      return asImp().dbegin ();
    }

    /** \brief returns dof iterator pointing behind the last degree of freedom of this discrete function 
        \return dof iterator pointing behind the last dof 
    */
    ConstDofIteratorType dend () const  
    {
      return asImp().dend ();
    }

    /** \brief return object of LocalFunction of the discrete function associated with given entity
        \param[in] entity Entity to focus view of discrete function 
        \return LocalFunction associated with entity
    */
    template <class EntityType>
    LocalFunctionType localFunction(const EntityType& entity) const {
      return asImp().localFunction(entity);
    }

    /** \brief evaluate Function f 
        \param[in] arg global coordinate
        \param[out] dest f(arg)
    */ 
    void evaluate (const DomainType & arg, RangeType & dest) const 
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(arg,dest));
    }

    /** \brief evaluate Function f
        \param[in] diffVariable derivation determizer 
        \param[in] arg global coordinate
        \param[out] dest f(arg)
    */ 
    template <int derivation>
    void evaluate  ( const FieldVector<deriType, derivation> &diffVariable, 
                     const DomainType& arg, RangeType & dest) const 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(diffVariable,arg,dest));
    }

  protected:
    /** \brief return pointer to new object of local function implementation 
        \return pointer to new object of local function implementation
    */
    LocalFunctionImp* newObject() const {
      return asImp().newObject();
    }

 protected:
    //! \brief Barton-Nackman trick 
    DiscreteFunctionType& asImp() 
    { 
      return static_cast<DiscreteFunctionType&>(*this); 
    }

    //! \brief Barton-Nackman trick 
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

    using DiscreteFunctionInterfaceType :: asImp;

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

    //! type of mapping base class for this discrete function 
    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

    //! Type of the discrete function (Barton-Nackman parameter)
    typedef typename DiscreteFunctionTraits::DiscreteFunctionType DiscreteFunctionType;
  
    //! Type of the local function
    typedef typename DiscreteFunctionTraits::LocalFunctionType LocalFunctionType;

    //! Type of the local function implementation 
    typedef typename DiscreteFunctionTraits::LocalFunctionImp LocalFunctionImp;
    //! Type of object to create from stack 
    typedef LocalFunctionImp ObjectType;
    //! Type of the dof iterator
    typedef typename DiscreteFunctionTraits::DofIteratorType DofIteratorType;
  
    //! Type of the const dof iterator
    typedef typename DiscreteFunctionTraits::ConstDofIteratorType ConstDofIteratorType;

    //! type of local function stack 
    typedef ObjectStack < DiscreteFunctionDefaultType > LocalFunctionStorageType;

    friend class ObjectStack < DiscreteFunctionDefaultType >;
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

    /** \brief print all degrees of freedom of this function to stream (for debugging purpose)
        \param[out] s std::ostream (e.g. std::cout)
    */
    void print(std::ostream & s) const;

    /** \brief check for NaNs
        \return if one  of the dofs is NaN <b>false</b> is returned, otherwise <b>true</b> 
    */
    inline bool dofsValid () const
    {
      const ConstDofIteratorType end = this->dend();
      for( ConstDofIteratorType it = this->dbegin(); it != end; ++it )
      {
        if( *it != *it )
          return false;
      }

      return true;
    }
    
    /** \brief set all degrees of freedom to zero
    */
    void clear();

    /** \brief Set all DoFs to a scalar value
        \param[in] s scalar value to assign for 
    */
    inline DiscreteFunctionType &assign ( const RangeFieldType s )
    {
      const DofIteratorType end = this->dend();
      for( DofIteratorType it = this->dbegin(); it != end; ++it )
        *it = s;
      return asImp();
    }

    /** \brief axpy operation
        \param[in] g discrete function that is added 
        \param[in] c scalar value to scale 
    */
    void addScaled(const DiscreteFunctionType& g, const RangeFieldType& c);

    /** \brief Evaluate a scalar product of the dofs of two DiscreteFunctions
        \param[in] g discrete function for evaluating scalar product with  
        \return returns the scalar product of the dofs 
    */
    RangeFieldType scalarProductDofs(const DiscreteFunctionType& g) const;

    /** \brief assign all degrees of freedom from given discrete function using the dof iterators 
        \param[in] g discrete function which is copied 
        \return reference to this (i.e. *this)
    */
    virtual DiscreteFunctionDefaultType& assign(const MappingType& g);

    /** \brief add all degrees of freedom from given discrete function using the dof iterators 
        \param[in] g discrete function which is added to this discrete function 
        \return reference to this (i.e. *this)
    */
    virtual DiscreteFunctionDefaultType& operator += (const MappingType& g);

    /** \brief substract all degrees of freedom from given discrete function using the dof iterators 
        \param[in] g discrete function which is substracted from this discrete function 
        \return reference to this (i.e. *this)
    */
    virtual DiscreteFunctionDefaultType& operator -= (const MappingType &g);
 
    /** \brief multiply all degrees of freedom with given scalar factor using the dof iterators 
        \param[in] scalar factor with which all dofs are scaled 
        \return reference to this (i.e. *this)
    */
    virtual DiscreteFunctionDefaultType& operator *=(const RangeFieldType &scalar);

    /** \brief devide all degrees of freedom with given scalar factor using the dof iterators 
        \param[in] scalar factor with which all dofs are devided  
        \return reference to this (i.e. *this)
    */
    virtual DiscreteFunctionDefaultType& operator /= (const RangeFieldType &scalar);

    /** \brief write discrete function to file with given filename using xdr encoding
        \param[in] filename name of file to which discrete function should be written using xdr 
        \return <b>true</b> if operation was successful 
    */
    bool write_xdr(std::string filename) const { return true; }

    /** \brief write discrete function to file with given filename using ascii encoding
        \param[in] filename name of file to which discrete function should be written using ascii 
        \return <b>true</b> if operation was successful 
    */
    bool write_ascii(std::string filename) const { return true; }

    /** \brief write discrete function to file with given filename using pgm encoding
        \param[in] filename name of file to which discrete function should be written using pgm 
        \return <b>true</b> if operation was successful 
    */
    bool write_pgm(std::string filename) const { return true; }

    /** \brief read discrete function from file with given filename using xdr decoding
        \param[in] filename name of file from which discrete function should be read using xdr 
        \return <b>true</b> if operation was successful 
    */
    bool read_xdr(std::string filename) const { return true; }
    /** \brief read discrete function from file with given filename using ascii decoding
        \param[in] filename name of file from which discrete function should be read using ascii 
        \return <b>true</b> if operation was successful 
    */
    bool read_ascii(std::string filename) const { return true; }
    /** \brief read discrete function from file with given filename using pgm decoding
        \param[in] filename name of file from which discrete function should be read using pgm 
        \return <b>true</b> if operation was successful 
    */
    bool read_pgm(std::string filename) const { return true; }

  protected: 
    //this methods are used by the LocalFunctionStorage class 

    /** \brief return reference for local function storage  
        \return reference to local function storage 
    */
    LocalFunctionStorageType& localFunctionStorage() const 
    { 
      return lfStorage_; 
    }

  private:    
    // the local function storage stack 
    mutable LocalFunctionStorageType lfStorage_;
  }; // end class DiscreteFunctionDefault 

} // end namespace Dune
#include "discretefunction.cc"
#include "discretefunctionadapter.hh"

#endif
