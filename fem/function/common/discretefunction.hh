#ifndef DUNE_DISCRETEFUNCTION_HH
#define DUNE_DISCRETEFUNCTION_HH

//- system includes 
#include <string>

//- Dune inlcudes 
#include <dune/grid/common/grid.hh>

//- local includes 
#include <dune/fem/version.hh>
#include <dune/fem/misc/debug.hh>
#include <dune/fem/storage/objectstack.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/localfunctionwrapper.hh>

namespace Dune
{

  /** @addtogroup DiscreteFunction 
      The DiscreteFunction is responsible for the dof storage. This can be
      done in various ways an is left to the user. The user has to derive his
      own implementation from the DiscreteFunctionDefault class. If some of
      the implementations in the default class are for ineffecient for the
      dof storage in the derived class these functions can be overloaded.

      \remarks 
      The interface for using a DiscreteFunction is defined by
      the class DiscreteFunctionInterface.
  
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

  template <class Traits>
  class DiscreteFunctionDefault;
  
  //----------------------------------------------------------------------
  //-
  //-  --DiscreteFunctionInterface
  //-
  //----------------------------------------------------------------------
  /** This is the interface of a discrete function which describes the
      features of a discrete function. 
      It contains a local function and a dof iterator which can 
      iterate over all dofs of one level. Via the method access the local
      dofs and basis functions can be accessed for a given entity.
      The DOF-Iterators are STL-like Iterators, i.e. they can be dereferenced
      giving the corresponding DOF.

      \interfaceclass
  */
  template< class TraitsImp >
  class DiscreteFunctionInterface
  : public Function< typename TraitsImp :: DiscreteFunctionSpaceType,
                     typename TraitsImp :: DiscreteFunctionType >,
    public IsDiscreteFunction, 
    public HasLocalFunction
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of the implementaton (Barton-Nackman)
    typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;

    //! type of associated discrete function space
    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

     //! type of the discrete function interface (this type)
    typedef DiscreteFunctionInterface< Traits > DiscreteFunctionInterfaceType;
 
  private:
    typedef DiscreteFunctionInterfaceType ThisType;
    typedef Function< DiscreteFunctionSpaceType, DiscreteFunctionType > BaseType;

  public:
    //! type of domain field, i.e. type of coordinate component
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! type of range field, i.e. dof type 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain, i.e. type of coordinates 
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range, i.e. result of evaluation 
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    //! type of jacobian, i.e. type of evaluated gradient 
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! Type of the underlying grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    // type of the local function storage
    typedef typename Traits :: LocalFunctionStorageType LocalFunctionStorageType;

    //! type of local functions
    typedef typename LocalFunctionStorageType :: LocalFunctionType LocalFunctionType;

    //! Type of the dof iterator used in the discrete function implementation.
    typedef typename Traits :: DofIteratorType DofIteratorType;

    //! Type of the constantdof iterator used in the discrete function implementation
    typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

    //! type of mapping base class for this discrete function 
    typedef Mapping<DomainFieldType, RangeFieldType,
                    DomainType, RangeType> MappingType;

    template< class Operation >
    struct CommDataHandle
    {
      typedef typename DiscreteFunctionSpaceType
        :: template CommDataHandle< DiscreteFunctionType, Operation > :: Type
        Type;
    };


  protected:
    using BaseType :: asImp;

  protected:
    /** \brief Constructor storing discrete function space
     *
     *  \param[in]  dfSpace  discrete function space 
     */
    inline explicit
    DiscreteFunctionInterface ( const DiscreteFunctionSpaceType &dfSpace )
    : BaseType( dfSpace )
    {}

  private:
    // prohibit copying
    DiscreteFunctionInterface ( const ThisType &other );
    // prohibit assignment
    ThisType &operator= ( const ThisType &other );
    
  public:
    /** \brief obtain the name of the discrete function 
     *
     *  \returns string holding name of discrete function
     */
    inline const std :: string &name () const
    {
      return asImp().name();
    }

    /** \brief obtain a local function for an entity (read-only)
     *
     *  \param[in]  entity  Entity to focus view of discrete function on
     *  \returns a local function associated with the entity
     */
    template< class EntityType >
    inline const LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return asImp().localFunction( entity );
    }
    
    /** \brief obtain a local function for an entity (read-write)
     *
     *  \param[in]  entity  Entity to focus view of discrete function
     *  \returns a local function associated with the entity
     */
    template< class EntityType >
    inline LocalFunctionType localFunction ( const EntityType &entity )
    {
      return asImp().localFunction( entity );
    }

    /** \brief set all degrees of freedom to zero
     */
    inline void clear()
    {
      asImp().clear();
    }

    /** \brief obtain total number of DoFs
     *
     *  The number of DoFs (degrees of freedom) can also be seen as the size
     *  of the discrete function, i.e., the size of the vector that forms this
     *  discrete funciton.
     *
     *  \returns total number of DoFs for this discrete function
     */
    inline int size() const
    {
      return asImp().size();
    }
    
    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
      return asImp().block( index );
    }
    
    inline DofBlockPtrType block ( unsigned int index )
    {
      return asImp().block( index );
    }
 
    /** \brief obtain an iterator pointing to the first DoF (read-only)
     *
     *  \returns a DoF iterator pointing to first DoF (degre of freedom)
     */
    inline ConstDofIteratorType dbegin () const
    {
      return asImp().dbegin ();
    }

    /** \brief obtain an iterator pointing behind the last DoF (read-only)
     *  
     *  \returns a DoF iterator pointing behind the last DoF (degree of freedom)
     */
    inline ConstDofIteratorType dend () const 
    {
      return asImp().dend ();
    }


    /** \brief obtain an iterator pointing to the first DoF (read-write)
     *
     *  \returns a DoF iterator pointing to first DoF (degre of freedom)
     */
    inline DofIteratorType dbegin () 
    {
      return asImp().dbegin ();
    }

    /** \brief obtain an iterator pointing behind the last DoF (read-write)
     *  
     *  \returns a DoF iterator pointing behind the last DoF (degree of freedom)
     */
    inline DofIteratorType dend () 
    {
      return asImp().dend ();
    }

    /** \brief allocate a pointer to a consecutive array storing the DoFs
     *
     *  To support external packages, it is often required to have the DoFs
     *  in a consecutive array. This function ensures this, making a copy if
     *  necessary.
     *
     *  \note The allocated pointer has to be freed by freeDofPointer.
     *
     *  \note Only one DoF pointer may be allocated at a time.
     *
     *  \returns a pointer to a consecutive copy of the DoF vector
     */
    inline RangeFieldType *allocDofPointer ()
    {
      return asImp().allocDofPointer();
    }

    /** \brief allocate a pointer to a consecutive array storing the DoFs
     *
     *  This method serves two purposes:
     *  - The user cannot know, if the DoF array returned by allocDofPointer
     *    has to be freed.
     *  - If the DoF array is just a copy, the DoFs shall be stored back into
     *    the discrete function.
     *
     *  \note The pointer must have been allocated by allocDofPointer.
     *
     *  \note Only one DoF pointer may be allocated at a time.
     *
     *  \param[in]  dofPointer  pointer to the dof array previously allocated
     *                          by allocDofPointer
     */
    inline void freeDofPointer( RangeFieldType *dofPointer )
    {
      asImp().freeDofPointer( dofPointer );
    }

    /** \brief axpy operation
     *
     *  Adds s * g to this discrete function.
     * 
     *  \param[in]  g  discrete function to add
     *  \param[in]  s  scalar value to scale g with
     */
    inline void addScaled( const DiscreteFunctionType &g,
                           const RangeFieldType &s )
    {
      asImp().addScaled( g, s );
    }

    /** \brief Scalar product between the DoFs of two discrete functions
     *
     *  \note This is a parallel scalar product, so do not sum over all
     *        processes after calling scalarProductDofs!
     *
     *  \note It is assumed that the discrete function has been communicated
     *        (i.e., every local DoF hat the value of the corresponding global
     *        DoF).
     *
     *  \param[in]  other  discrete function to evaluate the scalar product with
     *
     *  \returns the scalar product of the DoF-vectors
     */
    inline RangeFieldType
    scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const
    {
      return asImp().scalarProductDofs( other );
    }

    /** \brief print all DoFs to a stream (for debugging purposes)
     *
     *  \param[in]  out  stream to print to
     */
    inline void print( std :: ostream &out ) const
    {
      asImp().print( out );
    }

    /** \brief check for NaNs
     *  \returns if one of the DoFs is NaN \b false is returned, otherwise \b true
     */
    inline bool dofsValid () const
    {
      return asImp().dofsValid();
    }
    
    /** \brief assign the DoFs of another discrete function to this one
     * 
     *  \param[in]  g  discrete function which is copied
     */
    inline void assign( const DiscreteFunctionType &g )
    {
      asImp().assign( g );
    }

    /** \brief return reference to data handle object */
    template< class Operation >
    typename CommDataHandle< Operation > :: Type
    dataHandle( const Operation *operation )
    {
      return asImp().dataHandle( operation );
    }

    /** \brief add another discrete function to this one
     *
     *  \param[in]  g  discrete function to add 
     *
     *  \returns a reference to this discrete function (i.e. *this)
     */
    inline DiscreteFunctionType &operator+= ( const DiscreteFunctionType &g )
    {
      return asImp() += g;
    }

    /** \brief substract all degrees of freedom from given discrete function using the dof iterators 
        \param[in] g discrete function which is substracted from this discrete function 
        \return reference to this (i.e. *this)
    */
    template <class DFType>
    DiscreteFunctionType& operator -= (const DFType& g) 
    {
      return asImp().operator-=( g );
    }
 
    /** \brief multiply all DoFs by a scalar factor
     *
     *  \param[in]  scalar  factor to muliply all DoFs by
     *
     *  \returns a reference to this discrete function (i.e. *this)
     */
    inline DiscreteFunctionType &operator*= ( const RangeFieldType &scalar )
    {
      return asImp() *= scalar;
    }

    /** \brief devide all DoFs by a scalar factor
     *
     *  \param[in]  scalar  factor to divide all DoFs by
     *
     *  \returns a reference to this discrete function (i.e. *this)
     */
    inline DiscreteFunctionType &operator/= ( const RangeFieldType &scalar )
    {
      return asImp() /= scalar;
    }

    /** \brief read the discrete function from a stream
     *
     *  \param[in]  in  stream to read the discrete function from
     */
    template< class StreamTraits >
    inline void read ( InStreamInterface< StreamTraits > &in )
    {
      asImp().read( in );
    }

    /** \brief write the discrete function into a stream
     *
     *  \param[in]  out  stream to write the discrete function to
     */
    template< class StreamTraits >
    inline void write ( OutStreamInterface< StreamTraits > &out ) const
    {
      asImp().write( out );
    }
    
    /** \brief read discrete function from file using XDR decoding
     * 
     *  \param[in]  filename  name of file to read from
     *
     *  \return \b true if operation was successful
     */
    virtual bool read_xdr ( const std :: string filename ) = 0;
   
    /** \brief write discrete function to file using XDR encoding
     *
     *  \param[in]  filename  name of file to write to
     *
     *  \return \b true if operation was successful
     */
    virtual bool write_xdr ( const std :: string filename ) const = 0;
    
    /** \brief read discrete function from file using ASCII decoding
     *
     *  \param[in]  filename  name of file to read from 
     *
     *  \returns \b true if operation was successful
     */
    virtual bool read_ascii ( const std :: string filename ) = 0;

    /** \brief write discrete function to file using ASCII encoding
     *
     *  \param[in]  filename  name of the file to write to 
     *  
     *  \returns \b true if operation was successful
     */
    virtual bool write_ascii ( const std :: string filename ) const = 0;

    /** \brief write discrete function to file with given filename using pgm encoding
        \param[in] filename name of file to which discrete function should be written using pgm 
        \return \b true if operation was successful 
    */
    virtual bool DUNE_DEPRECATED write_pgm(const std::string filename) const = 0;

    /** \brief read discrete function from file with given filename using pgm decoding
        \param[in] filename name of file from which discrete function should be read using pgm 
        \return \b true if operation was successful 
    */
    virtual bool DUNE_DEPRECATED read_pgm(const std::string filename) const = 0;
    
    /** \brief Enable this discrete function for dof compression, 
         i.e. during grdi changes a dof compression 
         is done when the DofManagers compress is called. 
    */
    inline void enableDofCompression()
    {
      asImp().enableDofCompression();
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
  template< class DiscreteFunctionTraits >
  class DiscreteFunctionDefault
  : public DiscreteFunctionInterface< DiscreteFunctionTraits > 
  { 
  public:
    typedef DiscreteFunctionTraits Traits;

  private:
    typedef DiscreteFunctionDefault< Traits > ThisType;
    typedef DiscreteFunctionInterface< Traits > BaseType;

    //! type of the discrete function (Barton-Nackman parameter)
    typedef typename DiscreteFunctionTraits :: DiscreteFunctionType
      DiscreteFunctionType;

  private:
    typedef DiscreteFunctionInterface< DiscreteFunctionTraits >
      DiscreteFunctionInterfaceType;

    typedef DiscreteFunctionDefault< DiscreteFunctionTraits >
      DiscreteFunctionDefaultType;

    enum { myId_ = 0 };

    typedef ParallelScalarProduct< DiscreteFunctionInterfaceType >
      ScalarProductType;
  
  public:
    //! type of discrete function space
    typedef typename DiscreteFunctionTraits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
    //! type of domain vector
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vector
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  
    //! type of domain field (usually a float type)
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! type of range field (usually a float type)
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

    //! type of mapping base class for this discrete function 
    typedef Mapping< DomainFieldType, RangeFieldType, DomainType, RangeType>
      MappingType;

     //! type of the dof iterator
    typedef typename DiscreteFunctionTraits :: DofIteratorType DofIteratorType;
    //! type of the const dof iterator
    typedef typename DiscreteFunctionTraits :: ConstDofIteratorType
      ConstDofIteratorType;
 
    //! type of the local function factory
    typedef typename Traits :: LocalFunctionFactoryType LocalFunctionFactoryType;
    //! type of the local function storage
    typedef typename Traits :: LocalFunctionStorageType LocalFunctionStorageType;

    //! type of local functions
    typedef typename LocalFunctionStorageType :: LocalFunctionType LocalFunctionType;

    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

    template< class Operation >
    struct CommDataHandle
    : public BaseType :: template CommDataHandle< Operation >
    {};

  private:
    // the local function storage 
    mutable LocalFunctionStorageType lfStorage_;

    DebugLock dofPointerLock_;

  protected:
    std :: string name_;
    ScalarProductType scalarProduct_;

  protected:
    using BaseType :: asImp;

  protected:
    /** \brief Constructor storing discrete function space and local function
     *         factory
     *
     *  The discrete function space is passed to the interface class and the
     *  local function storage is initialized.
     * 
     *  \param[in]  name       name of the discrete function
     *  \param[in]  dfSpace    discrete function space 
     *  \param[in]  lfFactory  local function factory
     */
    inline DiscreteFunctionDefault ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &dfSpace,
                                     const LocalFunctionFactoryType &lfFactory );

  private:
    // prohibit copying
    inline DiscreteFunctionDefault ( const ThisType & );

  public:
    inline ~DiscreteFunctionDefault ();

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );
    
  public:
    // Default Implementations
    // -----------------------

    /** \copydoc Dune::DiscreteFunctionInterface::name() const */
    const std :: string &name () const;

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */
    template< class EntityType >
    inline const LocalFunctionType localFunction ( const EntityType &entity ) const;
    
    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
    template< class EntityType >
    inline LocalFunctionType localFunction ( const EntityType &entity );
  
    /** \copydoc Dune::DiscreteFunctionInterface::clear() */
    inline void clear();
    
    /** \copydoc Dune::DiscreteFunctionInterface::allocDofPointer
     *  
     *  \note The default implementation make a copy of the DoF vector using
     *        the DoF iterators.
     */
    inline RangeFieldType *allocDofPointer ();

    /** \copydoc Dune::DiscreteFunctionInterface::freeDofPointer
     *
     *  \note The default implementation make a copy of the DoF vector using
     *        the DoF iterators.
     */
    inline void freeDofPointer( RangeFieldType *dofPointer );

    /** \copydoc Dune::DiscreteFunctionInterface::addScaled */
    void addScaled ( const DiscreteFunctionType &g, const RangeFieldType &s );
    
    /** \copydoc Dune::DiscreteFunctionInterface::scalarProductDofs */
    inline RangeFieldType
    scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const;

    /** \copydoc Dune::DiscreteFunctionInterface::print */
    void print ( std :: ostream &out ) const;
    
    /** \copydoc Dune::DiscreteFunctionInterface::dofsValid */
    inline bool dofsValid () const;
 
    /** \copydoc Dune::DiscreteFunctionInterface::assign(const DiscreteFunctionType &g) */
    inline void assign ( const DiscreteFunctionType &g );
    
    /** \copydoc Dune::DiscreteFunctionInterface::dataHandle */
    template< class Operation >
    typename CommDataHandle< Operation > :: Type
    dataHandle ( const Operation *operation );
 
    /** \copydoc Dune::Function::evaluate(const DomainType &x,RangeType &ret) const
     *
     *  The default implementation just does
     *  \code
     *  FieldVector< deriType, 0 > diffVariable;
     *  asImp().evaluate( diffVariable, x, ret );
     *  \endcode
     */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::Function::evaluate(const FieldVector<deriType,diffOrder> &diffVariable,const DomainType &x,RangeType &ret) const */
    template< int diffOrder >
    inline void evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                           const DomainType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::DiscreteFunctionInterface::operator+=(const DiscreteFunctionType &g) */
    inline DiscreteFunctionType &operator+= ( const DiscreteFunctionType &g );

    /** \brief substract all degrees of freedom from given discrete function using the dof iterators 
        \param[in] g discrete function which is substracted from this discrete function 
        \return reference to this (i.e. *this)
    */
    template <class DFType>
    DiscreteFunctionType& operator -= (const DFType& g);
 
    /** \brief multiply all DoFs with a scalar factor
     *
     *  \param[in]  scalar  factor to multiply DoFs with
     *  
     *  \returns reference to this discrete function (i.e. *this)
     */
    inline DiscreteFunctionType &operator*= ( const RangeFieldType &scalar );

    /** \brief devide all DoFs by a scalar factor
     *
     *  \param[in]  scalar  factor with which all dofs are devided
     *  
     *  \returns reference to this discrete function (i.e. *this)
     */
    inline DiscreteFunctionType &operator/= ( const RangeFieldType &scalar );

    /** \copydoc Dune::DiscreteFunctionInterface::read */
    template< class StreamTraits >
    inline void read ( InStreamInterface< StreamTraits > &in );

    /** \copydoc Dune::DiscreteFunctionInterface::write */
    template< class StreamTraits >
    inline void write ( OutStreamInterface< StreamTraits > &out ) const;
    
    /** \copydoc Dune::DiscreteFunctionInterface::read_xdr
     * 
     *  \note The default implementation uses the read method to read the
     *        discrete function from an XDRFileInStream.
     */
    virtual bool read_xdr ( const std :: string filename );
 
    /** \copydoc Dune::DiscreteFunctionInterface::write_xdr
     *
     *  \note The default implementation uses the write method to write the
     *        discrete function to an XDRFileOutStream.
     */
    virtual bool write_xdr ( const std :: string filename ) const;
    
    /** \copydoc Dune::DiscreteFunctionInterface::read_ascii
     * 
     *  \note The default implementation uses the read method to read the
     *        discrete function from an ASCIIInStream.
     */
    virtual bool read_ascii ( const std :: string filename );
 
    /** \copydoc Dune::DiscreteFunctionInterface::write_ascii
     *
     *  \note The default implementation uses the write method to write the
     *        discrete function to an ASSCIIOutStream.
     */
    virtual bool write_ascii ( const std :: string filename ) const;

    /** \copydoc DiscreteFunctionInterface::read_pgm */
    virtual bool DUNE_DEPRECATED read_pgm ( const std :: string filename) const;

    /** \copydoc DiscreteFunctionInterface::write_pgm */
    virtual bool DUNE_DEPRECATED write_pgm ( const std :: string filename) const;

    /** \copydoc Dune::DiscreteFunctionInterface::enableDofCompression()
     *
     *  \note The default implementation does nothing.
     */
    inline void enableDofCompression ();


  public:
    // Non-Interface Methods
    // ---------------------
  
    inline bool operator== ( const DiscreteFunctionType &g ) const;
    
    inline bool operator!= ( const DiscreteFunctionType &g ) const;
  
    /** \brief obtain the local function storage
     *
     *  \returns a reference to the local function storage
     */
    inline LocalFunctionStorageType &localFunctionStorage () const;


  private:
    // Unimplemented Interface Methods
    // -------------------------------
    
    int size () const;

    ConstDofBlockPtrType block ( unsigned int index ) const;
    DofBlockPtrType block ( unsigned int index );

    ConstDofIteratorType dbegin () const;
    ConstDofIteratorType dend () const;
    
    DofIteratorType dbegin ();
    DofIteratorType dend ();
  }; // end class DiscreteFunctionDefault 


  template< class DiscreteFunction >
  class ManagedDiscreteFunction;
  
///@} 
} // end namespace Dune

#include "discretefunction_inline.hh"

#include "discretefunctionadapter.hh"
#endif
