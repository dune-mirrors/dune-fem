#ifndef DUNE_FEM_DISCRETEFUNCTION_HH
#define DUNE_FEM_DISCRETEFUNCTION_HH

#include <cassert>

#include <complex>
#include <memory>
#include <ostream>
#include <string>
#include <typeindex>

#include <dune/common/dynvector.hh>

#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/storage/referencevector.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    /** @addtogroup DiscreteFunction
     *  The DiscreteFunction is responsible for the dof storage. This can be
     *  done in various ways an is left to the user. The user has to derive his
     *  own implementation from the DiscreteFunctionDefault class. The implementations
     *  in the default class which are ineffecient for the dof storage in the derived
     *  class can be overloaded.
     *
     *  \remarks
     *  The interface for using a DiscreteFunction is defined by
     *  the class DiscreteFunctionInterface.
     *
     *  @{
     */

    /** \brief base class for determing whether a class is a discrete function or not */
    class IsDiscreteFunction
    {};

    /** \brief base class for determing whether a function has local functions or not */
    class HasLocalFunction
    {};


    template< class DiscreteFunction >
    struct DiscreteFunctionTraits;

    template< class Traits >
    class DiscreteFunctionDefault;

    //----------------------------------------------------------------------
    //-
    //-  --DiscreteFunctionInterface
    //-
    //----------------------------------------------------------------------
    /** This is the interface of a discrete function which describes the
     *  features of a discrete function.
     *  It contains a local function and a dof iterator which can
     *  iterate over all dofs of one level. Via the method access the local
     *  dofs and basis functions can be accessed for a given entity.
     *  The DOF-Iterators are STL-like Iterators, i.e. they can be dereferenced
     *  giving the corresponding DOF.
     *
     *  \interfaceclass
     */
    template< class Impl >
    class DiscreteFunctionInterface
    : public Fem::Function< typename DiscreteFunctionTraits< Impl >::DiscreteFunctionSpaceType::FunctionSpaceType, Impl >,
      public IsDiscreteFunction,
      public HasLocalFunction
    {
      typedef DiscreteFunctionInterface< Impl > ThisType;
      typedef Fem::Function< typename DiscreteFunctionTraits< Impl >::DiscreteFunctionSpaceType::FunctionSpaceType, Impl > BaseType;

    public:
      //! type of the traits
      typedef DiscreteFunctionTraits< Impl > Traits;

      //! type of the implementaton (Barton-Nackman)
      typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;

      //! type of associated discrete function space
      typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! type of associated function space
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

      //! type of the discrete function interface (this type)
      typedef DiscreteFunctionInterface< Impl > DiscreteFunctionInterfaceType;

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

      //! type of the underlying grid part
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename GridPartType::GridViewType GridView;

      //! type of the underlying grid
      typedef typename DiscreteFunctionSpaceType :: GridType GridType;

      //! type of local functions
      typedef typename Traits :: LocalFunctionType LocalFunctionType;

      //! type of the dof vector used in the discrete function implementation
      typedef typename Traits :: DofVectorType  DofVectorType;

      //! type of the dof iterator used in the discrete function implementation
      typedef typename Traits :: DofIteratorType DofIteratorType;

      //! type of the constantdof iterator used in the discrete function implementation
      typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

      typedef typename Traits :: DofType DofType;
      typedef typename Traits :: DofBlockType DofBlockType;
      typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
      typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
      typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

      //! type of mapping base class for this discrete function
      typedef typename BaseType :: MappingType MappingType;

      typedef typename DiscreteFunctionSpaceType::LocalBlockIndices BlockIndices;

      //! size of the dof blocks
      static constexpr std::size_t blockSize = Hybrid::size( BlockIndices() );

      template< class Operation >
      struct CommDataHandle
      {
        typedef typename DiscreteFunctionSpaceType
          :: template CommDataHandle< DiscreteFunctionType, Operation > :: Type
          Type;
      };

      //! type of entity local functions are defined on
      typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

    protected:
      using BaseType::asImp;

      //! default constructor
      DiscreteFunctionInterface () = default;

      DiscreteFunctionInterface ( const ThisType& ) = default;
      DiscreteFunctionInterface ( ThisType && ) = default;
    public:
      ThisType& operator= ( ThisType&& ) = delete;
      ThisType &operator= ( const ThisType& ) = delete;

      DofVectorType &dofVector()
      {
        return asImp().dofVector();
      }
      const DofVectorType &dofVector() const
      {
        return asImp().dofVector();
      }

      /** \brief obtain the name of the discrete function
       *
       *  \returns string holding name of discrete function
       */
      const std::string &name () const
      {
        return asImp().name();
      }

      /** \brief obtain the name of the discrete function
       *
       *  \returns string holding name of discrete function
       */
      std::string &name ()
      {
        return asImp().name();
      }

      /** \brief obtain an upper bound on the polynomial order of the underlying space.
       */
      const std::string &order () const
      {
        return asImp().order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous() const
      {
        return asImp().continuous();
      }

      /** \brief obtain a reference to the corresponding DiscreteFunctionSpace */
      const DiscreteFunctionSpaceType &space () const
      {
        return asImp().space();
      }

      /** \brief obtain a reference to the underlying grid part */
      const GridPartType &gridPart () const
      {
        return asImp().gridPart();
      }

      /** \brief obtain a local function for an entity (read-write)
       *
       *  \param[in]  entity  Entity to focus view of discrete function
       *  \returns a local function associated with the entity
       */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      LocalFunctionType localFunction ( const EntityType &entity )
      {
        return asImp().localFunction( entity );
      }

      /** \brief obtain a local function for an entity (read-write)
       *
       *  \param[in]  entity  Entity to focus view of discrete function
       *  \returns a local function associated with the entity
       */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      const LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        return asImp().localFunction( entity );
      }

      /** \brief obtain an uninitialized local function (read-write)
       *
       *  \note before calling any method of the local function initialize it passing an entity
       *
       *  \returns an uninitialized local function
       */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      LocalFunctionType localFunction ()
      {
        return asImp().localFunction();
      }


      /** \brief add scaled local Dofs to dof vector associated with the entity
       *
       *  \param[in]  entity    Entity to focus view of discrete function
       *  \param[in]  s         scaling factor
       *  \param[in]  localDofs the local dofs vector to be added
       */
      template< class LocalDofs >
      void addScaledLocalDofs ( const EntityType &entity, const RangeFieldType &s, const LocalDofs &localDofs )
      {
        asImp().addScaledLocalDofs( entity, s, localDofs );
      }

      /** \brief add local Dofs to dof vector associated with the entity
       *
       *  \param[in]  entity    Entity to focus view of discrete function
       *  \param[in]  localDofs the local dofs vector to be added
       */
      template< class LocalDofs >
      void addLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        asImp().addLocalDofs( entity, localDofs );
      }

      /** \brief set local Dofs to dof vector associated with the entity
       *
       *  \param[in]  entity    Entity to focus view of discrete function
       *  \param[in]  localDofs the local dofs vector to be set
       */
      template< class LocalDofs >
      void setLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        asImp().setLocalDofs( entity, localDofs );
      }

      /** \brief fill local Dofs to dof vector associated with the entity
       *
       *  \param[in]   entity    Entity to focus view of discrete function
       *  \param[out]  localDofs the local dofs vector to be set
       *
       *  \note localDofs should have sufficient size to store the dof values
       */
      template< class Vector >
      void getLocalDofs ( const EntityType &entity, Vector &localDofs ) const
      {
        asImp().getLocalDofs( entity, localDofs );
      }

      /** \brief obtain an uninitialized local function (read-write)
       *
       * \note before calling any method of the local function initialize it passing an entity
       *
       *  \returns an uninitialized local function
       */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      const LocalFunctionType localFunction () const
      {
        return asImp().localFunction();
      }

      /** \brief set all degrees of freedom to zero */
      void clear()
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
      int size() const
      {
        return asImp().size();
      }

      /** \brief obtain total number of blocks, i.e. size / blockSize.
       *
       *  The number of blocks of DoFs (degrees of freedom) can also be seen
       *  as the size of the discrete function divided by the blockSize.
       *
       *  \returns total number of DoFs blocks
       */
      int blocks() const
      {
        return asImp().blocks();
      }

      /** \brief obtain an iterator pointing to the first DoF (read-only)
       *
       *  \returns a DoF iterator pointing to first DoF (degre of freedom)
       */
      ConstDofIteratorType dbegin () const
      {
        return asImp().dbegin ();
      }

      /** \brief obtain an iterator pointing behind the last DoF (read-only)
       *
       *  \returns a DoF iterator pointing behind the last DoF (degree of freedom)
       */
      ConstDofIteratorType dend () const
      {
        return asImp().dend ();
      }

      /** \brief obtain an iterator pointing to the first DoF (read-write)
       *
       *  \returns a DoF iterator pointing to first DoF (degre of freedom)
       */
      DofIteratorType dbegin ()
      {
        return asImp().dbegin ();
      }

      /** \brief obtain an iterator pointing behind the last DoF (read-write)
       *
       *  \returns a DoF iterator pointing behind the last DoF (degree of freedom)
       */
      DofIteratorType dend ()
      {
        return asImp().dend ();
      }

      /** \brief axpy operation
       *
       *  Adds s * g to this discrete function.
       *
       *  \param[in]  s  scalar value to scale g with
       *  \param[in]  g  discrete function to add
       */
      void axpy( const RangeFieldType &s, const DiscreteFunctionInterfaceType &g )
      {
        asImp().axpy( s, g );
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
      template <class DFType>
      RangeFieldType scalarProductDofs ( const DiscreteFunctionInterface< DFType > &other ) const
      {
        return asImp().scalarProductDofs( other.asImp() );
      }

      /** \brief Squared small l^2 norm of all dofs
       *
       *  \note This is already parallel, so do not sum over all
       *        processes after calling scalarProductDofs!
       *
       *  \note It is assumed that the discrete function has been communicated
       *        (i.e., every local DoF hat the value of the corresponding global
       *        DoF).
       *
       *  \returns the squared norm of the DoF-vectors
       */
      typename Dune::FieldTraits< RangeFieldType >::real_type normSquaredDofs ( ) const
      {
        return asImp().normSquaredDofs( );
      }

      /** \brief print all DoFs to a stream (for debugging purposes)
       *
       *  \param[in]  out  stream to print to
       */
      void print( std :: ostream &out ) const
      {
        asImp().print( out );
      }

      /** \brief check for NaNs
       *  \returns if one of the DoFs is NaN \b false is returned, otherwise \b true
       */
      bool dofsValid () const
      {
        return asImp().dofsValid();
      }

      /** \brief assign the DoFs of another discrete function to this one
       *
       *  \param[in]  g  discrete function which is copied
       */
      template < class DFType >
      void assign( const DiscreteFunctionInterface< DFType > &g )
      {
        asImp().assign( g );
      }

      /** \brief return reference to data handle object */
      template< class Operation >
      typename CommDataHandle< Operation >::Type dataHandle( const Operation &operation )
      {
        return asImp().dataHandle( operation );
      }

      /** \brief do default communication of space for this discrete function */
      void communicate()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().communicate() );
      }

      template <class DFType>
      inline bool compare ( const DiscreteFunctionInterface<DFType> &g ) const
      {
        return asImp().compare( g );
      }

      /** \brief add another discrete function to this one
       *
       *  \param[in]  g  discrete function to add
       *
       *  \returns a reference to this discrete function (i.e. *this)
       */
      template < class DFType >
      DiscreteFunctionType &operator+=(const DiscreteFunctionInterface< DFType > &g)
      {
        return asImp().operator+=( g );
      }

      /** \brief substract all degrees of freedom from given discrete function using the dof iterators
       *
       *  \param[in] g discrete function which is substracted from this discrete function
       *
       *  \returns reference to this (i.e. *this)
       */
      template < class DFType >
      DiscreteFunctionType &operator-=(const DiscreteFunctionInterface< DFType > &g)
      {
        return asImp().operator-=( g );
      }

      /** \brief multiply all DoFs by a scalar factor
       *
       *  \param[in]  scalar  factor to muliply all DoFs by
       *
       *  \returns a reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator*= ( const RangeFieldType &scalar )
      {
        return asImp() *= scalar;
      }

      /** \brief devide all DoFs by a scalar factor
       *
       *  \param[in]  scalar  factor to divide all DoFs by
       *
       *  \returns a reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator/= ( const RangeFieldType &scalar )
      {
        return asImp() /= scalar;
      }

      /** \brief read the discrete function from a stream
       *
       *  \param[in]  in  stream to read the discrete function from
       *
       *  \note This call will automatically enable dof compression for this
       *        discrete function.
       */
      template< class StreamTraits >
      void read ( InStreamInterface< StreamTraits > &in )
      {
        asImp().read( in );
      }

      /** \brief write the discrete function into a stream
       *
       *  \param[in]  out  stream to write the discrete function to
       */
      template< class StreamTraits >
      void write ( OutStreamInterface< StreamTraits > &out ) const
      {
        asImp().write( out );
      }

      /** \brief Enable this discrete function for dof compression,
       *   i.e. during grid changes a dof compression
       *   is done when the DofManagers compress is called.
       */
      void enableDofCompression()
      {
        asImp().enableDofCompression();
      }

      //TODO: this needs to be revised, the definition should be in GridPart
      typedef LoadBalanceLeafData< ThisType > DefaultLoadBalanceContainsCheckType;
      DefaultLoadBalanceContainsCheckType defaultLoadBalanceContainsCheck() const
      {
        return DefaultLoadBalanceContainsCheckType( *this );
      }
    };



    //*************************************************************************
    //
    //  --DiscreteFunctionDefault
    //
    //! Default implementation of the discrete function. This class is
    //! responsible for the dof storage. Different implementations of the
    //! discrete function use different dof storage.
    //! The default implementation provides +=, -= and so on operators and
    //! a DofIterator access, which can run over all dofs in an efficient way.
    //! Furthermore with an entity you can access a local function to evaluate
    //! the discrete function by multiplying the dofs and the basefunctions.
    //!
    //*************************************************************************
    template< class Impl >
    class DiscreteFunctionDefault
    : public DiscreteFunctionInterface< Impl > ,
      public PersistentObject
    {
      typedef DiscreteFunctionDefault< Impl > ThisType;
      typedef DiscreteFunctionInterface< Impl > BaseType;

    public:
      typedef typename BaseType :: Traits Traits;

      //! type of the discrete function (Barton-Nackman parameter)
      typedef Impl DiscreteFunctionType;

      typedef typename BaseType::DiscreteFunctionInterfaceType DiscreteFunctionInterfaceType;

    private:
      typedef DiscreteFunctionDefault< Impl > DiscreteFunctionDefaultType;

      enum { myId_ = 0 };

    protected:
      typedef ParallelScalarProduct< DiscreteFunctionInterfaceType > ScalarProductType;

    public:
      //! type of discrete function space
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! type of the underlying grid part
      typedef typename BaseType::GridPartType GridPartType;

      //! type of domain vector
      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      //! type of range vector
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
      //! type of jacobian
      typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
      //! type of hessian
      typedef typename DiscreteFunctionSpaceType :: HessianRangeType HessianRangeType;

      //! type of domain field (usually a float type)
      typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
      //! type of range field (usually a float type)
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

      //! type of the dof iterator
      typedef typename Traits :: DofIteratorType DofIteratorType;
      //! type of the const dof iterator
      typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

      //! type of DofVector
      typedef typename Traits :: DofVectorType     DofVectorType;

      //! type of LocalDofVector
      typedef typename Traits :: LocalDofVectorType LocalDofVectorType;
      //! type of LocalDofVector
      typedef typename Traits :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      //! type of local functions
      typedef typename BaseType :: LocalFunctionType LocalFunctionType;
      typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;

      typedef typename BaseType :: DofBlockType DofBlockType;
      typedef typename BaseType :: ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType :: DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType :: ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType :: EntityType EntityType ;

      typedef typename BaseType :: DofType DofType;

      //! size type of the block vector
      typedef typename DofVectorType::SizeType SizeType;

      using BaseType::blockSize;

      template< class Operation >
      struct CommDataHandle
      : public BaseType :: template CommDataHandle< Operation >
      {};

    protected:
      using BaseType :: asImp;

      typedef TemporaryLocalFunction< DiscreteFunctionSpaceType > TemporaryLocalFunctionType;

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
      DiscreteFunctionDefault ( const std::string &name, const DiscreteFunctionSpaceType &dfSpace );

      DiscreteFunctionDefault ( std::string name, std::shared_ptr< const DiscreteFunctionSpaceType > dfSpace );

      DiscreteFunctionDefault ( const ThisType& );
      DiscreteFunctionDefault ( ThisType && other );

    public:
      ThisType& operator= ( ThisType&& ) = delete;
      ThisType &operator= ( const ThisType& ) = delete;

      // Default Implementations
      // -----------------------

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name() const */
      const std::string &name () const { return name_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name() */
      std::string &name () { return name_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::order() */
      constexpr int order() const
      {
        return space().order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::continuous */
      bool continuous() const
      {
        return space().continuous();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::space() const */
      const DiscreteFunctionSpaceType &space () const { return *dfSpace_; }

      /** \brief obtain a reference to the underlying grid part */
      const GridPartType &gridPart () const { return space().gridPart(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      LocalFunctionType localFunction ( const EntityType &entity ) { return LocalFunctionType( asImp(), entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      const LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( asImp(), entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction() */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      LocalFunctionType localFunction () { return LocalFunctionType( asImp() ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction() */
      [[deprecated("Use {Const,Temporary,Mutable}LocalFunction and LocalContribution instead!")]]
      const LocalFunctionType localFunction () const { return LocalFunctionType( asImp() ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::clear() */
      void clear() { dofVector().clear(); }

      DofVectorType &dofVector() { return asImp().dofVector(); }
      const DofVectorType &dofVector() const { return asImp().dofVector(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::blocks() */
      int blocks() const { return dofVector().size(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::block( unsigned int index ) */
      DofBlockPtrType block ( unsigned int index )
      {
        return dofVector().blockPtr( index );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::block( unsigned int index ) const */
      ConstDofBlockPtrType block ( unsigned int index ) const
      {
        return dofVector().blockPtr( index );
      }

      /** \brief Return the number of blocks in the block vector
       *
       *  \return Number of block in the block vector
       */
      SizeType size () const { return dofVector().size() * blockSize; }

      /** \brief Obtain the constant iterator pointing to the first dof
       *
       *  \return Constant iterator pointing to the first dof
       */
      ConstDofIteratorType dbegin () const { return dofVector().begin(); }

      /** \brief Obtain the non-constant iterator pointing to the first dof
       *
       *  \return Non-Constant iterator pointing to the first dof
       */
      DofIteratorType dbegin () { return dofVector().begin(); }

      /** \brief Obtain the constant iterator pointing to the last dof
       *
       *  \return Constant iterator pointing to the last dof
       */
      ConstDofIteratorType dend () const { return dofVector().end(); }

      /** \brief Obtain the non-constant iterator pointing to the last dof
       *
       *  \return Non-Constant iterator pointing to the last dof
       */
      DofIteratorType dend () { return dofVector().end(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy(const RangeFieldType &s,const DiscreteFunctionInterfaceType &g) */
      template <class DFType>
      void axpy ( const RangeFieldType &s, const DiscreteFunctionInterface< DFType > &g );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy(const RangeFieldType &s,const DiscreteFunctionInterfaceType &g) */
      void axpy ( const RangeFieldType &s, const DiscreteFunctionInterfaceType& g )
      {
        dofVector().axpy( s, g.dofVector() );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::scalarProductDofs */
      template <class DFType>
      RangeFieldType scalarProductDofs ( const DiscreteFunctionInterface< DFType > &other ) const
      {
        return scalarProduct_.scalarProductDofs( *this, other );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::normSquaredDofs */
      typename Dune::FieldTraits< RangeFieldType >::real_type normSquaredDofs ( ) const
      {
        return std::real( (*this).scalarProductDofs( *this ));
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::print */
      void print ( std :: ostream &out ) const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofsValid */
      inline bool dofsValid () const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
      template <class DFType>
      void assign ( const DiscreteFunctionInterface< DFType > &g );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
      void assign ( const DiscreteFunctionType &g )
      {
        dofVector() = g.dofVector();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dataHandle */
      template< class Operation >
      typename CommDataHandle< Operation >::Type dataHandle ( const Operation &operation );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::communicate() */
      void communicate()
      {
        // only call in single thread mode
        if( ! Fem :: MPIManager :: singleThreadMode() )
        {
          assert( Fem :: MPIManager :: singleThreadMode() );
          DUNE_THROW(InvalidStateException,"DiscreteFunctionInterface::communicate: only call in single thread mode!");
        }

        this->space().communicate( asImp() );
      }

      /** \copydoc Dune::Fem::Function::evaluate(const DomainType &x,RangeType &value) const */
      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        asImp().evaluateGlobal( x, [ &value ] ( const LocalCoordinateType &x, const TemporaryLocalFunctionType &localFunction )
                                              { localFunction.evaluate( x, value ); } );
      }

      /** \copydoc Dune::Fem::Function::jacobian(const DomainType &x,JacobianRangeType &jacobian) const */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        asImp().evaluateGlobal( x, [ &jacobian ] ( const LocalCoordinateType &x, const TemporaryLocalFunctionType &localFunction )
                                                 { localFunction.jacobian( x, jacobian ); } );

      }

      /** \copydoc Dune::Fem::Function::hessian (const DomainType &x,HessianRangeType &hessian) const */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        asImp().evaluateGlobal( x, [ &hessian ] ( const LocalCoordinateType &x, const TemporaryLocalFunctionType &localFunction )
                                                { localFunction.hessian( x, hessian ); } );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterface< DFType > &g) */
      template <class DFType>
      DiscreteFunctionType &operator+=(const DiscreteFunctionInterface< DFType > &g);

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterface< DFType > &g) */
      DiscreteFunctionType &operator+=(const DiscreteFunctionType &g)
      {
        dofVector() += g.dofVector();
        return asImp();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator-=(const DiscreteFunctionInterface< DFType > &g) */
      template <class DFType>
      DiscreteFunctionType &operator-=(const DiscreteFunctionInterface< DFType > &g);

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator-=(const DiscreteFunctionInterface< DFType > &g) */
      DiscreteFunctionType &operator-=(const DiscreteFunctionType& g)
      {
        dofVector() -= g.dofVector();
        return asImp();
      }

      /** \brief multiply all DoFs with a scalar factor
       *
       *  \param[in]  scalar  factor to multiply DoFs with
       *
       *  \returns reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator*= ( const RangeFieldType &scalar )
      {
        dofVector() *= scalar;
        return asImp();
      }

      /** \brief devide all DoFs by a scalar factor
       *
       *  \param[in]  scalar  factor with which all dofs are devided
       *
       *  \returns reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator/= ( const RangeFieldType &scalar )
      {
        return BaseType :: operator*=( RangeFieldType(1 ) / scalar );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::read */
      template< class StreamTraits >
      inline void read ( InStreamInterface< StreamTraits > &in );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::write */
      template< class StreamTraits >
      inline void write ( OutStreamInterface< StreamTraits > &out ) const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression()
       *
       *  \note The default implementation does nothing.
       */
      void enableDofCompression ()
      {}


      /** \copydoc Dune::Fem::DiscreteFunctionInterface::addScaledLocalDofs */
      template< class LocalDofs >
      void addScaledLocalDofs ( const EntityType &entity, const RangeFieldType &s, const LocalDofs &localDofs )
      {
        LeftAddScaled< const LocalDofs, const RangeFieldType > assignFunctor( localDofs, s );
        space().blockMapper().mapEach( entity, dofBlockFunctor( dofVector(), assignFunctor ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::addLocalDofs */
      template< class LocalDofs >
      void addLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        LeftAdd< const LocalDofs > assignFunctor( localDofs );
        space().blockMapper().mapEach( entity, dofBlockFunctor( dofVector(), assignFunctor ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::setLocalDofs */
      template< class LocalDofs >
      void setLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        LeftAssign< const LocalDofs > assignFunctor( localDofs );
        space().blockMapper().mapEach( entity, dofBlockFunctor( dofVector(), assignFunctor ) );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::getLocalDofs */
      template< class Vector >
      void getLocalDofs ( const EntityType &entity, Vector &localDofs ) const
      {
        AssignFunctor< Vector > assignFunctor( localDofs );
        space().blockMapper().mapEach( entity, dofBlockFunctor( dofVector(), assignFunctor ) );
      }

      template <class DFType>
      inline bool compare ( const DiscreteFunctionInterface< DFType> &g ) const;

      // Non-Interface Methods
      // ---------------------

      /** \brief obtain the local function storage
       *
       *  \returns a reference to the local function storage
       */
      LocalDofVectorAllocatorType &localDofVectorAllocator () const
      {
        return ldvAllocator_;
      }

      /** \brief Initiate the assemble of values using the LocalContribution concept
       *  \tparam AssembleOperation the specific operation (Add, Set, ...)
       */
      template< class AssembleOperation >
      void beginAssemble ()
      {
        const std::type_index id( typeid( AssembleOperation ) );
        if( assembleOperation_ != id )
        {
          if( assembleOperation_ != std::type_index( typeid( void ) ) )
            DUNE_THROW( InvalidStateException, "Another assemble operation in progress" );
          assembleOperation_ = id;
          assert( assembleCount_ == 0 );
          AssembleOperation::begin( asImp() );
        }
        ++assembleCount_;
      }

      /** \brief Finalize the assemble of values using the LocalContribution concept
       *  \tparam AssembleOperation the specific operation (Add, Set, ...)
       */
      template< class AssembleOperation >
      void endAssemble ()
      {
        const std::type_index id( typeid( AssembleOperation ) );
        if( assembleOperation_ != id )
          DUNE_THROW( InvalidStateException, "Assemble operation not in progress" );
        assert( assembleCount_ > 0 );
        if( --assembleCount_ == 0 )
        {
          AssembleOperation::end( asImp() );
          assembleOperation_ = std::type_index( typeid( void ) );
        }
      }

      //! get local Dofs and store a reference to it in the LocalDofVector
      void getLocalDofReferences ( const EntityType &entity, LocalDofVectorType &localDofs )
      {
        AssignVectorReference< LocalDofVectorType > assignFunctor( localDofs );
        space().blockMapper().mapEach( entity, dofBlockFunctor( dofVector(), assignFunctor ) );
      }

    protected:
      /** \copydoc Dune::PersistentObject::backup */
      virtual void backup() const
      {
        // get backup stream from persistence manager and write to it
        write( PersistenceManager :: backupStream() );
      }

      /** \copydoc Dune::PersistentObject::restore */
      virtual void restore()
      {
        // get restore stream from persistence manager and read from it
        read( PersistenceManager :: restoreStream() );
      }

      /** \copydoc Dune::PersistentObject::insertSubData */
      virtual void insertSubData();

      /** \copydoc Dune::PersistentObject::removeSubData */
      virtual void removeSubData();

      /** \brief evaluate functor in global coordinate */
      template< class Functor >
      void evaluateGlobal ( const DomainType &x, Functor functor ) const;

      // only PersistenceManager should call backup and restore
      friend class PersistenceManager;

      std::shared_ptr< const DiscreteFunctionSpaceType > dfSpace_;

      // the local function storage
      typename Traits :: LocalDofVectorStackType ldvStack_;
      mutable LocalDofVectorAllocatorType ldvAllocator_;

      mutable TemporaryLocalFunctionType localFunction_;

      std::string name_;
      ScalarProductType scalarProduct_;

      std::type_index assembleOperation_ = std::type_index( typeid( void ) );;
      std::size_t assembleCount_ = 0;
    }; // end class DiscreteFunctionDefault

    template< class ImplX, class ImplY >
    inline bool operator== ( const DiscreteFunctionInterface< ImplX > &x,
                             const DiscreteFunctionInterface< ImplY > &y )
    {
      return x.compare(y);
    }
    template< class ImplX, class ImplY >
    inline bool operator!= ( const DiscreteFunctionInterface< ImplX > &x,
                             const DiscreteFunctionInterface< ImplY > &y )
    {
      return !x.compare(y);
    }

    template< class DiscreteFunction >
    class ManagedDiscreteFunction;

    template< class DiscreteFunction >
    struct DiscreteFunctionTraits< ManagedDiscreteFunction< DiscreteFunction > >
    : public DiscreteFunctionTraits< DiscreteFunction > {};


    /** \class DiscreteFunctionTraits
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace, typename DofVector >
    struct DefaultDiscreteFunctionTraits
    {
      typedef DofVector DofVectorType;

      typedef DiscreteFunctionSpace                           DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainType  DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType   RangeType;

      typedef typename DofVectorType::IteratorType            DofIteratorType;
      typedef typename DofVectorType::ConstIteratorType       ConstDofIteratorType;
      typedef typename DofVectorType::DofBlockType            DofBlockType;
      typedef typename DofVectorType::ConstDofBlockType       ConstDofBlockType;
      typedef typename DofVectorType::DofBlockPtrType         DofBlockPtrType;
      typedef typename DofVectorType::ConstDofBlockPtrType    ConstDofBlockPtrType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType MapperType;
      typedef typename DofVectorType::FieldType DofType;

      typedef ThreadSafeValue< UninitializedObjectStack >         LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;
    };


  ///@}

  } // end namespace Fem

} // end namespace Dune

#include "discretefunction_inline.hh"

#include "gridfunctionadapter.hh"
#endif // #ifndef DUNE_FEM_DISCRETEFUNCTION_HH
