#ifndef DUNE_FEM_DISCRETEFUNCTION_HH
#define DUNE_FEM_DISCRETEFUNCTION_HH

// C++ includes
#include <string>
#include <iostream>
#include <fstream>

// dune-common inlcudes
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/debug.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
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
    {};

    /** Base class for determing whether a function has local functions or not.
    */
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
        features of a discrete function.
        It contains a local function and a dof iterator which can
        iterate over all dofs of one level. Via the method access the local
        dofs and basis functions can be accessed for a given entity.
        The DOF-Iterators are STL-like Iterators, i.e. they can be dereferenced
        giving the corresponding DOF.

        \interfaceclass
    */
    template< class Impl >
    class DiscreteFunctionInterface
    : public Function< typename DiscreteFunctionTraits< Impl >::DiscreteFunctionSpaceType::FunctionSpaceType, Impl >,
      public IsDiscreteFunction,
      public HasLocalFunction
    {
      typedef DiscreteFunctionInterface< Impl > ThisType;
      typedef Function< typename DiscreteFunctionTraits< Impl >::DiscreteFunctionSpaceType::FunctionSpaceType, Impl > BaseType;

    public:
      //! type of the traits
      typedef DiscreteFunctionTraits< Impl > Traits;

      //! type of the implementaton (Barton-Nackman)
      typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;

      //! type of associated discrete function space
      typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

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

      //! Type of the underlying grid
      typedef typename DiscreteFunctionSpaceType :: GridType GridType;

      //! type of local functions
      typedef typename Traits :: LocalFunctionType LocalFunctionType;

      //! Type of the dof iterator used in the discrete function implementation.
      typedef typename Traits :: DofIteratorType DofIteratorType;

      //! Type of the constantdof iterator used in the discrete function implementation
      typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

      typedef typename Traits :: DofType DofType;
      typedef typename Traits :: DofBlockType DofBlockType;
      typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
      typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
      typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

      //! type of mapping base class for this discrete function
      typedef typename BaseType :: MappingType MappingType;

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

      /** \brief default constructor */
      DiscreteFunctionInterface ()
      {}

    private:
      // prohibit copying and assignment
      DiscreteFunctionInterface ( const ThisType &other );
      ThisType &operator= ( const ThisType &other );

    public:
      /** \brief obtain the name of the discrete function
       *
       *  \returns string holding name of discrete function
       */
      const std::string &name () const
      {
        return asImp().name();
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
      LocalFunctionType localFunction ( const EntityType &entity )
      {
        return asImp().localFunction( entity );
      }

      /** \brief obtain a local function for an entity (read-write)
       *
       *  \param[in]  entity  Entity to focus view of discrete function
       *  \returns a local function associated with the entity
       */
      const LocalFunctionType localFunction ( const EntityType &entity ) const
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

      /** \brief obtain pointer to block of discrete function with block number index
       *
       *  \returns a DoFBlockPtrType pointing to block index
       */
      inline ConstDofBlockPtrType block ( unsigned int index ) const
      {
        return asImp().block( index );
      }

      /** \brief obtain pointer to block of discrete function with block number index (read-only)
       *
       *  \returns a DoFBlockPtrType pointing to block index (read-only)
       */
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

      /** \brief dump all DoFs into a file (for debugging purposes)
       *
       *  \param[in]  filename  name of file to dump into
       */
      inline void print( const std :: string &filename ) const
      {
        asImp().print( filename );
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
      void assign( const DiscreteFunctionInterfaceType &g )
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

      /** \brief do default communication of space for this discrete
                 function
      */
      inline void communicate()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().communicate() );
      }

      /** \brief add another discrete function to this one
       *
       *  \param[in]  g  discrete function to add
       *
       *  \returns a reference to this discrete function (i.e. *this)
       */
      DiscreteFunctionType &operator+= ( const DiscreteFunctionInterfaceType &g )
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

      /** \brief Enable this discrete function for dof compression,
           i.e. during grdi changes a dof compression
           is done when the DofManagers compress is called.
      */
      inline void enableDofCompression()
      {
        asImp().enableDofCompression();
      }

      // this needs to be revised, the definition should be in GridPart
      // further discussion needed
      typedef LoadBalanceLeafData< ThisType > DefaultLoadBalanceContainsCheckType;
      inline DefaultLoadBalanceContainsCheckType defaultLoadBalanceContainsCheck() const
      {
        return DefaultLoadBalanceContainsCheckType( *this );
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

      //! type of LocalDofVector
      typedef typename Traits :: LocalDofVectorType LocalDofVectorType;
      //! type of LocalDofVector
      typedef typename Traits :: LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      //! type of local functions
      typedef typename BaseType :: LocalFunctionType LocalFunctionType;

      typedef typename BaseType :: DofBlockType DofBlockType;
      typedef typename BaseType :: ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType :: DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType :: ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType :: EntityType EntityType ;

      typedef typename BaseType :: DofType DofType;

      template< class Operation >
      struct CommDataHandle
      : public BaseType :: template CommDataHandle< Operation >
      {};

    private:
      struct LocalFunctionEvaluateFunctor
      {
        typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
        typedef typename LocalFunctionType::RangeType RangeType;

        LocalFunctionEvaluateFunctor ( RangeType &value ) : value_( value ) {}

        void operator() ( const LocalCoordinateType &x, const LocalFunctionType &localFunction )
        {
          localFunction.evaluate( x, value_ );
        }

      private:
        RangeType &value_;
      };

      struct LocalFunctionJacobianFunctor
      {
        typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
        typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

        LocalFunctionJacobianFunctor ( JacobianRangeType &jacobian ) : jacobian_( jacobian ) {}

        void operator() ( const LocalCoordinateType &x, const LocalFunctionType &localFunction )
        {
          localFunction.jacobian( x, jacobian_);
        }

      private:
        JacobianRangeType &jacobian_;
      };

      struct LocalFunctionHessianFunctor
      {
        typedef typename LocalFunctionType::LocalCoordinateType LocalCoordinateType;
        typedef typename LocalFunctionType::HessianRangeType HessianRangeType;

        LocalFunctionHessianFunctor ( HessianRangeType &hessian ) : hessian_( hessian ) {}

        void operator() ( const LocalCoordinateType &x, const LocalFunctionType &localFunction )
        {
          localFunction.hessian( x, hessian_ );
        }

      private:
        HessianRangeType &hessian_;
      };

    protected:
      using BaseType :: asImp;

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
      DiscreteFunctionDefault ( const std::string &name,
                                const DiscreteFunctionSpaceType &dfSpace,
                                const LocalDofVectorAllocatorType &ldvAllocator );
    private:
      // prohibit copying and assignment
      inline DiscreteFunctionDefault ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    public:
      // Default Implementations
      // -----------------------

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name() const */
      const std::string &name () const { return name_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::space() const */
      const DiscreteFunctionSpaceType &space () const { return dfSpace_; }

      /** \brief obtain a reference to the underlying grid part */
      const GridPartType &gridPart () const { return space().gridPart(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      LocalFunctionType localFunction ( const EntityType &entity ) { return LocalFunctionType( asImp(), entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      const LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( asImp(), entity ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::clear() */
      void clear();

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::allocDofPointer
       *
       *  \note The default implementation make a copy of the DoF vector using
       *        the DoF iterators.
       */
      inline RangeFieldType *allocDofPointer ();

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::freeDofPointer
       *
       *  \note The default implementation make a copy of the DoF vector using
       *        the DoF iterators.
       */
      inline void freeDofPointer( RangeFieldType *dofPointer );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::axpy(const RangeFieldType &s,const DiscreteFunctionInterfaceType &g) */
      void axpy ( const RangeFieldType &s, const DiscreteFunctionInterfaceType &g );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::scalarProductDofs */
      inline RangeFieldType
      scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const { return scalarProduct_.scalarProductDofs( *this, other ); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::print */
      void print ( std :: ostream &out ) const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::print */
      void print ( const std :: string &filename ) const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofsValid */
      inline bool dofsValid () const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g) */
      void assign ( const DiscreteFunctionInterfaceType &g );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dataHandle */
      template< class Operation >
      typename CommDataHandle< Operation > :: Type
      dataHandle ( const Operation *operation );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::communicate() */
      void communicate()
      {
        assert( Fem :: ThreadManager :: singleThreadMode() );
        this->space().communicate( asImp() );
      }

      /** \copydoc Dune::Fem::Function::evaluate(const DomainType &x,RangeType &value) const */
      inline void evaluate ( const DomainType &x, RangeType &value ) const
      {
        LocalFunctionEvaluateFunctor functor( value );
        asImp().evaluateGlobal( x, functor );
      }

      /** \copydoc Dune::Fem::Function::jacobian(const DomainType &x,JacobianRangeType &jacobian) const */
      inline void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        LocalFunctionJacobianFunctor functor( jacobian );
        asImp().evaluateGlobal( x, functor );
      }

      /** \copydoc Dune::Fem::Function::hessian (const DomainType &x,HessianRangeType &hessian) const */
      inline void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        LocalFunctionHessianFunctor functor( hessian );
        asImp().evaluateGlobal( x, functor );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::operator+=(const DiscreteFunctionInterfaceType &g) */
      DiscreteFunctionType &operator+= ( const DiscreteFunctionInterfaceType &g );

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
      inline DiscreteFunctionType &operator/= ( const RangeFieldType &scalar ) { return BaseType :: operator*=( RangeFieldType(1 ) / scalar ); }

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
      inline void enableDofCompression () {}


    public:
      // Non-Interface Methods
      // ---------------------

      inline bool operator== ( const DiscreteFunctionType &g ) const;

      inline bool operator!= ( const DiscreteFunctionType &g ) const { return !(operator==( g )); }

      /** \brief obtain the local function storage
       *
       *  \returns a reference to the local function storage
       */
      inline LocalDofVectorAllocatorType &localDofVectorAllocator () const { return ldvAllocator_; }

      //! add local Dofs to dof vector
      template< class LocalDofs >
      void addLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        typedef LeftAdd< const LocalDofs > AssignFunctorType;
        AssignFunctorType assignFunctor( localDofs );

        DofBlockFunctor< DiscreteFunctionType, AssignFunctorType > functor( asImp(), assignFunctor );
        space().blockMapper().mapEach( entity, functor );
      }

      //! set local Dofs to dof vector
      template< class LocalDofs >
      void setLocalDofs ( const EntityType &entity, const LocalDofs &localDofs )
      {
        typedef LeftAssign< const LocalDofs > AssignFunctorType;
        AssignFunctorType assignFunctor( localDofs );

        DofBlockFunctor< DiscreteFunctionType, AssignFunctorType > functor( asImp(), assignFunctor );
        space().blockMapper().mapEach( entity, functor );
      }

      //! get local Dofs and store a reference to it in the LocalDofVector
      void getLocalDofs ( const EntityType &entity, LocalDofVectorType &localDofs )
      {
        typedef AssignVectorReference< LocalDofVectorType > AssignFunctorType;
        AssignFunctorType assignFunctor( localDofs );

        DofBlockFunctor< DiscreteFunctionType, AssignFunctorType > functor( asImp(), assignFunctor );
        space().blockMapper().mapEach( entity, functor );
      }

      //! get local Dofs and store the values  in LocalDofVector
      template< class A >
      void getLocalDofs ( const EntityType &entity, Dune::DynamicVector< DofType, A > &localDofs ) const
      {
        typedef AssignFunctor< Dune::DynamicVector< DofType, A > > AssignFunctorType;
        AssignFunctorType assignFunctor( localDofs );

        DofBlockFunctor< const DiscreteFunctionType, AssignFunctorType > functor( asImp(), assignFunctor );
        space().blockMapper().mapEach( entity, functor );
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

    private:
      const DiscreteFunctionSpaceType &dfSpace_;

      // the local function storage
      mutable LocalDofVectorAllocatorType ldvAllocator_;

      DebugLock dofPointerLock_;

    protected:
      std::string name_;
      ScalarProductType scalarProduct_;
    }; // end class DiscreteFunctionDefault


    template< class DiscreteFunction >
    class ManagedDiscreteFunction;

  ///@}

  } // end namespace Fem

} // end namespace Dune

#include "discretefunction_inline.hh"

#include "gridfunctionadapter.hh"
#endif // #ifndef DUNE_FEM_DISCRETEFUNCTION_HH
