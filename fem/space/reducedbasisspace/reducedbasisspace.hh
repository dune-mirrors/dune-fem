#ifndef DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH

#include "basefunctionset.hh"
#include "mapper.hh"
#include "commdatahandle.hh"

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune
{

  template< class BaseFunction >
  class ReducedBasisSpace;



  template< class BaseFunction >
  class ReducedBasisSpaceTraits
  {
  public:
    typedef BaseFunction BaseFunctionType;

  private:
    typedef ReducedBasisSpaceTraits< BaseFunctionType > ThisType;

  public:
    typedef ReducedBasisSpace< BaseFunctionType > DiscreteFunctionSpaceType;

    typedef typename BaseFunctionType :: FunctionSpaceType BaseFunctionSpaceType;

    typedef typename BaseFunctionSpaceType :: DomainType DomainType;
    typedef typename BaseFunctionSpaceType :: RangeType RangeType;
    typedef typename BaseFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    
    typedef typename BaseFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename BaseFunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef typename BaseFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
    typedef typename BaseFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseFunctionSpaceType :: GridType GridType;
    typedef typename BaseFunctionSpaceType :: IndexSetType IndexSetType;
    typedef typename BaseFunctionSpaceType :: IteratorType IteratorType;

    typedef ReducedBasisBaseFunctionSet< BaseFunctionType >
      BaseFunctionSetType;

    typedef typename BaseFunctionSetType :: BaseFunctionListType BaseFunctionListType;

    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType >
      MapperType;

    enum { localBlockSize = 1 };
    typedef MapperType BlockMapperType;

    template< class DiscreteFunction,
              class Operation = DFCommunicationOperation :: Add >
    struct CommDataHandle
    {
      // This is just a phony data handle doing nothing
      typedef ReducedBasisCommDataHandle< DiscreteFunction, Operation > Type;
      typedef Operation OperationType;
    };
  };



  /** \class   ReducedBasisSpace
   *  \ingroup RBSpace
   *  \brief   provides the space for reduced basis simulations 
   *
   *  The basis consists of discrete functions as basis functions. These
   *  discrete functions have an underlying arbitrary space. Consequently they
   *  inhert most of the structure from this space.
   *  
   *  Initially the space is empty and by using the add function you can bulid this
   *  space and discrete functions.
   *
   *  \param  BaseFunction  type of discrete function used to represend the
   *                        base functions
   */
  template< class BaseFunction >
  class ReducedBasisSpace
  : public DiscreteFunctionSpaceDefault< ReducedBasisSpaceTraits< BaseFunction > >
  {
  public:
    //! discrete function type of the base functions
    typedef BaseFunction BaseFunctionType;

    //! type of the traits
    typedef ReducedBasisSpaceTraits< BaseFunctionType > Traits;

  private:
    typedef ReducedBasisSpace< BaseFunctionType > ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  public:
    //! discrete function space, the base functions belong to
    typedef typename Traits :: BaseFunctionSpaceType BaseFunctionSpaceType;
    //! function space, the base functions belong to
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

    typedef typename Traits :: GridPartType GridPartType;
    typedef typename Traits :: GridType GridType;
    typedef typename Traits :: IndexSetType IndexSetType;

    //! type of the space's iterator
    typedef typename Traits :: IteratorType IteratorType;

    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    typedef typename Traits :: BaseFunctionListType BaseFunctionListType;

    //! type of the DoF mapper
    typedef typename Traits :: MapperType MapperType;

    enum { localBlockSize = Traits :: localBlockSize };
    //! type of the block mapper
    typedef typename Traits :: BlockMapperType BlockMapperType;

    enum { polynomialOrder = BaseFunctionSpaceType :: polynomialOrder };

  protected:
    BaseFunctionSpaceType &baseFunctionSpace_;
    BaseFunctionListType baseFunctionList_;
    mutable MapperType mapper_;

  public:
    //! default communication interface 
    static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

    //! default communication direction 
    static const CommunicationDirection defaultDirection =  ForwardCommunication;

    /** \brief contructor
     *
     *  \param[in]  baseFunctionSpace  DiscreteFunctionSpace containing the
     *                                 base functions belong to
     *  \param[in]  commInterface      communication interface 
     *  \param[in]  commDirection      communication direction 
     */
    inline explicit
    ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                        const InterfaceType commInterface = defaultInterface ,
                        const CommunicationDirection commDirection = defaultDirection );

    /** \brief contructor reading the base functions from a stream
     *
     *  \param[in]  baseFunctionSpace  DiscreteFunctionSpace containing the
     *                                 base functions belong to
     *  \param[in]  in                 stream to read the base functions from
     *  \param[in]  commInterface      communication interface 
     *  \param[in]  commDirection      communication direction 
     */
    template< class StreamTraits >
    inline
    ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                        InStreamInterface< StreamTraits > &in ,
                        const InterfaceType commInterface = defaultInterface ,
                        const CommunicationDirection commDirection = defaultDirection );

    inline ~ReducedBasisSpace ();


    // Implementation of DiscreteFunctionSpaceInterface
    // ------------------------------------------------
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous() const */
    inline bool continuous () const
    {
      return baseFunctionSpace_.continuous();
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::order() const */
    inline int order () const
    {
      return baseFunctionSpace_.order();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
    template< class EntityType >
    inline const BaseFunctionSetType baseFunctionSet( const EntityType &entity ) const
    {
      return BaseFunctionSetType( baseFunctionList_, entity );
    }

    //! get dimension of value
    inline int dimensionOfValue () const
    {
      return baseFunctionSpace_.dimensionOfValue;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper() const */
    inline MapperType &mapper () const
    {
      return mapper_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::blockMapper() const */
    inline BlockMapperType &blockMapper () const
    {
      return mapper_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::multipleBaseFunctionSets () const */
    inline bool multipleBaseFunctionSets () const
    {
      return true;
    }


    // ReducedBasisSpace Specific Methods
    // ----------------------------------

    /** \brief add a base function to the reduced basis space
     *
     *  \note After adding the base function, you can safely remove or
     *        overwrite it. The ReducedBasisSpace make itself a copy.
     * 
     *  \param[in]  baseFunction  base function to add to the reduced basis
     *                            space
     */
    inline void addBaseFunction ( const BaseFunctionType &baseFunction );

    /** \brief access a base function within the reduced basis space
     *
     *  \param[in]  i  number of the base function to access
     * 
     *  \returns a constant reference to the i-th base function
     */
    inline const BaseFunctionType &baseFunction ( unsigned int i ) const;

    inline const BaseFunctionSpaceType &baseFunctionSpace () const;

    /** \brief remove all base functions from the reduced basis space */
    inline void clear ();

    /** \brief crop base function set to the first n base functions
     *
     *  \param[in]  n  number of base functions to keep (must be less or equal
     *                 to the current number of base functions)
     */
    inline void crop ( unsigned int n );

    /** \brief obtain number of base functions within the reduced basis space */
    inline unsigned int numBaseFunctions () const;

    /** \brief project a discrete function over this space to the discrete
     *         function space of the base functions
     *
     *  \note This method expects the source discrete function to implement the
     *        dof method (which is not part of the DiscreteFunctionInterface).
     *
     *  \param[in]   sourceFunction  discrete function to be projected
     *  \param[out]  destFunction    discrete function to receive the projected
     *                               function
     */
    template< class DiscreteFunctionType >
    inline void project ( const DiscreteFunctionType &sourceFunction,
                          BaseFunctionType &destFunction ) const;

    /** \brief restrict a discrete function over the base space to the discrete
     *         function space of this base functions
     *
     *  \note This method expects the source discrete function to be a DiscreteFunction in the baseFunctionSpace
     *
     *  \param[in]   sourceFunction  discrete function to be restricted
     *  \param[out]  destFunction    discrete function to receive the restricted
     *                               function
     */
    template< class DiscreteFunctionType >
    inline void restrictFunction ( const BaseFunctionType &sourceFunction,
                                   DiscreteFunctionType &destFunction ) const;

    /** \brief restrict a bilinearform which operates in the high dimensional space to a lower dimensional matrix which 
     *         operates in the reducedbasisspace
     *
     *  \note This method expects the matrix of the FEM simulation
     *
     *  \param[in]   matrixOffline  matrix which corresponds the bilinearform
     *                              in the lagrange space
     *  \param[out]  matrixOnline   projected matrix to the lower dimensional
     *                              RB space
     */
    template< class MatrixOfflineType , class MatrixOnlineType >
    inline void restrictMatrix ( const MatrixOfflineType &matrixOffline,
                                 MatrixOnlineType &matrixOnline ) const;

    template< class StreamTraits >
    inline void read ( InStreamInterface< StreamTraits > &in );
   
    template< class StreamTraits >
    inline void write ( OutStreamInterface< StreamTraits > &out ) const;
  };
  
}

#include "reducedbasisspace_inline.hh"

#endif
