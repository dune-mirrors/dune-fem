#ifndef DUNE_FEM_ADAPTIVEFUNCTION_LOCALFUNCTION_HH
#define DUNE_FEM_ADAPTIVEFUNCTION_LOCALFUNCTION_HH

//- System includes
#include <string>

//- Dune includes
#include <dune/fem/space/common/dofstorage.hh>

#include <dune/fem/function/common/localfunction.hh>

namespace Dune
{
  
  //- Forward declarations
  template <class DiscreteFunctionSpaceImp>
  class AdaptiveDiscreteFunction;

  //- Forward declarations of Combined Space 
  template <class, int , DofStoragePolicy > 
  class CombinedSpace;



  /*  \class AdaptiveLocalFunction
   *  \brief Local function belonging to AdaptiveDiscreteFunction
   */
  template< class DiscreteFunctionSpaceImp >
  class AdaptiveLocalFunction
  : public LocalFunctionDefault
    < DiscreteFunctionSpaceImp,
      AdaptiveLocalFunction<DiscreteFunctionSpaceImp >
    >
  {
  public:
    //! type of  discrete function space this local function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef AdaptiveLocalFunction< DiscreteFunctionSpaceType > ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;
    
    friend class AdaptiveFunctionImplementation< DiscreteFunctionSpaceType >;

  public:
    //! type of discrete function this local function belongs to
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
      DiscreteFunctionType;
    
    //! Traits class with type definitions for AdaptiveDiscreteFunction and
    //! AdaptiveLocalFunction
    typedef AdaptiveDiscreteFunctionTraits< DiscreteFunctionSpaceType > Traits;

    //! Traits class of DiscreteFunctionSpaceType
    typedef typename DiscreteFunctionSpaceType :: Traits SpaceTraits;

    //! GridType
    typedef typename Traits :: GridType GridType;
    
    //! Function space type
    typedef typename SpaceTraits :: FunctionSpaceType FunctionSpaceType;
    
    //! Intrinsic data type for range field
    typedef typename Traits::RangeFieldType RangeFieldType;
    //! Vector type for the domain field
    typedef typename Traits::DomainType DomainType;
    //! Vector type for the range field
    typedef typename Traits::RangeType RangeType;
    //! Tensor type for the jacobian
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    //! Intrinsic data type for the degrees of freedom (dof)
    typedef RangeFieldType DofType;

    //! Container class type for the dofs
    typedef typename Traits::DofStorageType DofStorageType;

    //! The base function set of DiscreteFunctionSpaceType
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    //! Dimension of the range field
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };

    //! type of codim 0 entity
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

  private:  
    typedef typename GridType :: ctype ctype;
    enum { dim = GridType :: dimension };
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;

  public:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    inline AdaptiveLocalFunction( const DiscreteFunctionSpaceType &spc,
                                  DofStorageType& dofVec );
    
    //! Copy constructor
    inline AdaptiveLocalFunction( const ThisType &other );

    //! Destructor
    inline ~AdaptiveLocalFunction ();

    //- Operators
    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) */
    inline RangeFieldType &operator[] ( const int num );

    /** \copydoc Dune::LocalFunctionDefault::operator[](const int num) const */
    inline const RangeFieldType &operator[] ( const int num ) const;

    //- Methods

    /** \copydoc LocalFunctionInterface::numDofs */
    inline int numDofs() const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate( const DomainType &x,
                          RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian( const DomainType &x,
                          JacobianRangeType &ret ) const; 

    /** \copydoc Dune::LocalFunctionInterface::jacobian( const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian( const QuadratureType &quadrature,
                          const int quadPoint,
                          JacobianRangeType &ret ) const;
    
    /** \copydoc Dune::LocalFunctionDefault::baseFunctionSet */
    const BaseFunctionSetType& baseFunctionSet() const;

    /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor);

    /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const JacobianRangeType& factor);

    /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor1, const JacobianRangeType& factor2);

  private:
    //- Forbidden methods
    //! assignment operator
    ThisType& operator=(const ThisType& other);

  public:
    //! init local function 
    inline void init ( const EntityType &entity );

  private:
    inline
    void rightMultiply(const JacobianRangeType& factor, 
           const JacobianInverseType& jInv, 
                       JacobianRangeType& result) const;
    // return reference to actual entity
    const EntityType& en() const;

    //- Data members
    const DiscreteFunctionSpaceType& spc_;
  protected:
    DofStorageType& dofVec_;
    
    // vector holding pointer to local dofs 
    mutable MutableArray<RangeFieldType*> values_;
    
      //! number of local dofs 
    int numDofs_;
 
    mutable RangeType tmp_;
    mutable JacobianRangeType tmpGrad_;
    mutable JacobianRangeType factorInv_;

    mutable bool needCheckGeometry_;

    // base function set 
    mutable BaseFunctionSetType baseSet_;

    // actual entity
    mutable const EntityType* en_;
  }; // end class AdaptiveLocalFunction



  /*  \copydoc Dune::AdaptiveLocalFunction
   *
   *  Specialised version of AdaptiveLocalFunction for CombinedSpace
   */
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  class AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >
      : public LocalFunctionDefault<
    CombinedSpace<ContainedFunctionSpaceImp, N, p>,
    AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> > >
  {
  public:
    //- Friends
    friend class AdaptiveFunctionImplementation<
      CombinedSpace<ContainedFunctionSpaceImp, N, p> >;

  public:
    //- Public typedefs and enums
    typedef CombinedSpace<
      ContainedFunctionSpaceImp, N, p> DiscreteFunctionSpaceType;

  private:
    typedef AdaptiveLocalFunction< DiscreteFunctionSpaceType > ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType >
      BaseType;

  public:
    typedef AdaptiveDiscreteFunctionTraits<
      DiscreteFunctionSpaceType > Traits;
    typedef typename DiscreteFunctionSpaceType::Traits SpaceTraits;

    typedef typename Traits::GridType GridType;
    
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    
    typedef typename SpaceTraits::ContainedRangeType ContainedRangeType;
    typedef typename SpaceTraits::ContainedJacobianRangeType 
    ContainedJacobianRangeType;
    typedef typename SpaceTraits::BaseFunctionSetType BaseFunctionSetType;
       
    typedef typename Traits::RangeFieldType RangeFieldType;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::DofType DofType;
    
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef typename Traits::DofStorageType DofStorageType;

    typedef typename DiscreteFunctionSpaceType::GridType::
    template Codim<0>:: Entity EntityType;

  private:
    //- Typedefs
    typedef typename FieldVector<DofType*, N>::size_type SizeType;

    // store pointer to dofs belonging to local function 
    template <class Dummy, DofStoragePolicy policy>
    struct MapLocalDofs
    {
      template <class LocalDofVector>
      static void map(const DiscreteFunctionSpaceType& spc,
                      const EntityType& en, 
                      const int i, 
                      DofStorageType& dofVec,
                      LocalDofVector& values)
      {
        for (SizeType j = 0; j < N; ++j )
        {
          values[i][j] = &(dofVec[spc.mapToGlobal(en, i*N+j)]);
        } // end for j
      }
    };
    
    // specialisation for point based storage 
    template <class Dummy>
    struct MapLocalDofs<Dummy,PointBased>
    {
      template <class LocalDofVector>
      static void map(const DiscreteFunctionSpaceType& spc,
                      const EntityType& en, 
                      const int i, 
                      DofStorageType& dofVec,
                      LocalDofVector& values)
      {
        const int idx = spc.mapToGlobal(en, i*N );
        for (SizeType j = 0; j < N; ++j )
        {
          // only works for point base mappings 
          assert( (idx+ (int) j) == spc.mapToGlobal(en, i*N+j) );
          values[i][j] = &(dofVec[ idx+j ]);
        } // end for j
      }
    };
    
  private:  
    typedef typename GridType :: ctype ctype;
    enum { dim = GridType :: dimension };
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;

  public:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //- Public methods
    //- Constructors and destructors
    
    //! Constructor
    inline
    AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                          DofStorageType& dofVec);

    //! Copy constructor
    inline
    AdaptiveLocalFunction(const ThisType& other);
    
    //! Destructor
    ~AdaptiveLocalFunction();

    //- Operators
    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) const */
    inline const RangeFieldType &operator[]( const int num ) const;

    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) */
    inline RangeFieldType &operator[]( const int num );

    /** \copydoc Dune::LocalFunctionInterface::numDofs */
    inline int numDofs() const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint, 
                           RangeType &ret ) const;
    
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian ( const DomainType &x,
                           JacobianRangeType &ret ) const;
        
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian ( const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &ret ) const;

    //- Additional methods for specialisation
    /** \todo Please doc me!
     *
     *  \note Neither LocalFunctionInterface nor LocalFunctionDefault have this
     *        method.
     */
    inline
    void assign(int dofNum, const RangeType& dofs);
    
    /** \brief Number of contained scalar base functions */
    inline
    int numDifferentBaseFunctions() const;

    /** \copydoc Dune::LocalFunctionInterface::baseFunctionSet */
    inline const BaseFunctionSetType &baseFunctionSet() const;
   
     /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor);


    /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const JacobianRangeType& factor);

   
     /** \copydoc LocalFunctionDefault::axpy */
    template <class QuadratureType>
    inline void axpy(const QuadratureType&, const int qp, const RangeType& factor1, const JacobianRangeType& factor2);

  public:
    inline void init ( const EntityType &entity );

  private:
    //- Private methods
    // return reference to entity 
    inline
    const EntityType& en() const;

    // apply factor.rightmultiply(jInv) and strore in result 
    inline
    void rightMultiply(const JacobianRangeType& factor, const JacobianInverseType& jInv, 
                       JacobianRangeType& result) const;
  private:
    //- Member data
    const DiscreteFunctionSpaceType& spc_;

  protected:
    DofStorageType& dofVec_;
    
  private:
    // local dofs 
    typedef DofType* DofPointerType; 
    mutable MutableArray< DofPointerType[N] > values_;

    // number of local dofs 
    int numDofs_;
 
    // temporary variables 
    mutable ContainedRangeType cTmp_;
    mutable ContainedJacobianRangeType cTmpGradRef_;
    mutable ContainedJacobianRangeType cTmpGradReal_;
    mutable JacobianRangeType factorInv_;
    mutable RangeType tmp_;

    mutable bool needCheckGeometry_;
    mutable BaseFunctionSetType baseSet_;

    mutable const EntityType* en_;
  }; // end class AdaptiveLocalFunction (specialised for CombinedSpace)



  template< class DiscreteFunctionSpaceImp >
  class AdaptiveLocalFunctionFactory
  {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef AdaptiveLocalFunctionFactory< DiscreteFunctionSpaceType > ThisType;

    friend class AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >;

  public:
    typedef AdaptiveLocalFunction< DiscreteFunctionSpaceType > ObjectType;

    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
      DiscreteFunctionType;

  protected:
    DiscreteFunctionType &discreteFunction_;

  protected:
    inline explicit AdaptiveLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {
    }

  public:
    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_.space(),
                             discreteFunction_.dofVec_ );
    }
  };
  
} // end namespace Dune

#include "localfunction_inline.hh"

#endif
