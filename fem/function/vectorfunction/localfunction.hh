#ifndef DUNE_FEM_VECTORFUNCTION_LOCALFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_LOCALFUNCTION_HH

#include <dune/fem/function/common/localfunction.hh>

namespace Dune
{

  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  class VectorDiscreteFunction;



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  class VectorLocalFunction
  : public LocalFunctionDefault
    < DiscreteFunctionSpaceImp,
      VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    >
  {
  public:
    //! type fo the associated discrete function space
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef DofVectorImp DofVectorType;

  private:
    typedef VectorLocalFunction< DiscreteFunctionSpaceType, DofVectorType >
      ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType >
      BaseType;

  public:
    //! type of the discrete function, this local function belongs to
    typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >
      DiscreteFunctionType;

    //! type of base function sets
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    //! type of jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! field type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType
      DomainFieldType;
    //! field type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;

    enum
    {
      //! dimension of domain
      DimDomain = DiscreteFunctionSpaceType :: DimDomain,
      //! dimension of range
      DimRange = DiscreteFunctionSpaceType :: DimRange
    };

    //! type of DoFs
    typedef typename DofVectorType :: FieldType DofType;

  protected:
    //! type of the underlying grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity Codim0EntityType;

    //! type of geometry
    typedef typename Codim0EntityType :: Geometry GeometryType;

    //! type of geometry's jacobian inverse
    typedef FieldMatrix< typename GeometryType :: ctype,
                         GeometryType :: mydimension, GeometryType :: mydimension >
      GeometryJacobianInverseType;

  protected:
    DiscreteFunctionType &discreteFunction_;

    BaseFunctionSetType baseFunctionSet_;
    const GeometryType *geometry_;
    bool needCheckGeometry_;

    DynamicArray< DofType* > values_;

  public:
    //! constructor
    inline VectorLocalFunction ( DiscreteFunctionType &df )
    : discreteFunction_( df ),
      geometry_( 0 ),
      needCheckGeometry_( true ),
      values_( 0 )
    {
    }

    //! copy constructor
    inline VectorLocalFunction ( const ThisType &other )
    : discreteFunction_( other.discreteFunction_ ),
      baseFunctionSet_( other.baseFunctionSet_ ),
      geometry_( other.geometry_ ),
      needCheckGeometry_( other.needCheckGeometry_ ),
      values_( other.values_ )
    {
    }

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType &other );

  public:
    inline const DofType &operator[] ( const int dof ) const
    {
      return *(values_[ dof ]);
    }

    inline DofType &operator[] ( const int dof )
    {
      return *(values_[ dof ]);
    }

    inline const BaseFunctionSetType &baseFunctionSet () const
    {
      return baseFunctionSet_;
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const
     */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const;

    inline void init ( const Codim0EntityType &entity );

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const
     */
    inline void jacobian ( const DomainType &x,
                           JacobianRangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const
     */
    template< class QuadratureType >
    inline void jacobian ( const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &ret ) const;

    /** \copydoc Dune::LocalFunctionInterface::numDofs
     */
    inline int numDofs () const
    {
      return values_.size();
    }
  };



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  class VectorLocalFunctionFactory
  {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef DofVectorImp DofVectorType;

  private:
    typedef VectorLocalFunctionFactory
      < DiscreteFunctionSpaceType, DofVectorType >
      ThisType;
  
    friend class VectorDiscreteFunction
      < DiscreteFunctionSpaceType, DofVectorType >;

  public:
    typedef VectorLocalFunction< DiscreteFunctionSpaceType, DofVectorType >
      ObjectType;

    typedef VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >
      DiscreteFunctionType;

  protected:
    DiscreteFunctionType &discreteFunction_;

  protected:
    inline explicit VectorLocalFunctionFactory ( DiscreteFunctionType &df )
    : discreteFunction_( df )
    {
    }

  public:
    ObjectType *newObject () const
    {
      return new ObjectType( discreteFunction_ );
    }
  };

}

#include "localfunction_inline.hh"

#endif
