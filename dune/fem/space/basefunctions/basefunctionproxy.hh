#ifndef DUNE_BASEFUNCTIONPROXY_HH
#define DUNE_BASEFUNCTIONPROXY_HH

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

namespace Dune
{

  /** \addtogroup BaseFunction
   *  \{
   */
    
  template < class BaseFunctionSetImp > 
  class SimpleBaseFunctionProxy; 

  template <class BaseFunctionSetImp>
  struct SimpleBaseFunctionProxyTraits  
  {
    //! Export function space type
    typedef typename BaseFunctionSetImp :: Traits :: FunctionSpaceType FunctionSpaceType;

    //! Exact type of the base function
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
  };

  /** \class SimpleBaseFunctionProxy
   *  \brief proxy object for base function sets
   */
  template< class BaseFunctionSetImp > 
  class SimpleBaseFunctionProxy
  : public BaseFunctionSetDefault< SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > >
  {
  public:
    typedef SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > Traits;
    
  private:
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;

  public:
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

    enum { dimRange = FunctionSpaceType :: dimRange };
    enum { dimDomain = FunctionSpaceType :: dimDomain };

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  protected:
    // base function set 
    const BaseFunctionSetImp *baseSet_; 

  public:
    inline SimpleBaseFunctionProxy ()
    : baseSet_( NULL )
    {
    }
      
    inline SimpleBaseFunctionProxy ( const BaseFunctionSetImp* baseSet )
    : baseSet_( baseSet ) 
    {
    }

    //! Constructor creating empty local function 
    inline SimpleBaseFunctionProxy( const ThisType& org )
    : baseSet_( org.baseSet_ )
    {
    }

    //! asignment operator 
    SimpleBaseFunctionProxy &operator= ( const ThisType& org )
    {
      baseSet_ = org.baseSet_; 
      return *this;
    }   

    //! destructor 
    ~SimpleBaseFunctionProxy() 
    { 
      // unset pointer to be sure its gone 
      baseSet_ = 0;
    }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions() const 
    {
      return baseFunctionSet().numBaseFunctions(); 
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return baseFunctionSet().geometryType();
    }
   
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, diffVariable, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const PointType &x,RangeType &phi) const */
    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetDefault::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */
    template< class PointType >
    inline void jacobian( const int baseFunction,
                          const PointType &x,
                          JacobianRangeType &phi ) const
    {
      baseFunctionSet().jacobian( baseFunction, x, phi );
    }

    template< class PointType, 
              class LocalDofVectorType >
    inline void evaluateAll ( const PointType &x,
                              const LocalDofVectorType& dofs,
                              RangeType& ret) const
    {
      baseFunctionSet().evaluateAll( x, dofs, ret );
    }

    template< int diffOrder, 
              class PointType, 
              class LocalDofVectorType >
    inline void evaluateAll ( const FieldVector<int, diffOrder>& diffVariable,
                              const PointType &x,
                              const LocalDofVectorType& dofs,
                              RangeType& ret) const
    {
      baseFunctionSet().evaluateAll( diffVariable, x, dofs, ret );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    inline void jacobianAll ( const PointType &x,
                              const GeometryJacobianInverseType& gjit,   
                              const LocalDofVectorType& dofs,
                              GlobalJacobianRangeType& ret) const
    {
      baseFunctionSet().jacobianAll( x, gjit, dofs, ret );
    }

#if 0
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           const RangeType &psi ) const
    {
      return baseFunctionSet().evaluateSingle( baseFunction, x, psi );
    }

     /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const PointType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    inline RangeFieldType evaluateGradientSingle ( const int baseFunction,
                                                   const EntityType &entity,
                                                   const PointType &x,
                                                   const JacobianRangeType &psi ) const
    {
      return baseFunctionSet().evaluateGradientSingle( baseFunction, entity, x, psi );
    }
#endif

    template< class PointType, class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const RangeType& factor,
                       LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, factor, dofs );
    }

    template< class QuadratureType, 
              class LocalDofVectorType,
              class RangeVectorType>
    inline void 
    evaluateRanges ( const QuadratureType& quad,
                     const LocalDofVectorType& dofs,
                     RangeVectorType& rangeVector) const
    {
      baseFunctionSet().evaluateRanges( quad, dofs, rangeVector );
    }

    template< class QuadratureType, 
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    inline void 
    evaluateJacobians( const QuadratureType& quad,
                       const Geometry& geometry, 
                       const LocalDofVectorType& dofs,
                       JacobianRangeVectorType& jacVector) const
    {
      baseFunctionSet().evaluateJacobians( quad, geometry, dofs, jacVector );
    }

    template< class QuadratureType, 
              class RangeVectorType,
              class LocalDofVectorType >
    inline void axpyRanges ( const QuadratureType& quad,
                             const RangeVectorType& factorVec,
                             LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpyRanges( quad, factorVec, dofs );
    }

    template< class QuadratureType, 
              class Geometry,
              class JacobianRangeVectorType,
              class LocalDofVectorType >
    inline void axpyJacobians( const QuadratureType& quad,
                               const Geometry& geometry,
                               const JacobianRangeVectorType& factorVec,
                               LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpyJacobians( quad, geometry, factorVec, dofs );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const GeometryJacobianInverseType& gjit,   
                       const GlobalJacobianRangeType& factor,
                       LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, gjit, factor, dofs );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const GeometryJacobianInverseType& gjit,
                       const RangeType& factor1,
                       const GlobalJacobianRangeType& factor2,
                       LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, gjit, factor1, factor2, dofs );
    }

  protected:
    inline const BaseFunctionSetImp &baseFunctionSet () const
    {
      assert( baseSet_ );
      return *baseSet_;
    }
  }; // end SimpleBaseFunctionProxy 

  /** \} */

} // end namespace Dune

#endif
