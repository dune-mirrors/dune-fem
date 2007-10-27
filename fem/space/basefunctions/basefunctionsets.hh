#ifndef DUNE_BASEFUNCTIONSETS_HH
#define DUNE_BASEFUNCTIONSETS_HH

//- Dune includes 
#include <dune/common/fvector.hh>

//- local includes 
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/common/dofstorage.hh>

namespace Dune
{
  
  /** \addtogroup BaseFunction
   *  \{
   */
  
  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet;



  //! Traits class for standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StandardBaseFunctionSetTraits
  {
    //! Export function space type
    typedef FunctionSpaceImp FunctionSpaceType;
    //! Type of the base function storage policy
    typedef StorageImp<FunctionSpaceType> StorageType;
    //! Exact type of the base function
    typedef StandardBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    //! Factory type for the corresponding base functions (polymorphic)
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;
  };



  /** \class StandardBaseFunctionSet
   *  \brief standard base function set
   */
  template< class FunctionSpaceImp, template< class > class StorageImp >
  class StandardBaseFunctionSet
  : public BaseFunctionSetDefault
    < StandardBaseFunctionSetTraits< FunctionSpaceImp, StorageImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

    typedef StandardBaseFunctionSetTraits< FunctionSpaceType, StorageImp > Traits;

  private:
    typedef StandardBaseFunctionSet< FunctionSpaceType, StorageImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;
    
  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: RangeFieldType DofType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename Traits :: FactoryType FactoryType;

  public:
    // use evaluate of default implementation 
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! Constructor
    inline explicit StandardBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory )
      //diffVar0_( 0 ),
      //tmp_( 0 ),
      //jTmp_( 0 )
    {
    }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions () const
    {
      return storage_.numBaseFunctions();
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return storage_.geometryType();
    }
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const DomainType &x,RangeType &phi) const
     */ 
    template< int diffOrd >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const QuadratureType &quadrature,const int quadPoint,RangeType &phi) const
     */
    template< int diffOrd, class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      storage_.evaluate( baseFunction, diffVariable, quadrature, quadPoint, phi );
    }

#if 0
    /** \copydoc Dune::BaseFunctionSetDefault::evaluateSingle(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,const RangeType &psi) const */ 
    template< class QuadratureType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const QuadratureType &quadrature,
                                           const int quadPoint,
                                           const RangeType &psi ) const;
      
    /** \copydoc Dune::BaseFunctionSetDefault::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &psi) const */ 
    template <class EntityType, class QuadratureType>
    inline RangleFieldType evaluateGradientSingle ( const int baseFunction,
                                                    const EntityType &entity,
                                                    const QuadratureType &quadrature,
                                                    const int quadPoint,
                                                    const JacobianRangeType &psi ) const;
#endif
   
  private:
    StandardBaseFunctionSet( const StandardBaseFunctionSet& );

  protected:
    typename Traits::StorageType storage_;
    
    //mutable FieldVector<int, 0> diffVar0_;
    //mutable RangeType tmp_;
    //mutable JacobianRangeType jTmp_;
  };



  //- VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VectorialBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };



  //! \brief Special base function implementation that takes advantage
  //! of the vectorial structure of the base functions.
  //! This base function can be used in conjunction with scalar basefunctions
  //! \f$ \phi_i \f$ which are extended to vectorial base functions like 
  //! \f$ \Phi_j = \phi_i e_k \f$, where \f$ e_k = [ \kronecker_ik ]_i \f$.
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet : 
    public BaseFunctionSetDefault<VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  private:
    typedef BaseFunctionSetDefault<
      VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> > BaseType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
  public:  
    typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType::JacobianRangeType 
      ScalarJacobianRangeType;
  private:  
    typedef BaseFunctionFactory<ScalarFunctionSpaceType> FactoryType;
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;
    typedef StorageImp<ScalarFunctionSpaceType> StorageType;
 
  public:
    typedef VectorialBaseFunctionSetTraits<FunctionSpaceImp,StorageImp> Traits;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename FunctionSpaceType::RangeFieldType DofType;
        
  public:
    //! Constructor
    inline explicit VectorialBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory ),
      util_( FunctionSpaceType::DimRange ),
      tmp_( 0 ),
      jTmp_( 0 ) // changed to integer in case of integer func-space
    {
    }

    // use evaluate of default implementation 
    using BaseType :: evaluate;
    using BaseType :: jacobian;

    GeometryType geometryType () const 
    { 
      return storage_.geometryType(); 
    }
    
    inline
    int numBaseFunctions() const;

    inline 
    int numDifferentBaseFunctions() const;

    template <class QuadratureType>
    inline
    void evaluateScalar(int baseFunct, 
			const QuadratureType & quad, int p, 
			ScalarRangeType& phi) const;
    inline
    void evaluateScalar(int baseFunct, 
			const DomainType& xLocal, 
			ScalarRangeType& phi) const;

    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    inline
    void evaluate(int baseFunct, 
                  const FieldVector<deriType, diffOrd> &diffVariable, 
                  QuadratureType & quad, 
                  int quadPoint, RangeType & phi ) const;

    inline
    void jacobianScalar(int baseFunct, const DomainType& xLocal, 
                  ScalarJacobianRangeType& gradPhi) const;

    template <class QuadratureImp>
    inline
    void jacobianScalar(int baseFunct, QuadratureImp& quad, int quadPoint,
                  ScalarJacobianRangeType& gradPhi) const;
    inline
    void jacobian(int baseFunct, const DomainType& xLocal, 
                  JacobianRangeType& gradPhi) const;

    template <class QuadratureImp>
    inline
    void jacobian(int baseFunct, QuadratureImp& quad, int quadPoint,
                  JacobianRangeType& gradPhi) const;

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const DomainType &x,const RangeType &psi) const */
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const DomainType &x,
                                           const RangeType &psi ) const
    {
      ScalarRangeType phi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );

      FieldVector< deriType, 0 > diffVar;
      storage_.evaluate( scalarBaseFunction, diffVar, x, phi );
      return psi[ util_.component( baseFunction ) ] * phi[ 0 ];
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,const RangeType &psi) const */
    template< class QuadratureType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const QuadratureType &quadrature,
                                           const int quadPoint,
                                           const RangeType &psi ) const
    {
      ScalarRangeType phi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );

      FieldVector< deriType, 0 > diffVar;
      storage_.evaluate
        ( scalarBaseFunction, diffVar, quadrature, quadPoint, phi );
      return psi[ util_.component( baseFunction ) ] * phi[ 0 ];
    }
     
    template <class Entity>
    inline
    DofType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;
    
    template <class Entity, class QuadratureType>
    inline
    DofType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const QuadratureType& quad, int quadPoint,
                                   const JacobianRangeType& factor) const;

  private:
    inline VectorialBaseFunctionSet ( const VectorialBaseFunctionSet & );
      
    StorageType storage_;
    DofConversionUtility<PointBased> util_;

    mutable FieldVector<int, 0> diffVar0_;
    mutable FieldVector<int, 1> diffVar1_;
    mutable ScalarRangeType tmp_;
    mutable ScalarJacobianRangeType jTmp_;
  };

  /** \} */

} // end namespace Dune

#include "basefunctionsets.cc"
#endif
