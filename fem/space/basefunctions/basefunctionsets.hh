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
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */ 
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */ 
    template< class PointType >
    inline void jacobian ( const int baseFunction,
                           const PointType &x,
                           JacobianRangeType &phi ) const
    {
      storage_.jacobian( baseFunction, x, phi );
    }
    
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
  //! \f$ \Phi_j = \phi_i e_k \f$, where \f$ e_k = [ \delta_{ik} ]_i \f$.
  template< class FunctionSpaceImp, template< class > class StorageImp >
  class VectorialBaseFunctionSet
  : public BaseFunctionSetDefault
    < VectorialBaseFunctionSetTraits< FunctionSpaceImp, StorageImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

    typedef VectorialBaseFunctionSetTraits< FunctionSpaceType, StorageImp > Traits;
    
  private:
    typedef VectorialBaseFunctionSet< FunctionSpaceType, StorageImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;
    
  public:
    typedef typename FunctionSpaceType :: ScalarFunctionSpaceType
      ScalarFunctionSpaceType;
    
  private:  
    typedef BaseFunctionFactory< ScalarFunctionSpaceType > FactoryType;
    typedef typename FactoryType :: BaseFunctionType BaseFunctionType;
    typedef StorageImp< ScalarFunctionSpaceType > StorageType;

  public:
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    enum
    {
      dimDomain = FunctionSpaceType :: dimDomain,
      dimRange = FunctionSpaceType :: dimRange
    };
    
    typedef RangeFieldType DofType;
    
    typedef typename ScalarFunctionSpaceType :: RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType :: JacobianRangeType
      ScalarJacobianRangeType;
    
  protected:
    StorageType storage_;
    DofConversionUtility<PointBased> util_;

  public:
    using BaseType :: evaluate;
    using BaseType :: evaluateSingle;
    using BaseType :: evaluateGradientSingle;
    using BaseType :: jacobian;
       
  public:
    //! Constructor
    inline explicit VectorialBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory ),
      util_( FunctionSpaceType::dimRange )
    {
    }

  private:
    VectorialBaseFunctionSet ( const ThisType & );
    
  public:
    inline GeometryType geometryType () const
    { 
      return storage_.geometryType(); 
    }
    
    inline int numBaseFunctions () const
    {
      return numDifferentBaseFunctions() * dimRange;
    }

    inline int numDifferentBaseFunctions () const
    {
      return storage_.numBaseFunctions();
    }

    template< int diffOrd, class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const FieldVector< int, diffOrd > &diffVariable,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }
    
    template< class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      FieldVector< int, 0 > diffVar;
      evaluateScalar( baseFunction, diffVar, x, phi );
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      ScalarRangeType tmp;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      evaluateScalar( scalarBaseFunction, diffVariable, x, tmp );

      phi = 0;
      phi[ util_.component( baseFunction ) ] = tmp[ 0 ];
    }

    template< class PointType >
    inline void jacobianScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarJacobianRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.jacobian( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */
    template< class PointType >
    inline void jacobian ( const int baseFunction,
                           const PointType &x, 
                           JacobianRangeType &phi ) const
    {
      ScalarJacobianRangeType tmp;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      jacobianScalar( scalarBaseFunction, x, tmp );

      phi = 0;
      phi[ util_.component( baseFunction )] = tmp[ 0 ];
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           const RangeType &psi ) const
    {
      ScalarRangeType phi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );

      evaluateScalar( scalarBaseFunction, x, phi );
      return psi[ util_.component( baseFunction ) ] * phi[ 0 ];
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const PointType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    inline RangeFieldType evaluateGradientSingle( const int baseFunction,
                                                  const EntityType &entity,
                                                  const PointType &x,
                                                  const JacobianRangeType &psi ) const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;

      const GeometryType &geometry = entity.geometry();
      const GeometryJacobianType &jacobianInverseTransposed
        = geometry.jacobianInverseTransposed( coordinate( x ) );
      
      ScalarJacobianRangeType gradPhi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      jacobianScalar( scalarBaseFunction, x, gradPhi );
    
      DomainType gradScaled( 0 );
      jacobianInverseTransposed.umv( gradPhi[ 0 ], gradScaled );
      return gradScaled * psi[ util_.component( baseFunction ) ];
    }
  };

  /** \} */

} // end namespace Dune

#endif
