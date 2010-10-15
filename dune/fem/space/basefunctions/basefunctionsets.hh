#ifndef DUNE_FEM_BASEFUNCTIONSETS_HH
#define DUNE_FEM_BASEFUNCTIONSETS_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

#include <dune/fem/space/basefunctions/vectorialbasefunctionset.hh>

//#define USE_BASEFUNCTIONSET_OPTIMIZED
//#define USE_BASEFUNCTIONSET_CODEGEN
//#define BASEFUNCTIONSET_CODEGEN_GENERATE

#ifdef USE_BASEFUNCTIONSET_CODEGEN 
#define USE_BASEFUNCTIONSET_OPTIMIZED
#endif


#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basefunctions/codegen.hh>
#endif

#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
#include <dune/fem/space/basefunctions/evaluatecaller.hh>
#endif

namespace Dune
{
  /** \addtogroup BaseFunction
   *  \{
   */
  
  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;



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
    
  protected:  
    typedef typename Traits::StorageType StorageType; 
  public:
    enum { dimRange = FunctionSpaceType :: dimRange };  
      
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
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */ 
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, diffOrd > &diffVariable,
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
    StorageType storage_;
  };

  /** \} */

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BASEFUNCTIONSETS_HH
