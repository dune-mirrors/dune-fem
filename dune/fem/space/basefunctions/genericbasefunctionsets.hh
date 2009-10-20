#ifndef DUNE_GENERICBASEFUNCTIONSETS_HH
#define DUNE_GENERICBASEFUNCTIONSETS_HH

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
  template <class FiniteElementImp, class FunctionSpaceImp>
  class GenericBaseFunctionSet;



  //! Traits class for standard base function set
  template <class FiniteElementImp, class FunctionSpaceImp>
  struct GenericBaseFunctionSetTraits
  {
    //! Export finite element type
    typedef FiniteElementImp                                         FiniteElementType;
    //! Export function space type
    typedef FunctionSpaceImp                                         FunctionSpaceType;
    //! Exact type of the base function
    typedef GenericBaseFunctionSet< FiniteElementType,
                                    FunctionSpaceType >              BaseFunctionSetType;
  };



  /** \class StandardBaseFunctionSet
   *  \brief standard base function set
   */
  template< class FiniteElementImp, class FunctionSpaceImp >
  class GenericBaseFunctionSet
  : public BaseFunctionSetDefault
    < GenericBaseFunctionSetTraits< FiniteElementImp, FunctionSpaceImp > >
  {
  public:
    typedef FiniteElementImp                                         FiniteElementType;
    typedef typename FiniteElementType :: Traits :: LocalBasisType   LocalBasisType;
    typedef FunctionSpaceImp                                         FunctionSpaceType;

    typedef GenericBaseFunctionSetTraits< FiniteElementType,
                                          FunctionSpaceType >        Traits;

  private:
    typedef GenericBaseFunctionSet< FiniteElementType,
                                    FunctionSpaceType >              ThisType;
    typedef BaseFunctionSetDefault< Traits >                         BaseType;
    
  public:
    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;
    static const int dimRange = RangeType :: dimension;
    typedef typename FunctionSpaceType :: RangeFieldType             DofType;
    typedef typename LocalBasisType::Traits::JacobianType JacobianRangeType;
    // typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  public:
    // use evaluate of default implementation 
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! Constructor
    inline explicit GenericBaseFunctionSet ( const FiniteElementType & finiteElement )
    : finiteElement_(finiteElement),
      localBasis_(finiteElement_.localBasis())
    { }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions () const
    {
      return localBasis_.size();
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return finiteElement_.type();
    }
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */ 
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      std :: cerr << "derivates of order >= 2 not supported." << std :: endl;
      abort();
    }

    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, 0 > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      static std::vector<RangeType> tmpphi;
      localBasis_.evaluateFunction( coordinate(x), tmpphi );
      phi = tmpphi[baseFunction];
    }

    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, 1 > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      static std::vector<JacobianRangeType> tmpphi;
      localBasis_.evaluateJacobian( coordinate(x), tmpphi );
      for( int j = 0 ; j < dimRange ; ++j )
      {
        phi[j] = tmpphi[baseFunction][j][diffVariable[0]];
      }
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */ 
    template< class PointType >
    inline void jacobian ( const int baseFunction,
                           const PointType &x,
                           JacobianRangeType &phi ) const
    {
      static std::vector<JacobianRangeType> tmpphi;
      localBasis_.evaluateJacobian( coordinate(x), tmpphi );
      phi = tmpphi[baseFunction];
    }
    
  private:
    GenericBaseFunctionSet( const GenericBaseFunctionSet& );

  protected:
    const FiniteElementType &finiteElement_;

  private:
    const LocalBasisType    &localBasis_;
    
  };


} // end namespace Dune

#endif // #ifndef DUNE_GENERICBASEFUNCTIONSETS_HH
