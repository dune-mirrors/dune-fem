#ifndef DUNE_MAPPING_HH
#define DUNE_MAPPING_HH

//- system includes 
#include <vector>

//- local includes 
#include "vectorspace.hh"

namespace Dune{

//! type of derivative component chooser 
typedef int deriType;

//! type of derivative specializer 
template <int dim>
struct DiffVariable
{
  typedef FieldVector<deriType, dim> Type;
};

/** @defgroup Mapping Mapping
  \ingroup OperatorCommon
  Mappings in Dune always map from one vector space into another vector space.
  Mapping are the base class for Operators and Functions. 
  Operators work on vector spaces containing Functions (i.e. Domain and Range are Functions). In
  contrary Functions work on vector spaces containing real numbers (i.e.
  \f$R^n\f$). For both Mapping the base interface class. Furthermore,
  Mapping provided a machinery to combine different mapping linearly. 

  \remarks 
  The interface for Mappings is defined by the class Mapping.

  @{
 */

/** \brief A mapping from one vector space into another
    This class describes a general mapping from the domain vector space into 
    the range vector space.
    It can also be used to construct linear combinations of mappings.

    This two-sided character has the following consequence: when you address
    an object of type mapping or any of its descendants through a reference
    or pointer of type Mapping, the linear combination defined for that mapping
    is evaluated. On the other hand, if you address through a reference of the
    type of any of its descendants (notably Operator and Function), you'll
    get the functionality specific for that type.
*/
template<typename DFieldType,typename RFieldType, class DType, class RType>
class Mapping //: public Vector < RFieldType > 
{
public:
  
  /** \brief domain vector space (for usage in derived classes) 
      This can either be for example a discrete function (in case of
      operators) or a FieldVector (in case of discrete functions)
  */

  typedef DType DomainType;
  /** \brief range vector space (for usage in derived classes) 
      This can either be for example a discrete function (in case of
      operators) or a FieldVector (in case of discrete functions)
  */
  typedef RType  RangeType;
  
  /** \brief type of field the domain vector space, i.e. double */
  typedef DFieldType DomainFieldType;
  
  /** \brief type of field the range vector space, i.e. double */
  typedef RFieldType RangeFieldType;
  
  //! type of this class 
  typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;

  //! create Mappiung with empty linear combination  
  Mapping( ) 
  {
    lincomb_.push_back( term( *this, 1 ) );
  }
  
  //! delete linear combination if necessary  
  virtual ~Mapping () 
  {
  }

  /** \brief add mapping 
      \param m mapping to add 
      \returns new object mapping 
  */
  virtual MappingType operator + (const MappingType &m) const ;
  
  /** \brief substract mapping 
      \param m mapping to substract  
      \returns new object mapping 
  */
  virtual MappingType operator - (const MappingType &m) const ;
  
  /** \brief scale mapping with factor 
      \param f factor with which mapping is scaled 
      \returns new object mapping 
  */
  virtual MappingType operator * (const RangeFieldType &f) const  ;
  
  /** \brief devide  mapping by factor 
      \param f factor with which mapping is devided 
      \returns new object mapping 
  */
  virtual MappingType operator / (const RangeFieldType &f) const  ;
 
  /** \brief assignment of mapping mapping 
      \param m mapping which is copied  
      \returns reference to mapping  
  */
  virtual MappingType& operator  = (const MappingType &m) ;

  //! apply the whole linear combination which was created with the
  //! operators above, using the apply method of the combined operators  
  void operator() (const DomainType &Arg, RangeType &Dest ) const 
  {
    //Dest.clear();
    
    int count = 0;   
    for ( typename std::vector<term>::const_iterator it = lincomb_.begin(); it != lincomb_.end(); it++ ) 
    {
      if ( count == 0 ) {
        it->v_->apply( Arg, Dest );
        if ( it->scalar_ != 1. ) {
          Dest *= it->scalar_;
        } 
      } else {
        RangeType tmp( Dest );
        it->v_->apply( Arg, tmp );
        if ( it->scalar_ == 1. ) {
          Dest += tmp;
        } else if ( it->scalar_ == -1. ) {
          Dest -= tmp;
        } else {
          tmp *= it->scalar_;
          Dest += tmp;
        }
      }
      count++;
    }
  }
private:
  virtual void apply (const DomainType &Arg, RangeType &Dest) const {
    operator()(Arg, Dest);
  }

  struct term {
    term() : v_(NULL), scalar_(1.0), scaleIt_(false) { }

    term(const MappingType &mapping, RangeFieldType scalar ) : v_(&mapping), scalar_(scalar), scaleIt_( true ) {
      if ( scalar_ == 1. ) {
        scaleIt_ = false;
      }
    }

    const MappingType *v_;
    RangeFieldType scalar_;
    bool scaleIt_;
  };

  std::vector<term> lincomb_;
};

#include "mapping_imp.cc"

/** @} end documentation group */
} // end namespace Dune 
#endif
