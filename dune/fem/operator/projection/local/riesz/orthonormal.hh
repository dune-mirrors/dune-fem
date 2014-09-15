#ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_ORTHONORMAL_HH
#define DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_ORTHONORMAL_HH

#include <cassert>

#include <algorithm>
#include <utility>

#include <dune/geometry/referenceelements.hh>

#include "localrieszprojection.hh"

namespace Dune
{

  namespace Fem
  {

    // OrthonormalLocalRieszProjection
    // -------------------------------

    template< class BasisFunctionSet >
    class OrthonormalLocalRieszProjection
    : public LocalRieszProjection< BasisFunctionSet, OrthonormalLocalRieszProjection< BasisFunctionSet > >
    {
      typedef OrthonormalLocalRieszProjection< BasisFunctionSet > ThisType;
      typedef LocalRieszProjection< BasisFunctionSet, OrthonormalLocalRieszProjection< BasisFunctionSet > > BaseType;

    public:
      /** \copydoc Dune::Fem::LocalRieszProjection::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    private:
      typedef typename BasisFunctionSet::FunctionSpaceType::RangeFieldType RangeFieldType;

    public:
      /** \name Construction
       *  \{
       */

      explicit OrthonormalLocalRieszProjection ( const BasisFunctionSetType &basisFunctionSet )
        : basisFunctionSet_( std::forward< BasisFunctionSetType >( basisFunctionSet ) ),
          factor_( ratio( basisFunctionSet.entity().geometry() ) )
      {}

      explicit OrthonormalLocalRieszProjection ( BasisFunctionSetType &&basisFunctionSet )
        : basisFunctionSet_( std::forward< BasisFunctionSetType >( basisFunctionSet ) ),
          factor_( ratio( basisFunctionSet.entity().geometry() ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      OrthonormalLocalRieszProjection ( const ThisType & ) = default;

      OrthonormalLocalRieszProjection ( ThisType &&other )
        : basisFunctionSet_( std::move( other.basisFunctionSet_ ) ),
          factor_( other.factor_ )
      {}

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType &&other )
      {
        basisFunctionSet_ = std::move( other.basisFunctionSet_ );
        factor_( other.factor_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::LocalRieszProjection::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet () const
      {
        return basisFunctionSet_;
      }

      /** \copydoc Dune::Fem::LocalRieszProjection::apply */
      template< class F, class LocalDofVector >
      void apply ( const F &f, LocalDofVector &dofs ) const
      {
        assert( f.size() == dofs.size() );
        std::transform( f.begin(), f.end(), dofs.begin(), [this]( RangeFieldType value ){ return factor_*value; } );
      }

      /** \} */

    private:
      template< class Geometry >
      static RangeFieldType ratio ( const Geometry &geometry )
      {
        assert( geometry.affine() );
        const auto &referenceElement
          = Dune::ReferenceElements< typename Geometry::ctype, Geometry::mydimension >::general( geometry.type() );
        return referenceElement.volume()/geometry.volume();
      }

      BasisFunctionSetType basisFunctionSet_;
      RangeFieldType factor_;
    };

  } // namespace Fem

} // namepsace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_ORTHONORMAL_HH
