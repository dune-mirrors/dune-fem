#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PIOLATRANSFORMATION_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_PIOLATRANSFORMATION_HH

#include <algorithm>

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/diagonalmatrix.hh>

#include <dune/fem/common/fmatrixcol.hh>

namespace Dune
{

  namespace Fem
  {


    // usefull implementations
    // -----------------------

    template< class Mat >
    double determinante ( const Mat &m )
    { return m.det(); }

    template< class F, int d, int l >
    double determinante ( const FieldMatrix< F, d, l > &m )
    { return m.determinant(); }

    template< class F, int d >
    double determinante ( const DiagonalMatrix< F, d > &m )
    { return m.determinant(); }


    // internal forward declaration
    // ----------------------------

    template< class Geometry, int dimRange >
    class InversePiolaTransformation;


    // PiolaTransformation
    // -------------------

    template< class Geometry, int dimRange >
    class PiolaTransformation
    {
      typedef PiolaTransformation< Geometry, dimRange > ThisType;

      static const int dimDomain = Geometry::GlobalCoordinate::dimension;
      typedef typename Geometry::JacobianTransposed JacobianTransposed;

      static_assert( dimRange % dimDomain == 0, "PiolaTransformation can only be applied if dimRange is a multiple of dimDomain." );

      typedef typename FieldTraits< JacobianTransposed >::real_type real_type;
      static const int blocks = dimRange / dimDomain;

    public:
      typedef InversePiolaTransformation< Geometry, dimRange > InverseTransformationType;

      template< class Point >
      PiolaTransformation ( const Geometry &geo, const Point &p )
        : gjt_( geo.jacobianTransposed( p ) ),
        detInv_( 1.0 / determinante( gjt_ ) )
      {}

      template< class F >
      FieldVector< F, dimRange > apply ( const FieldVector< F, dimRange > &d ) const
      {
        FieldVector< F, dimRange > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < blocks; ++r )
        {
          std::copy_n( d.begin() + r*dimDomain, dimDomain, arg.begin() );
          gjt_.mtv( arg, dest );
          std::copy_n( dest.begin(), dimDomain, ret.begin() + r*dimDomain );
        }
        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldVector< F, dimRange > apply_t ( const FieldVector< F, dimRange > &d ) const
      {
        FieldVector< F, dimRange > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < blocks; ++r )
        {
          std::copy_n( ret.begin() + r*dimDomain, dimDomain, arg.begin() );
          gjt_.mv( arg, dest );
          std::copy_n( dest.begin(), dimDomain, ret.begin() + r*dimDomain );
        }
        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldMatrix< F, dimRange, dimDomain > apply ( const FieldMatrix< F, dimRange, dimDomain > &d ) const
      {
        FieldMatrix< F, dimRange, dimDomain > ret( d );
        FieldVector< F, dimDomain > arg, dest;

        for( std::size_t r = 0; r < dimDomain; ++r )
        {
          FieldMatrixColumn< FieldMatrix< F, dimRange, dimDomain > > col( ret, r );
          for( std::size_t b = 0; b < blocks; ++b )
          {
            std::copy_n( col.begin() + b*dimDomain, dimDomain, arg.begin() );
            gjt_.mv( arg, dest );
            std::copy_n( dest.begin(), dimDomain, col.begin() + b*dimDomain );
          }
        }

        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldMatrix< F, dimRange, dimDomain > apply_t ( const FieldMatrix< F, dimRange, dimDomain > &d ) const
      {
        FieldMatrix< F, dimRange, dimDomain > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < dimDomain; ++r )
        {
          FieldMatrixColumn< FieldMatrix< F, dimRange, dimDomain > > col( ret, r );
          for( std::size_t b = 0; b < blocks; ++b )
          {
            std::copy_n( col.begin() + b*dimDomain, dimDomain, arg.begin() );
            gjt_.mtv( arg, dest );
            std::copy_n( dest.begin(), dimDomain, col.begin() + b*dimDomain );
          }
        }
        ret *= detInv_;
        return ret;
      }

    protected:
      JacobianTransposed gjt_;
      real_type detInv_;
    };



    // InverseTransformationType
    // -------------------------

    template< class Geometry, int dimRange >
    class InversePiolaTransformation
    {
      typedef InversePiolaTransformation< Geometry, dimRange > ThisType;

      static const int dimDomain = Geometry::GlobalCoordinate::dimension;
      typedef typename Geometry::JacobianInverseTransposed JacobianInverseTransposed;

      static_assert( dimRange % dimDomain == 0, "PiolaTransformation can only be applied if dimRange is a multiple of dimDomain." );

      typedef typename FieldTraits< JacobianInverseTransposed >::real_type real_type;
      static const int blocks = dimRange / dimDomain;

    public:
      typedef PiolaTransformation< Geometry, dimRange > InverseTransformationType;

      template< class Point >
      InversePiolaTransformation ( const Geometry &geo, const Point &p )
        : gjit_( geo.jacobianInverseTransposed( p ) ),
        detInv_( 1.0 / determinante( gjit_ ) )
      {}

      template< class F >
      FieldVector< F, dimRange > apply ( const FieldVector< F, dimRange > &d ) const
      {
        FieldVector< F, dimRange > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < blocks; ++r )
        {
          std::copy_n( d.begin() + r*dimDomain, dimDomain, arg.begin() );
          gjit_.mtv( arg, dest );
          std::copy_n( dest.begin(), dimDomain, ret.begin() + r*dimDomain );
        }
        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldVector< F, dimRange > apply_t ( const FieldVector< F, dimRange > &d ) const
      {
        FieldVector< F, dimRange > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < blocks; ++r )
        {
          std::copy_n( ret.begin() + r*dimDomain, dimDomain, arg.begin() );
          gjit_.mv( arg, dest );
          std::copy_n( dest.begin(), dimDomain, ret.begin() + r*dimDomain );
        }
        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldMatrix< F, dimRange, dimDomain > apply ( const FieldMatrix< F, dimRange, dimDomain > &d ) const
      {
        FieldMatrix< F, dimRange, dimDomain > ret( d );
        FieldVector< F, dimDomain > arg, dest;

        for( std::size_t r = 0; r < dimDomain; ++r )
        {
          FieldMatrixColumn< FieldMatrix< F, dimRange, dimDomain > > col( ret, r );
          for( std::size_t b = 0; b < blocks; ++b )
          {
            std::copy_n( col.begin() + b*dimDomain, dimDomain, arg.begin() );
            gjit_.mv( arg, dest );
            std::copy_n( dest.begin(), dimDomain, col.begin() + b*dimDomain );
          }
        }

        ret *= detInv_;
        return ret;
      }

      template< class F >
      FieldMatrix< F, dimRange, dimDomain > apply_t ( const FieldMatrix< F, dimRange, dimDomain > &d ) const
      {
        FieldMatrix< F, dimRange, dimDomain > ret( d );
        FieldVector< F, dimDomain > arg, dest;
        for( std::size_t r = 0; r < dimDomain; ++r )
        {
          FieldMatrixColumn< FieldMatrix< F, dimRange, dimDomain > > col( ret, r );
          for( std::size_t b = 0; b < blocks; ++b )
          {
            std::copy_n( col.begin() + b*dimDomain, dimDomain, arg.begin() );
            gjit_.mtv( arg, dest );
            std::copy_n( dest.begin(), dimDomain, col.begin() + b*dimDomain );
          }
        }
        ret *= detInv_;
        return ret;
      }

    protected:
      JacobianInverseTransposed gjit_;
      real_type detInv_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PIOLATRANSFORMATION_HH
