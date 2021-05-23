#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH

#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/fem/misc/threads/threadsafevalue.hh>

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/space/basisfunctionset/transformed.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/combinedspace/interpolation.hh>
#include <dune/fem/space/common/localinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {
      // LocalFunctionWrapper
      // --------------------
      template< class LocalFunction, class BasisFunctionSet >
      struct LocalFunctionWrapper
      {
        struct Traits
        {
          typedef typename LocalFunction::DomainType DomainType;
          typedef typename LocalFunction::RangeType RangeType;
        };
        typedef typename LocalFunction::DomainType DomainType;
        typedef typename LocalFunction::RangeType RangeType;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSet &bset ) : lf_( lf ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          lf_.evaluate( x, y );
        }
        template< class Arg >
        typename Traits::RangeType operator()(const Arg &x) const
        {
          typename Traits::RangeType y;
          evaluate(x,y);
          return y;
        }

      private:
        const LocalFunction &lf_;
      };
      // LocalFunctionWrapper
      // --------------------
      template< class LocalFunction, class Entity, class ShapeFunctionSet, class Transformation >
      struct LocalFunctionWrapper< LocalFunction, TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > >
      {
        typedef TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > BasisFunctionSetType;

        struct Traits
        {
          typedef typename LocalFunction::DomainType DomainType;
          typedef typename LocalFunction::RangeType RangeType;
        };
        typedef typename LocalFunction::DomainType DomainType;
        typedef typename LocalFunction::RangeType RangeType;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSetType &bset )
          : lf_( lf ), geometry_( lf_.entity().geometry() ), bset_( bset ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          typename Traits::RangeType help;
          lf_.evaluate( x, help );
          typename Transformation::InverseTransformationType transf( geometry_, x );
          y = transf.apply( help );
        }
        template< class Arg >
        typename Traits::RangeType operator() ( const Arg &x ) const
        {
          typename Traits::RangeType y;
          evaluate(x,y);
          return y;
        }

      private:
        const LocalFunction &lf_;
        const typename Entity::Geometry geometry_;
        const BasisFunctionSetType &bset_;
      };

    } // namespace Impl


    // LocalFiniteElementInterpolation
    // -------------------------------

    template< class Space, class LocalInterpolation, bool scalarSFS >
    class LocalFiniteElementInterpolation;

    // a vector valued shape basis function set taken from
    // dune-localfunction - the vector value 'blow up' is not yet supported
    // for this case
    template< class Space, class LocalInterpolation >
    class LocalFiniteElementInterpolation<Space,LocalInterpolation,false>
    {
      typedef LocalFiniteElementInterpolation< Space, LocalInterpolation,false > ThisType;

    public:
      typedef typename Space::BasisFunctionSetType BasisFunctionSetType;
      typedef LocalInterpolation LocalInterpolationType;

    private:
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef std::size_t size_type;

      template< class LocalFunction >
      using LocalFunctionWrapper = Impl::LocalFunctionWrapper< LocalFunction, BasisFunctionSetType >;

    public:
      LocalFiniteElementInterpolation()
        : basisFunctionSet_()
        , localInterpolation_()
      {
      }

      explicit LocalFiniteElementInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                 const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
      : basisFunctionSet_( basisFunctionSet ),
        localInterpolation_( localInterpolation )
      {}

      ThisType& operator=( const ThisType& other )
      {
        basisFunctionSet_ = other.basisFunctionSet_;
        localInterpolation_.emplace( other.localInterpolation() );
        return *this;
      }

      void unbind()
      {
        // basisFunctionSet_ = BasisFunctionSetType();
      }

      template< class LocalFunction, class Dof >
      void operator() ( const LocalFunction &localFunction, std::vector< Dof > &localDofVector ) const
      {
        LocalFunctionWrapper< LocalFunction > wrapper( localFunction, basisFunctionSet() );
        localInterpolation().interpolate( wrapper, localDofVector );
      }

      template< class LocalFunction, class Dof, class Allocator >
      void operator() ( const LocalFunction &localFunction, Dune::DynamicVector< Dof, Allocator > &localDofVector ) const
      {
        (*this)(localFunction,localDofVector.container() );
      }

      template< class LocalFunction, class DiscreteFunction, template< class > class Assembly >
      void operator() ( const LocalFunction &localFunction, LocalContribution< DiscreteFunction, Assembly > &localContribution ) const
      {
        (*this)(localFunction,localContribution.localDofVector());
      }

      BasisFunctionSetType basisFunctionSet () const { return basisFunctionSet_; }
      const LocalInterpolationType &localInterpolation () const
      {
        assert( localInterpolation_.has_value() );
        return *localInterpolation_;
      }

    private:
      BasisFunctionSetType    basisFunctionSet_;
      std::optional< LocalInterpolationType > localInterpolation_;
    };

    namespace Impl
    {
      template <int dimRange>
      struct RangeConverter
      {
        RangeConverter ( std::size_t range ) : range_( range ) {}

        template< class T >
        FieldVector< T, 1 > operator() ( const FieldVector< T, dimRange > &in ) const
        {
          return in[ range_ ];
        }

        template< class T, int j >
        FieldMatrix< T, 1, j > operator() ( const FieldMatrix< T, dimRange, j > &in ) const
        {
          // return in[ range_ ]; // implicit conversion fails
          FieldMatrix<T,1,j> ret;
          ret[0] = in[range_];
          return ret;
        }

      protected:
        std::size_t range_;
      };
      template <class DofVector, class DofAlignment>
      struct SubDofVectorWrapper
        : public SubDofVector< DofVector, DofAlignment >
      {
        typedef SubDofVector< DofVector, DofAlignment > BaseType;

        SubDofVectorWrapper( DofVector& dofs, int coordinate, const DofAlignment &dofAlignment )
          : BaseType( dofs, coordinate, dofAlignment )
        {}

        //! do nothing on clear/resize since it's done in apply of this class
        void clear() {}
        void resize( const unsigned int) {}
      };

    }

    // a scalar valued shape basis function set taken from
    // dune-localfunction - for this the vector value 'blow up' is supported
    // as with other spaces
    template< class Space, class LocalInterpolation >
    class LocalFiniteElementInterpolation<Space,LocalInterpolation,true>
    {
      typedef LocalFiniteElementInterpolation< Space, LocalInterpolation,true > ThisType;

    public:
      typedef typename Space::BasisFunctionSetType BasisFunctionSetType;
      typedef LocalInterpolation LocalInterpolationType;

    private:
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;
      // typedef typename Space::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      static const int dimRange = FunctionSpaceType::dimRange;
      static const int dimR = Space::FunctionSpaceType::dimRange;

      typedef std::size_t size_type;

      typedef VerticalDofAlignment< BasisFunctionSetType, RangeType> DofAlignmentType;

    public:
      LocalFiniteElementInterpolation ()
        : basisFunctionSet_(),
          localInterpolation_(),
          dofAlignment_()
      {}

      explicit LocalFiniteElementInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                 const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
      : basisFunctionSet_( basisFunctionSet ),
        localInterpolation_( localInterpolation ),
        dofAlignment_( basisFunctionSet_ )
      {}

      LocalFiniteElementInterpolation ( const ThisType& other )
      : basisFunctionSet_( other.basisFunctionSet_ ),
        localInterpolation_( other.localInterpolation_ ),
        dofAlignment_( basisFunctionSet_ )
      {}

      ThisType& operator=( const ThisType& other )
      {
        basisFunctionSet_ = other.basisFunctionSet_;
        localInterpolation_.emplace( other.localInterpolation() );
        return *this;
      }

      void unbind()
      {
        // basisFunctionSet_ = BasisFunctionSetType();
      }

      template< class LocalFunction, class Vector>
      void operator() ( const LocalFunction &localFunction, Vector &localDofVector ) const
      {
        // clear dofs before something is adedd
        // localDofVector.clear(); // does not exist on DynVector so use 'fill' instead
        std::fill(localDofVector.begin(),localDofVector.end(),0);
        for( std::size_t i = 0; i < dimR; ++i )
        {
          Impl::SubDofVectorWrapper< Vector, DofAlignmentType > subLdv( localDofVector, i, dofAlignment_ );
          (*this)(
              localFunctionConverter( localFunction, Impl::RangeConverter<dimR>(i) ),
              subLdv, PriorityTag<42>()
              );
        }
      }

      template< class LocalFunction, class Dof, class Allocator >
      void operator() ( const LocalFunction &localFunction, Dune::DynamicVector< Dof, Allocator > &localDofVector ) const
      {
        (*this)(localFunction,localDofVector.container() );
      }

      template< class LocalFunction, class DiscreteFunction, template< class > class Assembly >
      void operator() ( const LocalFunction &localFunction, LocalContribution< DiscreteFunction, Assembly > &localContribution ) const
      {
        (*this)(localFunction,localContribution.localDofVector());
      }

      BasisFunctionSetType basisFunctionSet () const { return basisFunctionSet_; }
      const LocalInterpolationType &localInterpolation () const
      {
        assert( localInterpolation_.has_value() );
        return *localInterpolation_;
      }

    private:
      template< class LocalFunction, class Vector>
      auto operator() ( const LocalFunction &localFunction, Vector &localDofVector, PriorityTag<1> ) const
      -> void_t< decltype(
          std::declval<LocalInterpolationType>().interpolate(
            localFunction, localDofVector)) >
      {
        localInterpolation().interpolate( localFunction, localDofVector );
      }
      template< class LocalFunction, class Vector>
      void operator() ( const LocalFunction &localFunction, Vector &localDofVector, PriorityTag<0> ) const
      {
        std::vector<typename Vector::value_type> tmp(basisFunctionSet_.size());
        localInterpolation().interpolate( localFunction, tmp);
        for (unsigned int i=0;i<tmp.size();++i)
          localDofVector[i] = tmp[i];
      }
      BasisFunctionSetType basisFunctionSet_;
      std::optional< LocalInterpolationType > localInterpolation_;
      DofAlignmentType dofAlignment_;
    };

    template <class DiscreteFunctionSpace >
    class LocalFEInterpolationWrapper
      : public LocalInterpolationWrapper< DiscreteFunctionSpace >
    {
      typedef LocalInterpolationWrapper< DiscreteFunctionSpace > BaseType;
    protected:
      using BaseType :: interpolation_;
      typedef std::vector< typename DiscreteFunctionSpace::RangeFieldType >
        TemporarayDofVectorType;
      mutable ThreadSafeValue< TemporarayDofVectorType > tmpDofs_;

      using BaseType::interpolation;
    public:
      LocalFEInterpolationWrapper( const DiscreteFunctionSpace& space )
        : BaseType( space )
      {}

      using BaseType :: bind;
      using BaseType :: unbind;

      template< class LocalFunction, class Dof >
      void operator() ( const LocalFunction &localFunction, std::vector< Dof > &localDofVector ) const
      {
        // interpolation only works for std::vector storage
        interpolation()( localFunction, localDofVector );
      }

      template< class LocalFunction, class Dof, class Allocator >
      void operator() ( const LocalFunction &localFunction, Dune::DynamicVector< Dof, Allocator > &localDofVector ) const
      {
        // DynamicVector internally stores a std::vector
        (*this)(localFunction, localDofVector.container() );
      }


      template< class LocalFunction, class DiscreteFunction, template< class > class Assembly >
      void operator() ( const LocalFunction &localFunction, LocalContribution< DiscreteFunction, Assembly > &localContribution ) const
      {
        // LocalContribution internally stores a std::vector
        (*this)(localFunction, localContribution.localDofVector());
      }

      /** \brief general interpolation, needs to copy from std::vector to dofs */
      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        const int size = dofs.size();
        TemporarayDofVectorType& tmpDofs = *tmpDofs_;
        tmpDofs.resize( size );
        (*this)(localFunction, tmpDofs );

        // copy to dof vector
        for( int i=0; i<size; ++i )
        {
          dofs[ i ] = tmpDofs[ i ];
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH
