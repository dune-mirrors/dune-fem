#ifndef DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH
#define DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH

#include <cstddef>
#include <utility>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/misc/engineconcept.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    // DofAlignment
    // ------------

    /** \brief Interface documentation for Dof alignment classes
     *         used in VectorialBasisFunctionSet.
     *
     *  \note This interface class is implemented by 
     *        HorizontalDofAlignment and 
     *        VerticalDofAlignment
     *
     *  \tparam  Imp  implementation type (CRTP)
     */
    template< class Imp >
    class DofAlignment
    : public EngineDefault< Imp >
    {
      typedef EngineDefault< Imp > BaseType;

    protected:
      using BaseType::asImp;

    public:
      //! \brief global Dof type
      typedef std::size_t GlobalDofType;
      /** \brief local Dof type consists of coordinate number and Dof 
       *         number in scalar basis function set
       */
      typedef std::pair< int, std::size_t > LocalDofType;

      /** \brief map local to global Dof
       *
       *  \note methods localDof and globalDof must be inverse
       *
       *  \param[in]  localDof  local Dof
       *
       *  \returns global Dof
       */
      GlobalDofType globalDof ( const LocalDofType &localDof ) const
      {
        return asImp().globalDof( localDof );
      }

      /** \brief map global to local Dof
       *
       *  \note methods localDof and globalDof must be inverse
       *
       *  \param[in]  globalDof  global Dof
       *
       *  \returns local Dof 
       */
      LocalDofType localDof ( const GlobalDofType &globalDof ) const
      {
        return asImp().globalDof( localDof );
      }
    };



    // HorizontalDofAlignment
    // ----------------------

    /** \brief Implementation of DofAlignment.
     *
     *  This Dof alignment uses consecutive enumeration for
     *  each coordinate in vectial basis function set.
     *
     *  \tparam  ScalarBasisFunctionSet  scalar basis function set
     *  \tparam  Range                   range
     */
    template< class ScalarBasisFunctionSet, class Range >
    class HorizontalDofAlignment
    : public DofAlignment< HorizontalDofAlignment< ScalarBasisFunctionSet, Range > >
    {
      typedef HorizontalDofAlignment< ScalarBasisFunctionSet, Range > ThisType;
      typedef DofAlignment< HorizontalDofAlignment< ScalarBasisFunctionSet, Range > > BaseType;

    public:
      typedef typename BaseType::GlobalDofType GlobalDofType;
      typedef typename BaseType::LocalDofType LocalDofType;

      HorizontalDofAlignment () {}

      explicit HorizontalDofAlignment ( const ScalarBasisFunctionSet &scalarBasisFunctionSet )
      : scalarSize_( scalarBasisFunctionSet.size() )
      {}

      // need to implement copy constructor because of basis class
      HorizontalDofAlignment ( const ThisType &other ) : scalarSize_( other.scalarSize_ ) {}

      // need to implement assignment operator because of basis class
      ThisType &operator= ( const ThisType &other )
      {
        scalarSize_ = other.scalarSize_;
        return *this;
      }

      /** @copydoc Dune::Fem::DofAlignment::globalDof */
      GlobalDofType globalDof ( const LocalDofType &localDof ) const
      {
        return GlobalDofType( localDof.first*scalarSize_ + localDof.second );
      }

      /** @copydoc Dune::Fem::DofAlignment::localDof */
      LocalDofType localDof ( const GlobalDofType &globalDof ) const
      {
         return LocalDofType( globalDof / scalarSize_, globalDof % scalarSize_ ); 
      }

    private:
      std::size_t scalarSize_;
    };



    // VerticalDofAlignment
    // --------------------

    /** \brief Implementation of DofAlignment.
     *
     *  This Dof uses the same Dof enumeration as VectorialShapeFunctionSet.
     *
     *  \tparam  ScalarBasisFunctionSet  scalar basis function set
     *  \tparam  Range                   range
     */
    template< class ScalarBasisFunctionSet, class Range >
    class VerticalDofAlignment
    : public DofAlignment< VerticalDofAlignment< ScalarBasisFunctionSet, Range > >
    {
      typedef VerticalDofAlignment< ScalarBasisFunctionSet, Range > ThisType;
      typedef DofAlignment< VerticalDofAlignment< ScalarBasisFunctionSet, Range > > BaseType;

      static const int dimRange = Range::dimension;

    public:
      typedef typename BaseType::GlobalDofType GlobalDofType;
      typedef typename BaseType::LocalDofType LocalDofType;

      VerticalDofAlignment () {}

      explicit VerticalDofAlignment ( const ScalarBasisFunctionSet & ) {}

      // need to implement copy constructor because of basis class
      VerticalDofAlignment ( const ThisType & ) {}

      // need to implement assignment operator because of basis class
      ThisType &operator= ( const ThisType & ) { return *this; }

      /** @copydoc Dune::Fem::DofAlignment::globalDof */
      GlobalDofType globalDof ( const LocalDofType &localDof ) const
      {
        return GlobalDofType( localDof.first + localDof.second*dimRange );
      }

      /** @copydoc Dune::Fem::DofAlignment::localDof */
      LocalDofType localDof ( const GlobalDofType &globalDof ) const
      {
        return LocalDofType( globalDof % dimRange, globalDof / dimRange ); 
      }
    };



    // SubDofVector
    // ------------

    /** \brief Extract Sub dof vector for single coordinate.
     *
     *  \tparam  DofVector     Dof vector
     *  \tparam  DofAlignment  implementation of DofAlignment
     */
    template< class DofVector, class DofAlignment >
    class SubDofVector;

    template< class DofVector, class ScalarBasisFunctionSet, class Range >
    class SubDofVector< DofVector, HorizontalDofAlignment< ScalarBasisFunctionSet, Range > >
    {
      typedef SubDofVector< DofVector, HorizontalDofAlignment< ScalarBasisFunctionSet, Range > > ThisType;

      typedef typename Range::value_type RangeFieldType;

      typedef HorizontalDofAlignment< ScalarBasisFunctionSet, Range > DofAlignmentType;
      typedef typename DofAlignmentType::LocalDofType LocalDofType;

    public:
      typedef RangeFieldType value_type;

      SubDofVector( const DofVector &dofs, int coordinate, const DofAlignmentType &dofAlignment )
      : dofs_( &(dofs[ dofAlignment.globalDof( LocalDofType( coordinate, 0 ) ) ] ) )
      {}

      const RangeFieldType &operator[] ( std::size_t i ) const
      {
        return dofs_[ i ];
      }

    private:
      const RangeFieldType *dofs_;
    };

    template< class DofVector, class ScalarBasisFunctionSet, class Range >
    class SubDofVector< DofVector, VerticalDofAlignment< ScalarBasisFunctionSet, Range > >
    {
      typedef SubDofVector< DofVector, VerticalDofAlignment< ScalarBasisFunctionSet, Range > > ThisType;

      typedef typename Range::value_type RangeFieldType;

      typedef VerticalDofAlignment< ScalarBasisFunctionSet, Range > DofAlignmentType;
      typedef typename DofAlignmentType::LocalDofType LocalDofType;

    public:
      typedef RangeFieldType value_type;

      SubDofVector( const DofVector &dofs, int coordinate, const DofAlignmentType &dofAlignment )
      : dofs_( dofs ),
        coordinate_( coordinate ),
        dofAlignment_( dofAlignment )
      {}

      const RangeFieldType &operator[] ( std::size_t i ) const
      {
        return dofs_[ dofAlignment_.globalDof( LocalDofType( coordinate_, i ) ) ];
      }

    private:
      const DofVector &dofs_;
      int coordinate_;
      DofAlignmentType dofAlignment_;
    };



    // VectorialBasisFunctionSet
    // -------------------------

    /** \brief Builds a vectorial basis function set 
     *         from given scalar basis function set.
     *
     *  \tparam  ScalarBasisFunctionSet  scalar basis function set
     *  \tparam  Range                   range vector
     *  \tparam  DofAlignment            type of dof alignment
     */
    template< class ScalarBasisFunctionSet, class Range, template< class, class > class DofAlignment = VerticalDofAlignment >
    class VectorialBasisFunctionSet
    {
      typedef VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range, DofAlignment > ThisType;

    public:
      typedef ScalarBasisFunctionSet ScalarBasisFunctionSetType;

      typedef typename ScalarBasisFunctionSetType::EntityType EntityType;
      typedef typename ScalarBasisFunctionSetType::ReferenceElementType ReferenceElementType;

    private:
      typedef typename ScalarBasisFunctionSetType::FunctionSpaceType ScalarFunctionSpaceType;
      static const int dimRange = Range::dimension;

    public:
      typedef typename ToNewDimRangeFunctionSpace< ScalarFunctionSpaceType, dimRange >::Type FunctionSpaceType;
       
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef DofAlignment< ScalarBasisFunctionSet, Range > DofAlignmentType;

    private:
      typedef typename DofAlignmentType::GlobalDofType GlobalDofType;
      typedef typename DofAlignmentType::LocalDofType LocalDofType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      struct EvaluateAll;
      struct JacobianAll;
      struct HessianAll;

    public:
      VectorialBasisFunctionSet () {}

      explicit VectorialBasisFunctionSet ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet )
      : scalarBasisFunctionSet_( scalarBasisFunctionSet ),
        dofAlignment_( scalarBasisFunctionSet_ )
      {}

      int order () const { return scalarBasisFunctionSet().order(); }

      std::size_t size () const { return dimRange*scalarBasisFunctionSet().size(); }

      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return scalarBasisFunctionSet().type(); }

      const ReferenceElementType &referenceElement () const { return scalarBasisFunctionSet().referenceElement(); }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        axpy< EvaluateAll >( x, valueFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        axpy< JacobianAll >( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        evaluateAll< EvaluateAll >( x, dofs, value );
      }

      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        evaluateAll< EvaluateAll >( x, values );
      }

      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        evaluateAll< JacobianAll >( x, dofs, jacobian );
      }

      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        evaluateAll< JacobianAll >( x, jacobians );
      }

      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        evaluateAll< HessianAll >( x, dofs, hessian );
      }

      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        evaluateAll< HessianAll >( x, hessians );
      }

      const EntityType &entity () const { return scalarBasisFunctionSet().entity(); }

      DofAlignmentType dofAlignment () const { return dofAlignment_ ; }

      const ScalarBasisFunctionSetType &scalarBasisFunctionSet () const
      {
        return scalarBasisFunctionSet_;
      }

    private:
      template< class Evaluate, class Point, class DofVector >
      void axpy ( const Point &x, const typename Evaluate::Vector &factor, DofVector &dofs ) const
      {
        const std::size_t size = scalarBasisFunctionSet().size();
        std::vector< typename Evaluate::Scalar > scalars( size );
        Evaluate::apply( scalarBasisFunctionSet(), x, scalars );

        for( int r = 0; r < dimRange; ++r )
        {
          for( std::size_t i = 0; i < size; ++i )
          {
            const GlobalDofType globalDof = dofAlignment_.globalDof( LocalDofType( r, i ) );
            dofs[ globalDof ] += factor[ r ] * scalars[ i ][ 0 ];
          }
        }
      }

      template< class Evaluate, class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, typename Evaluate::Vector &vector ) const
      {
        typename Evaluate::Scalar scalar;
        for( int r = 0; r < dimRange; ++r )
        {
          SubDofVector< DofVector, DofAlignmentType > subDofs( dofs, r, dofAlignment_ );
          Evaluate::apply( scalarBasisFunctionSet(), x, subDofs, scalar );
          vector[ r ] = scalar[ 0 ];
        }
      }

      template< class Evaluate, class Point, class VectorArray >
      void evaluateAll ( const Point &x, VectorArray &vectorials ) const
      {
        const std::size_t size = scalarBasisFunctionSet().size();
        std::vector< typename Evaluate::Scalar > scalars( size );
        Evaluate::apply( scalarBasisFunctionSet(), x, scalars );

        typedef typename Evaluate::Vector Vector;

        for( int r = 0; r < dimRange; ++r )
        {
          for( std::size_t i = 0; i < size; ++i )
          {
            const GlobalDofType globalDof = dofAlignment_.globalDof( LocalDofType( r, i ) );
            Vector &vector = vectorials[ globalDof ];
            vector = Vector( typename Vector::value_type( 0 ) );
            vector[ r ] = scalars[ i ][ 0 ];
          }
        }
      }

      ScalarBasisFunctionSetType scalarBasisFunctionSet_;
      DofAlignmentType dofAlignment_;
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::EvaluateAll
    // -----------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range, template< class, class > class DofAlignment >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range, DofAlignment >::EvaluateAll
    {
      typedef typename ScalarFunctionSpaceType::RangeType Scalar;
      typedef RangeType Vector;

      template< class Point, class SubDofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const SubDofVector &dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.evaluateAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.evaluateAll( x, scalars );
      }
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::JacobianAll
    // -----------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range, template< class, class > class DofAlignment >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range, DofAlignment >::JacobianAll
    {
      typedef typename ScalarFunctionSpaceType::JacobianRangeType Scalar;
      typedef JacobianRangeType Vector;

      template< class Point, class SubDofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const SubDofVector &dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.jacobianAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.jacobianAll( x, scalars );
      }
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::HessianAll
    // ----------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range, template< class, class > class DofAlignment >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range, DofAlignment >::HessianAll
    {
      typedef typename ScalarFunctionSpaceType::HessianRangeType Scalar;
      typedef HessianRangeType Vector;

      template< class Point, class SubDofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const SubDofVector &dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.hessianAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.hessianAll( x, scalars );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH
