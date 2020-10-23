#ifndef DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH
#define DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class >
    struct DiscreteFunctionTraits;

    class IsDiscreteFunction;



    namespace Assembly
    {

      namespace Global
      {

        // AddBase
        // -------

        template< class DiscreteFunction >
        struct AddBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        {
          typedef typename DiscreteFunction::DofType DofType;

          static void begin ( DiscreteFunction &df )
          {
            typedef typename DiscreteFunction::DiscreteFunctionSpaceType::LocalBlockIndices LocalBlockIndices;

            // clear slave DoFs
            auto &dofVector = df.dofVector();
            for( const auto &slaveDof : df.space().slaveDofs() )
              Hybrid::forEach( LocalBlockIndices(), [ &dofVector, &slaveDof ] ( auto &&j ) { dofVector[ slaveDof ][ j ] = DofType( 0 ); } );
          }

          static void end ( DiscreteFunction &df ) { df.space().communicate( df, DFCommunicationOperation::Add() ); }
        };



        // SetBase
        // -------

        template< class DiscreteFunction >
        struct SetBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        {
          static void begin ( DiscreteFunction &df ) {}
          static void end ( DiscreteFunction &df ) { df.space().communicate( df, DFCommunicationOperation::Copy() ); }
        };

      } // namespace Global



      // AddBase
      // -------

      template< class DiscreteFunction >
      struct AddBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      {
        typedef typename DiscreteFunction::DofType DofType;

        typedef Global::Add< DiscreteFunction > GlobalOperationType;

        template< class Entity, class LocalDofVector >
        void begin ( const Entity &entity, const DiscreteFunction &df, LocalDofVector &localDofVector ) const
        {
          std::fill( localDofVector.begin(), localDofVector.end(), DofType( 0 ) );
        }

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.addLocalDofs( entity, localDofVector );
        }
      };



      // AddScaledBase
      // -------------

      template< class DiscreteFunction >
      struct AddScaledBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        : public AddBase< DiscreteFunction >
      {
        AddScaledBase ( typename DiscreteFunction::DofType factor ) : factor_( std::move( factor ) ) {}

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.addScaledLocalDofs( entity, factor_, localDofVector );
        }

      private:
        typename DiscreteFunction::DofType factor_;
      };



      // SetBase
      // -------

      template< class DiscreteFunction >
      struct SetBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      {
        typedef typename DiscreteFunction::DofType DofType;

        typedef Global::Set< DiscreteFunction > GlobalOperationType;

        template< class Entity, class LocalDofVector >
        void begin ( const Entity &entity, const DiscreteFunction &df, LocalDofVector &localDofVector ) const
        {
          std::fill( localDofVector.begin(), localDofVector.end(), DofType( 0 ) );
        }

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.setLocalDofs( entity, localDofVector );
        }
      };


      // SetSelectedBase
      // ---------------

      template< class DiscreteFunction >
      struct SetSelectedBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      {
        typedef typename DiscreteFunction::DofType DofType;

        typedef Global::Set< DiscreteFunction > GlobalOperationType;

        template< class Entity, class LocalDofVector >
        void begin ( const Entity &entity, const DiscreteFunction &df, LocalDofVector &localDofVector ) const
        {
          df.getLocalDofs ( entity, localDofVector );
        }

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.setLocalDofs( entity, localDofVector );
        }
      };

    } // namespace Assembly

  } // namespace Fem



  // DenseMatVecTraits for LocalContribution for Discrete Functions
  // --------------------------------------------------------------

  template< class DiscreteFunction, template< class > class AssemblyOperation >
  struct DenseMatVecTraits< Fem::LocalContribution< DiscreteFunction, AssemblyOperation, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > > >
  {
    typedef Fem::LocalContribution< DiscreteFunction, AssemblyOperation > derived_type;
    typedef typename DiscreteFunction::DofType value_type;
    typedef std::vector< value_type > container_type;
    typedef typename container_type::size_type size_type;
  };



  // FieldTraits for LocalContribution for Discrete Functions
  // --------------------------------------------------------

  template< class DiscreteFunction, template< class > class AssemblyOperation >
  struct FieldTraits< Fem::LocalContribution< DiscreteFunction, AssemblyOperation, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > > >
  {
    typedef typename FieldTraits< typename DiscreteFunction::DofType >::field_type field_type;
    typedef typename FieldTraits< typename DiscreteFunction::DofType >::real_type real_type;
  };



  namespace Fem
  {

    // LocalContribution for Discrete Functions
    // ----------------------------------------

    template< class DiscreteFunction, template< class > class AssemblyOperation >
    class LocalContribution< DiscreteFunction, AssemblyOperation, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      : public Dune::DenseVector< LocalContribution< DiscreteFunction, AssemblyOperation > >
    {
      typedef LocalContribution< DiscreteFunction, AssemblyOperation > ThisType;
      typedef Dune::DenseVector< LocalContribution< DiscreteFunction, AssemblyOperation > > BaseType;

    public:
      typedef DiscreteFunction DiscreteFunctionType;
      typedef AssemblyOperation< typename DiscreteFunctionTraits< DiscreteFunctionType >::DiscreteFunctionType > AssemblyOperationType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename DiscreteFunctionType::DofType DofType;

      typedef typename DiscreteFunctionType::RangeType RangeType;
      typedef typename RangeType::field_type RangeFieldType;
      typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;

      typedef std::vector< DofType > LocalDofVectorType;
      typedef typename LocalDofVectorType::size_type SizeType;

      typedef typename BasisFunctionSetType::EntityType EntityType;

      template< class... Args >
      explicit LocalContribution ( DiscreteFunctionType &discreteFunction, Args &&... args )
        : discreteFunction_( discreteFunction ),
          localDofVector_( discreteFunction.space().maxNumDofs() ),
          assemblyOperation_( std::forward< Args >( args )... ),
          bound_(false)
      {
        discreteFunction.template beginAssemble< typename AssemblyOperationType::GlobalOperationType >();
      }

      LocalContribution ( const ThisType & ) = delete;
      LocalContribution ( ThisType && ) = delete;

      ~LocalContribution () { discreteFunction().template endAssemble< typename AssemblyOperationType::GlobalOperationType >(); }

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      const DiscreteFunctionType &discreteFunction () const { return discreteFunction_; }
      DiscreteFunctionType &discreteFunction () { return discreteFunction_; }

      /**
       * \brief access to local dofs (read-only)
       *
       * \param[in]  i  local dof number
       * \return reference to dof
       **/
      const DofType &operator[] ( SizeType i ) const { return localDofVector()[ i ]; }

      /**
       * \brief access to local dofs (read-write)
       *
       * \param[in]  i  local DoF number
       * \return reference to DoF
       **/
      DofType &operator[] ( SizeType i ) { return localDofVector()[ i ]; }

      /** \brief set all DoFs to zero **/
      void clear ()
      {
        std::fill( localDofVector().begin(), localDofVector().end(), DofType( 0 ) );
      }

      /**
       * \brief axpy operation for local contribution
       *
       * Denoting the DoFs of the local contribution by \f$u_i\f$ and the basis
       * functions by \f$\varphi_i\f$, this function performs the following
       * operation:
       * \f[
       * u_i = u_i + factor \cdot \varphi_i( x )
       * \f]
       *
       * \param[in]  x       point to evaluate basis functions in
       * \param[in]  factor  axpy factor
       **/
      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor )
      {
        basisFunctionSet().axpy( x, factor, localDofVector() );
      }

      /**
       * \brief axpy operation for local contribution
       *
       * Denoting the DoFs of the local contribution by \f$u_i\f$ and the basis
       * functions by \f$\varphi_i\f$, this function performs the following
       * operation:
       * \f[
       * u_i = u_i + factor \cdot \nabla\varphi_i( x )
       * \f]
       *
       * \param[in]  x       point to evaluate jacobian of basis functions in
       * \param[in]  factor  axpy factor
       **/
      template< class PointType >
      void axpy ( const PointType &x, const JacobianRangeType &factor )
      {
        basisFunctionSet().axpy( x, factor, localDofVector() );
      }

      /**
       * \brief axpy operation for local contribution
       *
       * Denoting the DoFs of the local contribution by \f$u_i\f$ and the basis
       * functions by \f$\varphi_i\f$, this function performs the following
       * operation:
       * \f[
       * u_i = u_i + factor1 \cdot \varphi_i( x ) + factor2 \cdot \nabla\varphi_i( x )
       * \f]
       *
       * \param[in]  x        point to evaluate basis functions in
       * \param[in]  factor1  axpy factor for \f$\varphi( x )\f$
       * \param[in]  factor2  axpy factor for \f$\nabla\varphi( x )\f$
       **/
      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
      {
        basisFunctionSet().axpy( x, factor1, factor2, localDofVector() );
      }

      /**
       * \brief obtain the order of this local contribution
       *
       * The order of a local contribution refers to the polynomial order required
       * to integrate it exactly.
       *
       * \note It is not completely clear what this value should be, e.g., for
       *       bilinear basis functions.
       *
       * \returns order of the local contribution
       **/
      int order () const { return basisFunctionSet().order(); }

      /**
       * \brief obtain the basis function set for this local contribution
       * \returns reference to the basis function set
       */
      const BasisFunctionSetType &basisFunctionSet () const { return basisFunctionSet_; }

      /**
       * \brief obtain the entity, this local contribution lives on
       * \returns reference to the entity
       **/
      const EntityType &entity () const { return basisFunctionSet().entity(); }

      SizeType size () const { return localDofVector().size(); }

      /**
       * \brief evaluate all basisfunctions for all quadrature points, multiply with the given factor and
       *        add the result to the local coefficients
       **/
      template< class QuadratureType, class VectorType >
      void axpyQuadrature ( const QuadratureType &quad, const VectorType &values )
      {
        basisFunctionSet().axpy( quad, values, localDofVector() );
      }

      /**
       * \brief evaluate all basisfunctions for all quadrature points, multiply with the given factor and
       *        add the result to the local coefficients
       **/
      template< class QuadratureType, class RangeVectorType, class JacobianRangeVectorType >
      void axpyQuadrature ( const QuadratureType &quad, const RangeVectorType& rangeVector, const JacobianRangeVectorType& jacobianVector )
      {
        basisFunctionSet().axpy( quad, rangeVector, jacobianVector, localDofVector() );
      }

      void bind ( const EntityType &entity )
      {
        bound_ = true;
        basisFunctionSet_ = discreteFunction().space().basisFunctionSet( entity );
        localDofVector().resize( basisFunctionSet().size() );
        assemblyOperation_.begin( entity, discreteFunction(), localDofVector() );
      }

      void unbind ()
      {
        if (bound_)
        {
          assemblyOperation_.end( entity(), localDofVector(), discreteFunction() );
          basisFunctionSet_ = BasisFunctionSetType();
        }
      }

      /** \brief return const reference to local DoF vector **/
      const LocalDofVectorType &localDofVector () const { return localDofVector_; }
      /** \brief return mutable reference to local DoF vector **/
      LocalDofVectorType &localDofVector () { return localDofVector_; }

    private:
      DiscreteFunctionType &discreteFunction_;
      LocalDofVectorType localDofVector_;
      BasisFunctionSetType basisFunctionSet_;
      AssemblyOperationType assemblyOperation_;
      bool bound_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH
