#ifndef SPACE_P1BUBBLE_HH
#define SPACE_P1BUBBLE_HH

#include <vector>
#include <dune/common/exceptions.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/localkey.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/space/basisfunctionset/default.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>

namespace Dune
{

  namespace Fem
  {
    constexpr int factorial(unsigned int n) { return n <= 1 ? 1 : n * factorial(n-1); }

    // The following would normalize the max of the bubble to 1
    // std::pow( dimDomain + 1.0, dimDomain + 1.0 );

    // This normalizes the bubble to have mean-value 1 on the refernce element.
    constexpr int bubbleNormalization(int dimDomain) { return factorial(2*dimDomain + 1); }

    template< int dim >
    struct BubbleElementLocalKeyMap
    {
      //! [Constructor of LocalKey tripple]
      BubbleElementLocalKeyMap ( int vertices )
      {
        for( int i = 0; i < vertices; ++i )
          map_.emplace_back( i, dim, 0 );
        map_.emplace_back( 0, 0, 0 );
      }
      //! [Constructor of LocalKey tripple]

      std::size_t size() const { return map_.size(); }

      LocalKey& localKey ( std::size_t i ) { return map_[ i ]; }
      const LocalKey& localKey ( std::size_t i ) const { return map_[ i ]; }

    private:
      std::vector< LocalKey > map_;
    };

    struct BubbleElementDofMapperCodeFactory
    {
      // return the shape functions for a given reference element. If this
      // is not possible an empty DofMapperCode is returned.
      template< class RefElement,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().size( 0 ) ) >, int >::value, int > = 0,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().type( 0, 0 ) ) >, GeometryType >::value, int > = 0 >
      DofMapperCode operator() ( const RefElement &refElement ) const
      {
        static const int dim = RefElement::dimension;
        if( refElement.type().isSimplex() )
          return compile( refElement, BubbleElementLocalKeyMap< dim >(dim+1) );
        if( refElement.type().isCube() && refElement.type().dim() == 2)
          return compile( refElement, BubbleElementLocalKeyMap< dim >(pow(2,dim)) );
        else
          return DofMapperCode();
      }
    };

    // SimplexBubbleElementShapeFunctionSet
    // ----------------------------------

    template< class FunctionSpace >
    class SimplexBubbleElementShapeFunctionSet
    {
      typedef SimplexBubbleElementShapeFunctionSet< FunctionSpace > ThisType;
    public:
      static const int dimDomain =  FunctionSpace :: dimDomain;
      static const int polynomialOrder = dimDomain + 1;
      static const int numShapeFunctions = dimDomain + 2;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpace :: DomainType DomainType;
      typedef typename FunctionSpace :: RangeType RangeType;
      static_assert( RangeType::dimension == 1, "This is a scalar shapefunction set" );
      typedef typename FunctionSpace :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpace :: HessianRangeType HessianRangeType;

      SimplexBubbleElementShapeFunctionSet () {}

      template<class GeometryType >
      SimplexBubbleElementShapeFunctionSet ( const GeometryType& gt )
      {
        if( !gt.isSimplex() )
          DUNE_THROW( NotImplemented, "Wrong geometry type for Cube2DBubblelementShapeFunctionSet." );
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );
        RangeType phi(1), phi0(1);
        for( int i=0; i< dimDomain; ++i )
        {
          functor( i+1, RangeType( xRef[ i ] ) );
          phi0[ 0 ] -= xRef[ i ];
          phi[ 0 ] *= xRef[ i ] ;
        }

        phi[ 0 ] *= phi0[ 0 ] / bubbleNormalization(dimDomain); //  std::pow( ( dimDomain + 1.0 ), dimDomain + 1.0 );
        functor( 0, phi0 );
        functor( dimDomain +1, phi );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );

        JacobianRangeType jac(0), jac0( -1 );
        RangeType phi0( 1 );

        functor( 0, jac0 );

        for( int i=0; i< dimDomain; ++i )
        {
          phi0[ 0 ] -= xRef[ i ];

          for( int j=1; j < dimDomain; ++j )
            jac0[ 0 ][ (i+j)%dimDomain ] *= xRef[ i ];

          jac[ 0 ][ i ] = 1;
          functor( i+1, jac );
          jac[ 0 ][ i ] = 0;
        }

        for( int i=0; i< dimDomain; ++i )
          jac0[ 0 ][ i ] *= -(phi0[ 0 ] - xRef[ i ]);
        jac0[ 0 ] *= 1.0 / bubbleNormalization(dimDomain); // ::pow( dimDomain + 1.0, dimDomain + 1.0 );
        functor( dimDomain +1, jac0 );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DUNE_THROW( NotImplemented, "NotImplemented" );
        DomainType xRef = coordinate( x );
        HessianRangeType hes;
        functor( 0, hes );
        functor( 1, hes );
        functor( 2, hes );
        functor( 3, hes );
      }

      int order () const { return dimDomain + 1; }

      std::size_t size () const { return dimDomain +2; }
    };

    // 2d cube element taken from
    // http://arxiv.org/pdf/1507.04417v1.pdf
    template< class FunctionSpace >
    class Cube2DBubbleElementShapeFunctionSet
    {
      typedef Cube2DBubbleElementShapeFunctionSet< FunctionSpace > ThisType;
    public:
      static const int dimDomain =  FunctionSpace :: dimDomain;
      static const int polynomialOrder = dimDomain + 1;
      static const int numShapeFunctions = dimDomain + 2;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpace :: DomainType DomainType;
      typedef typename FunctionSpace :: RangeType RangeType;
      static_assert( RangeType::dimension == 1, "This is a scalar shapefunction set" );
      typedef typename FunctionSpace :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpace :: HessianRangeType HessianRangeType;

      //! [Main methods for shape functions]
      Cube2DBubbleElementShapeFunctionSet () {}

      template<class GeometryType >
      Cube2DBubbleElementShapeFunctionSet ( const GeometryType& gt )
      {
        if( !gt.isCube() )
          DUNE_THROW( NotImplemented, "Wrong geometry type for Cube2DBubblelementShapeFunctionSet." );
        if( gt.dim() != 2 )
          DUNE_THROW( NotImplemented, "2DCubeBubbleElementShapeFunctionSet only implemented for dimension = 2." );
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      //! [Main methods for shape functions]
      {
        DomainType xRef = coordinate( x );
        RangeType phi(1), phi0(1);
        // (1-x)(1-y)  -> grad = ( (y-1),(x-1) )
        phi0[0] = (1.-xRef[0])*(1.-xRef[1]);
        functor( 0, phi0 );
        // x(1-y)      -> grad = ( (1-y),-x )
        phi[0] = xRef[0]*(1.-xRef[1]);
        functor( 1, phi );
        // (1-x)y      -> grad = ( (-y,(1-x) )
        phi[0] = (1.-xRef[0])*xRef[1];
        functor( 2, phi );
        // xy          -> grad = ( (y,x) )
        phi[0] = xRef[0]*xRef[1];
        functor( 3, phi );
        // 64 xy phi_0(x,y)^2       -> grad = 128 phi_0 xy grad_0 +
        //                                     64 phi_0^2 (y,x)
        phi[0] = 64.*phi0[0]*phi0[0]*xRef[0]*xRef[1];
        functor( 4, phi );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );

        JacobianRangeType jac, jac0;
        RangeType phi0;
        phi0[0] = (1.-xRef[0])*(1.-xRef[1]);
        jac0[0] = {xRef[1]-1,xRef[0]-1};

        functor( 0, jac0 );
        jac[0] = {1-xRef[1],xRef[0]};
        functor( 1, jac );
        jac[0] = {xRef[1],1-xRef[0]};
        functor( 2, jac );
        jac[0] = {xRef[1],xRef[0]};
        functor( 3, jac );
        // 64 xy phi_0(x,y)^2       -> grad = 128 phi_0 xy grad_0 +
        //                                     64 phi_0^2 (y,x)
        jac0[0] *= 128.*phi0*xRef[0]*xRef[1];
        jac[0] = {xRef[1],xRef[0]};
        jac[0] *= 64.*phi0*phi0;
        jac[0] += jac0[0];
        functor( 4, jac );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DUNE_THROW( NotImplemented, "NotImplemented" );
        DomainType xRef = coordinate( x );
        HessianRangeType hes;
        functor( 0, hes );
        functor( 1, hes );
        functor( 2, hes );
        functor( 3, hes );
      }

      int order () const { return 4; }

      std::size_t size () const { return 5; }
    };

    template< class FunctionSpace, class GridPart, class Storage >
    class BubbleElementSpace;

    // BubbleElementSpaceTraits
    // ----------------------

    template< class FunctionSpace, class GridPart, class Storage >
    struct BubbleElementSpaceTraits
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
    public:
      typedef typename FunctionSpaceType :: ScalarFunctionSpaceType ScalarFunctionSpaceType;

      // defined in shapefunctionset.hh
      // First check if we have a grid with only one element type
      // (hybrid // grids are not supported with this space yet)
      // If it only has one get the topologyId of that element type
      // (note simplex=0 and cube=2^dim-1)
      static_assert( Dune::Capabilities::hasSingleGeometryType<typename GridPartType::GridType>::v,
                     "bubble elements only implemented for grids with a single geometry type" );
      static const unsigned int topologyId =
        Dune::Capabilities::hasSingleGeometryType<typename GridPartType::GridType>::topologyId;
      // now choose either the simplex or the cube implementation of the
      // bubble space
      typedef SimplexBubbleElementShapeFunctionSet< ScalarFunctionSpaceType > ScalarSimplexShapeFunctionSetType;
      typedef Cube2DBubbleElementShapeFunctionSet< ScalarFunctionSpaceType > ScalarCubeShapeFunctionSetType;
      typedef typename std::conditional< topologyId == 0,
              ScalarSimplexShapeFunctionSetType, ScalarCubeShapeFunctionSetType >::type
                ScalarShapeFunctionSetType;
      // now extend it to a vector valued funcion space
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;
      // and add some default functionality
      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
      // finished with the shape function set

      static const int localBlockSize = FunctionSpaceType::dimRange;
      static const int polynomialOrder = ScalarShapeFunctionSetType::polynomialOrder;

      typedef IndexSetDofMapper< GridPartType > BlockMapperType;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };

    template< class FunctionSpace, class GridPart, class Storage >
    struct LocalBubbleElementInterpolation
    {
      typedef BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > TraitsType;

      typedef typename FunctionSpace::DomainType DomainType;
      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeType RangeType;
      static const int dimDomain = FunctionSpace::dimDomain;
      static const int dimRange = FunctionSpace::dimRange;

      LocalBubbleElementInterpolation ()
        : points_( dimDomain + 2, DomainType( 0.0 ) )
      {
        for ( int i = 0; i < dimDomain; ++i )
          points_[ i + 1 ][ i ] = 1.0;

        points_[ dimDomain +1 ] = DomainType( 1.0 / ( dimDomain + 1.0 ) );
      }

      LocalBubbleElementInterpolation ( const LocalBubbleElementInterpolation & ) = default;
      LocalBubbleElementInterpolation ( LocalBubbleElementInterpolation && ) = default;

      //! [Evaluation of local interpolations]
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &lf, LocalDofVector &ldv ) const
      //! [Evaluation of local interpolations]
      {
        static_assert( TraitsType::topologyId == 0, "p1Bubble interpolation is only implemented for simplicial grids.");

        RangeType mean(0);

        // standard Lagrange P1 interpolation on the vertices
        for (int dof = 0; dof < dimDomain + 1; ++dof) {
          const auto& x = points_[dof];
          RangeType phi;
          lf.evaluate( x, phi );
          for( int i = 0; i < dimRange; ++i ) {
            ldv[ dimRange * dof + i ] = phi[ i ];
            mean[i] += phi[i];
          }
        }

        // Now for the bubble. If we assume that the bubble is used
        // for the sole purpose to implement the so called
        // Mini-Element, then we should have a look at the proof for
        // its stability. There it show that the bubble DOF is used in
        // order to ensured that the interpolant of a function f has
        // the same mean-value on each element as f.
        RangeType centerValue;
        lf.evaluate(points_[dimDomain+1], centerValue);
        int dof = dimDomain + 1;
        for( int i = 0; i < dimRange; ++i ) {
          // explicit mid-point rule
          ldv[ dimRange * dof + i ] = centerValue[i]  - mean[i] / (dimDomain + 1.0);
        }
      }

      template <class Entity>
      void bind( const Entity& entity ) {}
      void unbind() {}

    private:
      std::vector< DomainType > points_;
      mutable DomainFieldType entityVolume_;
    };

    // BubbleElementSpace
    // ----------------

    //! [Class definition for new space]
    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    class BubbleElementSpace
    : public DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > >
    //! [Class definition for new space]
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

    public:
      typedef BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > Traits;
      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      static const int polynomialOrder = Traits::polynomialOrder;

      typedef typename BaseType::GridPartType GridPartType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      // type of local interpolation
      typedef LocalBubbleElementInterpolation< FunctionSpace, GridPart, Storage > InterpolationType;

      // static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      static const InterfaceType defaultInterface = GridPartType::indexSetInterfaceType;
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      BubbleElementSpace ( GridPartType &gridPart,
                         const InterfaceType commInterface = defaultInterface,
                         const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_(  gridPart, BubbleElementDofMapperCodeFactory() )
      {
      }

      ~BubbleElementSpace ()
      {
      }

      bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there
        return blockMapper().contains( codim );
      }

      bool continuous () const { return true; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType & intersection ) const
      {
        // forward to the subsapces
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template<class Entity>
      int order ( const Entity &entity ) const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, ShapeFunctionSetType( entity.geometry().type() ) );
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const { return blockMapper_; }

      // Non-interface methods

      /** \brief return local interpolation object for LocalInterpolation
       */
      InterpolationType interpolation () const
      {
        return InterpolationType();
      }

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       *  \note this method is needed to call the global function inteprolate( ... )
       */
      [[deprecated]]
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return interpolation();
      }

    protected:
      mutable BlockMapperType blockMapper_;
    };

    template< class FunctionSpace,
              class GridPart,
              class Storage,
              class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< BubbleElementSpace < FunctionSpace, GridPart, Storage >, NewFunctionSpace >
    {
      typedef BubbleElementSpace< NewFunctionSpace, GridPart, Storage > Type;
    };

    template< class FunctionSpace, class GridPart, class Storage >
    class DefaultLocalRestrictProlong < BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    : public EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    {
      typedef typename FunctionSpace::DomainType DomainType;
      typedef typename FunctionSpace::RangeType RangeType;
      static const int dimDomain = FunctionSpace::dimDomain;
      static const int dimRange = FunctionSpace::dimRange;

      typedef EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > > BaseType;

      typedef BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > SpaceTraitsType;

      public:
      DefaultLocalRestrictProlong( const BubbleElementSpace< FunctionSpace, GridPart, Storage > &space )
      : BaseType(), points_( dimDomain + 2, DomainType( 0.0 ) ), childCounter_(-1)
      {
        for ( int i = 0; i < dimDomain; ++i )
          points_[ i + 1 ][ i ] = 1.0;

        points_[ dimDomain +1 ] = DomainType( 1.0 / ( dimDomain + 1.0 ) );
      }

      template< class LFParent, class LFChild, class LocalGeometry >
      void restrictLocal ( LFParent &lfParent,
                           const LFChild &lfChild,
                           const LocalGeometry &geometryInParent,
                           const bool initialize ) const
      {
        static_assert( SpaceTraitsType::topologyId == 0, "p1Bubble restriction and prolongation is only implemented for simplicial grids.");

        // Lagrange P1 on a simplex: restriction just means: forget

        // now for the bubble, if used as Stokes discretization the
        // bubble should adjust the local mean-value. Mmmh. Ok. Let's
        // just do something which is not too far off the mark but
        // much simpler ... the value of the center DOF of the parent
        // is just the mean value of two center DOFs of the
        // children. This way refinement and coarsening which directly
        // follow each other result in just the same function
        // (although the fine-grid function is not identical to the
        // course grid function)
        int dof = dimDomain + 1;
        if (initialize) {
          childCounter_ = 1;
          for( int coordinate = 0; coordinate < dimRange; ++coordinate ) {
            lfParent[ dimRange * dof + coordinate ] = lfChild[ dimRange * dof + coordinate ];
          }
        } else {
          ++ childCounter_;
          for( int coordinate = 0; coordinate < dimRange; ++coordinate ) {
            lfParent[ dimRange * dof + coordinate ] += lfChild[ dimRange * dof + coordinate ];
          }
        }
      }

      template< class LFParent >
      void restrictFinalize ( LFParent &lfParent ) const
      {
        int dof = dimDomain + 1;
        for( int coordinate = 0; coordinate < dimRange; ++coordinate ) {
          lfParent[ dimRange * dof + coordinate ] /= (double)childCounter_;
        }
#ifndef NDEBUG
        childCounter_ = -1;
#endif
      }

      template< class LFParent, class LFChild, class LocalGeometry >
      void prolongLocal ( const LFParent &lfParent, LFChild &lfChild,
                          const LocalGeometry &geometryInParent,
                          bool initialize ) const
      {
        static_assert( SpaceTraitsType::topologyId == 0, "p1Bubble restriction and prolongation is only implemented for simplicial grids.");

        // P1 Lagrange prolongation: need to install coordinates into the child LF
        for ( int dof = 0; dof < dimDomain + 1; ++dof) {
          const DomainType &pointInChild = points_[dof];
          const DomainType pointInParent = geometryInParent.global( pointInChild );
          RangeType phi;
          lfParent.evaluate( pointInParent, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate ) {
            lfChild[ dimRange * dof + coordinate ] = phi[ coordinate ];
          }
        }

        // but now for the bubble: just copy the DOF value to the children
        int dof = dimDomain + 1;
        for( int coordinate = 0; coordinate < dimRange; ++coordinate ) {
          lfChild[ dimRange * dof + coordinate ] = lfParent[ dimRange * dof + coordinate ];
        }
      }

      bool needCommunication () const { return true; }
    private:
      std::vector< DomainType > points_;
      mutable int childCounter_;
    };

  } // namespace Fem

} // namespace Dune

#endif // SPACE_P1BUBBLE_HH
