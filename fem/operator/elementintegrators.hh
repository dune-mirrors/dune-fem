/**************************************************************************
**       Title: general element integrators to be used in elliptic problems
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: several implementations of element-matrix and 
**              element-rhs-integrators, 
**              which can be accumulated for solving a general 
**              elliptic problems (or simplifications of it) with the FEOp 
**              class:
**
**   - div(a(x)*grad(u)) + div (b(x)*u) + c(x)*u = f(x)       in Omega
**                                             u = g_D(x)   in \Gamma_D
**                      (a(x)*grad(u) -b(x)*u) n = g_N(x)   in \Gamma_N
**            (a(x)*grad(u) -b(x)*u) n + alpha*u = g_R(x)   in \Gamma_R
**
**************************************************************************/
#ifndef DUNE_ELEMENTINTEGRATOR_HH
#define DUNE_ELEMENTINTEGRATOR_HH

#include <config.h>
#include <dune/grid/io/file/dgfparser/gridtype.hh>
// class required for lagrange-node-access with point() method:
// #include "lagrangedofhandler.hh"
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/operator/model/linearellipticmodel.hh>

namespace Dune 
{
/*! @defgroup CGFiniteElement Continuous Finite Elements
 *  @ingroup OperatorCommon
 * Description: several implementations of element-matrix and 
**              element-rhs-integrators, 
**              which can be accumulated for solving a general 
**              elliptic problems (or simplifications of it) with the FEOp 
**              class:
** \f{eqnarray*}
**   - div(a(x)*grad(u)) + div (b(x)*u) + c(x)*u &=& f(x)     \quad\mbox{in}\quad \Omega    \\
**                                             u &=& g_D(x)   \quad\mbox{in}\quad \Gamma_D \\
**                      (a(x)*grad(u) -b(x)*u) n &=& g_N(x)   \quad\mbox{in}\quad \Gamma_N \\
**            (a(x)*grad(u) -b(x)*u) n + alpha*u &=& g_R(x)   \quad\mbox{in}\quad \Gamma_R 
** \f}
** @{
**************************************************************************/



  /*! \class ElementMatrixIntegratorInterface
   *  \brief ElementMatrixIntegratorInterface is the interface for the 
   *         definition of a class providing local element matrices in a FEM 
   *         problem.
   *
   *  This is the interface for a general ElementMatrixIntegrator, i.e. a class
   *  providing a local element matrix. The interface does not contain 
   *  implementations, but indicates the minimum functionality to be defined by
   *  a ElementMatrixIntegrator class. Default implementations are provided in
   *  the DefaultElementMatrixIntegrator class, so the latter should be used 
   *  as base class for deriving own classes, which then can be used in the 
   *  FEOp finite element operator.
   *
   *  The class requires two template arguments: a traits and the derived type 
   *  for barton-nackman.
   *  The Traits class is not instantiated at any time
   *  only provides types. 
   *
   *  An FEOp using this local element matrix integrator therefore must iterate 
   *  over the grid, determine element matrices and accumulate them in a large 
   *  matrix.
   *
   *  An example and minimum requirements on the TraitsImp class can be found
   *  in elementintegratortraits.hh
   *
   *  An example of the use is given in dune/fem/examples/elliptic
   */
  template< class TraitsImp, class ModelImp, class ElementMatrixIntegratorImp >
  class ElementMatrixIntegratorInterface;
  
  template< class TraitsImp, class ModelImp, class ElementMatrixIntegratorImp >
  class ElementMatrixIntegratorInterface
    < TraitsImp, 
      LinearEllipticModelInterface< typename ModelImp :: FunctionSpaceType,
                                    ModelImp,
                                    typename ModelImp :: Properties >,
      ElementMatrixIntegratorImp
    >
  {
  public:
    typedef TraitsImp                                TraitsType;
    typedef LinearEllipticModelInterface< typename ModelImp :: FunctionSpaceType,
                                          ModelImp,
                                          typename ModelImp :: Properties >
      ModelType;
    //typedef ModelImp                                 ModelType;
    typedef typename TraitsType :: EntityType        EntityType;
    typedef typename TraitsType :: ElementMatrixType ElementMatrixType;

    /*!
     *  addElementMatrix: Interface method that adds a multiple of a local
     *                    element matrix on an entity.
     *
     *   Instead of pure computation and storing in the matrix, a scaled addition
     *   is performed on the current matrix, such that different
     *   ElementMatrices can be combined for more complex problem
     *
     *   This method has to be provided by derived classes.
     *
     *   \param[in] entity reference to the entity
     *   \param     matrix instance of the local matrix implementation 
     *   \param[in] coef   scaling coefficient
     */
    inline void addElementMatrix ( EntityType &entity,
                                   ElementMatrixType& matrix,
                                   double coefficient = 1 ) // const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().addElementMatrix( entity, matrix, coefficient ) );
    }

    inline const ModelType &model () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().model() );
      return asImp().model();
    }

  protected:
    inline const ElementMatrixIntegratorImp &asImp () const 
    { 
      return static_cast< const ElementMatrixIntegratorImp& >( *this );
    }

    inline ElementMatrixIntegratorImp &asImp ()
    { 
      return static_cast< ElementMatrixIntegratorImp& >( *this ); 
    }
  }; // end of ElementMatrixIntegratorInterface class



  /** \class DefaultElementMatrixIntegrator
   *
   *  \brief The DefaultElementMatrixIntegrator class implements some 
   *         default functionality.
   *
   *  The class can be used for deriving own ElementMatrixIntegrator classes. 
   *  Mainly it provides model access and storage functionality. It does not 
   *  implement the addElementMatrix method, so this still has to be implemented 
   *  in derived classes. But the class provides building blocks for use in 
   *  future addElementMatrix methods. In particular methods 
   *  addMassElementMatrix, addDiffusiveFluxElementMatrix,
   *  addConvectiveFluxElementMatrix.
   *
   *  So a derivation of an own class is very simple, e.g. done in
   *  dune/fem/examples/elliptic/elliptic.cc 
   */
  template <class TraitsImp, class ModelImp, class ElementMatrixIntegratorImp>
  class DefaultElementMatrixIntegrator
  : public ElementMatrixIntegratorInterface< TraitsImp, 
                                             typename ModelImp :: LinearEllipticModelInterfaceType,
                                             ElementMatrixIntegratorImp >
  {
  public:
    typedef DefaultElementMatrixIntegrator
      < TraitsImp, ModelImp, ElementMatrixIntegratorImp >
      ThisType;
    typedef ElementMatrixIntegratorInterface
      < TraitsImp,
        typename ModelImp :: LinearEllipticModelInterfaceType,
        ElementMatrixIntegratorImp >
      BaseType;
  
  public:
    typedef typename BaseType :: TraitsType TraitsType;
    typedef typename BaseType :: ModelType ModelType;
    typedef typename TraitsType::ElementMatrixType    ElementMatrixType;
    typedef typename TraitsType::EntityType           EntityType;
    typedef typename TraitsType::EntityPointerType    EntityPointerType;

    //! derived typedefs for ommitting long type-dereferencing in the methods
    typedef typename TraitsType :: ElementQuadratureType ElementQuadratureType;
    typedef typename TraitsType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename TraitsType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename TraitsType::RangeType RangeType;
    typedef typename TraitsType::DomainType DomainType;
    typedef typename TraitsType::JacobianRangeType JacobianRangeType;
    typedef typename TraitsType::IntersectionIteratorType 
                     IntersectionIteratorType;
    typedef typename TraitsType::IntersectionQuadratureType 
                     IntersectionQuadratureType;
    typedef typename TraitsType::GridPartType GridPartType;
    //! used grid type 
    typedef typename GridPartType :: GridType GridType; 

    //! dimension of world 
    enum { dimworld = GridType :: dimensionworld };

    /**  constructor: initialize ElementMatrixIntegrator with given model
     *
     *   Default implementation storing a reference to the 
     *   model with the data functions. The constructor addditionally 
     *   allocates storage for temporary basis-function gradients, which are
     *   required in case of diffusive contributions.
     *
     *   perhaps extension by a further parameter (global matrix) might be useful,
     *   if a local-matrix interface is used, which directly accesses the 
     *   global matrix?
     *
     *   \param[in] model model providing the (continuous) data
     *
     *   \param{in] dfSpace discrete function space
     *
     *   \param[in] verbose optional verbosity flag
     */
    DefaultElementMatrixIntegrator ( ModelType& model,
                                     const DiscreteFunctionSpaceType &dfSpace,
                                     int verbose = 0 )
    : model_( model ),
      discreteFunctionSpace_( dfSpace ),
      verbose_( verbose )
    {
      // determine number of basis functions on first entity 
      if (verbose_)
          std::cout << "entered constructor of " 
                    << " DefaultElementMatrixIntegrator\n";
      
      if (verbose_)
          std::cout << "got discrete functionspace\n";

      EntityPointerType ep = dfSpace.begin();
      EntityType& entity = *ep;     

      if (verbose_)
          std::cout << "got first entity\n";
      
//            const typename EntityType::Geometry& geo = entity.geometry();

      if (verbose_)
          std::cout << "successful got geometry\n";
      
      const BaseFunctionSetType& baseSet = 
        dfSpace.baseFunctionSet( entity );

      if (verbose_)
          std::cout << "got base function set\n";
      
      numBaseFunctions_ = baseSet.numBaseFunctions();
      gradPhiPtr_ = new JacobianRangeType[ numBaseFunctions_ ];

      if (verbose_)
          std::cout << "allocated temporary storage for gradients\n";
    }
    
    /*!  access function for model
     *
     *   Default implementation is return of the stored reference.
     *
     *   \return reference to the locally stored model reference
     */
    inline const ModelType& model () const
    {
      return model_;
    }

    /*!  access function for model
     *
     *   Default implementation is return of the stored reference.
     *
     *   \return reference to the locally stored model reference
     */
    inline ModelType& model () 
    {
      return model_;
    }

    inline const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    /*! 
     *   destructor: free of temporary memory for basis-fct-gradients 
     *               in mygrad
     */
    ~DefaultElementMatrixIntegrator ()
    {
      if( verbose_ )
        std::cout << "entered destructor of DefaultElementMatrixIntegrator"
                  << std :: endl;
      delete[] gradPhiPtr_;
    }
    
    inline void addElementMatrix ( EntityType &entity,
                                   ElementMatrixType& matrix,
                                   double coefficient = 1 ) // const
    {
      addDiffusiveFluxElementMatrix( entity, matrix, coefficient );
      if( ModelType :: Properties :: hasConvectiveFlux )
        addConvectiveFluxElementMatrix( entity, matrix, coefficient );
      if( ModelType :: Properties :: hasMass )
        addMassElementMatrix( entity, matrix, coefficient );
      
      if( entity.hasBoundaryIntersections() )
      {
        if( ModelType :: Properties :: hasGeneralizedNeumannValues )
          addGeneralizedNeumannElementMatrix( entity, matrix, coefficient );
        if( ModelType :: Properties :: hasRobinValues )
          addRobinElementMatrix( entity, matrix, coefficient );
      }
    }

    /*!
     *  addDiffusiveFluxElementMatrix: accumulate diffusive contributions
     *
     *  The method is used for the diffusive flux of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[
     *     L_ij :=  \int_\entity   [a     grad(phi_j) ]^T  grad(phi_i) 
     *  \f]
     *  The model class is assumed to have a diffusiveFlux() method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param[in]  entity the entity over which the intrgration is performed
     *
     *  \param[out] matrix reference ot the local element matrix storage to be
     *              increased
     *
     *  \param[in]  coefficient an optional weighting coefficient, which is
     *              multiplied to the increase before matrix addition
     */
    template< class EntityImp, class ElementMatrixImp >
    void addDiffusiveFluxElementMatrix( const EntityImp& entity,
                                        ElementMatrixImp& matrix, 
                                        double coefficient = 1 ) const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;

      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &discreteFunctionSpace = this->discreteFunctionSpace();

      const GeometryType &geometry = entity.geometry();
      
      // get local basis
      const BaseFunctionSetType& baseSet
        =  discreteFunctionSpace.baseFunctionSet( entity );
      int numBaseFunctions = baseSet.numBaseFunctions();

      // assert that allocated space for gradPhiPtr is sufficient!!
      assert( numBaseFunctions <= numBaseFunctions_ );
          
      // assert that matrix allocation is sufficient
      assert( matrix.rows() >= numBaseFunctions );
      assert( matrix.cols() >= numBaseFunctions );
          
      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree );
      const int numQuadraturePoints = quadrature.nop();
      for ( int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const CoordinateType &x = quadrature.point( pt );

        const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
        const double volume = geometry.integrationElement( x );
            
        for( int i = 0; i < numBaseFunctions; ++i )
        {
          JacobianRangeType &gradPhi = gradPhiPtr_[ i ];
          baseSet.jacobian( i, quadrature, pt, gradPhi );
          // multiply with transposed of the jacobian inverse 
          gradPhi[ 0 ] = FMatrixHelp :: mult( inv, gradPhi[ 0 ] );
        }
            
        // evaluate diffusiveFlux for all gradients of basis functions
        const double factor = coefficient * quadrature.weight( pt ) * volume;
        for( int j = 0; j < numBaseFunctions; ++j )
        {
          JacobianRangeType psi;
          model.diffusiveFlux( entity, quadrature, pt, gradPhiPtr_[ j ], psi );
          for( int i = 0; i < numBaseFunctions; ++i )
            matrix.add( i, j, factor * (psi[ 0 ] * gradPhiPtr_[ i ][ 0 ]) );
        }            
      }
    } // end addDiffusiveFluxElementMatrix

    /** \brief accumulate convective contributions
     *
     *  The method is used for the convective flux of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[                
     *     L_ij :=  \int_\entity   [-  b   phi_j]^T         grad(phi_i) 
     *  \f]
     *  The model class is assumed to have a convectiveFlux() method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param[in]  entity      entity over which the intrgration is performed
     *
     *  \param[out] matrix      matrix to be updated
     *
     *  \param[in]  coefficient optional weighting coefficient by which the
     *                          update is multiplied
     */
    template< class EntityType, class ElementMatrixType >
    void addConvectiveFluxElementMatrix ( EntityType &entity,
                                          ElementMatrixType &matrix,
                                          double coefficient = 1 ) // const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;

      assert( ModelType :: Properties :: hasConvectiveFlux );
      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace(); 

      const GeometryType &geometry = entity.geometry();
      
      const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
      const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
            
      assert( matrix.rows() >= numBaseFunctions );
      assert( matrix.cols() >= numBaseFunctions );
            
      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree );
      const int numQuadraturePoints = quadrature.nop();
      for( int pt=0; pt < numQuadraturePoints; ++pt )
      {
        const CoordinateType &x = quadrature.point( pt );
        // calc Jacobian inverse before volume is evaluated 
        const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
        const double volume = geometry.integrationElement( x );
              
        const double factor = coefficient * quadrature.weight( pt ) * volume;
        for( int i = 0; i < numBaseFunctions; ++i ) 
	      { 
          // evaluate gradient of base function
          JacobianRangeType gradPhi;
          baseFunctionSet.jacobian( i, quadrature, pt, gradPhi );  
		      gradPhi[ 0 ] = FMatrixHelp :: mult ( inv, gradPhi[ 0 ] );
		
          for( int j = 0; j < numBaseFunctions; ++j )
          {
            // evaluate base function
            RangeType phi;
            baseFunctionSet.evaluate( j, quadrature, pt, phi );
            phi *= factor;
            
            // evaluate convectiveFlux
            DomainType flux;
            model.convectiveFlux( entity, quadrature, pt, phi, flux );
          }
        }
      } // end loop over quadrature points            
    } // end method addConvectiveFluxElementMatrix
    
    /** \brief accumulate mass contributions
     *
     *  The method is used for the mass term of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[
     *     L_ij :=  \int_\entity   c          phi_i        phi_j
     *  \f]
     *  The model class is assumed to have a mass() method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param mat reference ot the local element matrix storage to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before matrix addition
     */
    template< class EntityType, class ElementMatrixType >
    void addMassElementMatrix ( EntityType &entity,
                                ElementMatrixType &matrix, 
                                double coefficient = 1 ) // const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;

      assert( ModelType :: Properties :: hasMass );
      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();

      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
      const int numBaseFunctions = baseFunctionSet.numBaseFunctions();

      //assert( numBaseFunctions <= numBaseFunctions_);
            
      assert( matrix.rows() >= numBaseFunctions );
      assert( matrix.cols() >= numBaseFunctions );
           
      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree );
      const int numQuadraturePoints = quadrature.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
	    {
        const CoordinateType &x = quadrature.point( pt );
	      const double volume = geometry.integrationElement( x );
	      
	      const double factor = coefficient * quadrature.weight( pt ) * volume;
	      for( int i = 0; i < numBaseFunctions; ++i )
        {
	        RangeType phi_i;
          baseFunctionSet.evaluate( i, quadrature, pt, phi_i );

          for( int j = 0; j < numBaseFunctions; ++j )
          {
	          RangeType phi_j;
            baseFunctionSet.evaluate( j, quadrature, pt, phi_j );
            phi_j *= factor;
            
    	      RangeType mass;
            model.mass( entity, quadrature, pt, mass );
            matrix.add( i, j, mass[ 0 ] * (phi_i[ 0 ] * phi_j[ 0 ]) );
          }
        }
      } // end loop over quadrature points
    } // end method addMassElementMatrix
    
    /*! \brief accumulate generalized Neumann boundary contributions
     *
     *  The method is used for the Neumann boundary of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[
     *     L_ij :=  +  \int_\Gamma_Neu alphaFunction (entity, iquad, pt,it,ret )  phi_i        phi_j     
     *  \f]
     *  
     *  The model class is assumed to implement the following methods
     *  \code
     *  generalizedNeumannAlpha( intersection, quadrature, point )
     *  boundaryType()
     *  \endcode
     *
     *  \param[in] entity entity over which the intrgration is performed
     *
     *  \param     matrix local element matrix storage to be updated
     *
     *  \param[in] coef   optional weighting coefficient by which the update
     *                    is multiplied
     */
    template< class EntityType, class ElementMatrixType >
    void addGeneralizedNeumannElementMatrix ( EntityType& entity,
                                              ElementMatrixType &matrix,
                                              double coef = 1 ) //  const
    {
      assert( ModelType :: Properties :: hasNeumannValues );
      const ModelType &model = this->model(); 

	    // for all intersections check whether boundary
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace(); 
      
      GridPartType &gridPart = dfSpace.gridPart();
           
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        // check whether this is a boundary
        if( !it.boundary() )
          continue;
        // check for generalized neumann boundary type
        if( model.boundaryType( it ) != ModelType :: GeneralizedNeumann )
          continue;
                   
        const BaseFunctionSetType &baseFunctionSet
          = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();     

        // integrate over intersection
        IntersectionQuadratureType 
          quadrature( it, TraitsType :: quadDegree, IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop();
        for ( int pt = 0; pt < numQuadraturePoints; ++pt ) 
        {  
		      // the following seems wrong...
          const double volume = it.intersectionGlobal().integrationElement(quadrature.localPoint(pt));
          const double alpha = model.generalizedNeumannAlpha( it, quadrature, pt );
          const double factor = coef * alpha * quadrature.weight( pt ) * volume;
                                            
          for(int i = 0; i < numBaseFunctions; ++i )
          {
            RangeType phi_i;
            baseFunctionSet.evaluate( i, quadrature, pt, phi_i );
            for( int j = 0; j < numBaseFunctions; ++j ) 
            {
			        RangeType phi_j;
      			  baseFunctionSet.evaluate( j, quadrature, pt, phi_j );
              matrix.add( i, j, factor * (phi_i[ 0 ] * phi_j[ 0 ]) );
			      }
          } 
        } // end loop over quadraturepoints
      } // end loop over intersections
    } // end method addGeneralizedNeumannElementMatrix

    /*! addRobinElementMatrix: accumulate Robin boundary contributions
     *
     *  The method is used for the robin boundary of general elliptic problems.
     *  The following matrix is computed, where i,j run over the local dofs
     *  of base functions, which have support on an entity.
     *  \f[
     *     L_ij :=  +  \int_\Gamma_R alpha      phi_i        phi_j     
     *  \f]
     *  The model class is assumed to have an alpha() and a
     *  boundaryType() member method.
     *
     *  the method must be a template method, such that the model requirements
     *  are only mandatory, if the method is instantiated.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param mat reference ot the local element matrix storage to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before matrix addition
     */
    template< class EntityType, class ElementMatrixType >
    void addRobinElementMatrix( EntityType &entity, 
                                ElementMatrixType &matrix,
                                double coefficient = 1 ) //  const
    {
      assert( ModelType :: Properties :: hasRobinValues );
      const ModelType &model = this->model();
      
      // for all intersections check whether boundary
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();
      
      const GridPartType &gridPart = dfSpace.gridPart();
            
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        // check for boundary
        if( !it.boundary() )
          continue;

        // check for robin boundary values
        if( model.boundaryType( it ) != ModelType :: Robin )
          continue;

        // integrate over intersection
        IntersectionQuadratureType
          quadrature( it, TraitsType :: quadDegree,
                      IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop();
                    
        const BaseFunctionSetType &baseFunctionSet
          = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
                    
        for ( int pt = 0; pt < numQuadraturePoints; ++pt ) 
        {
          // the following seems wrong...
          const double volume = it.intersectionGlobal().integrationElement(quadrature.localPoint(pt));
          const double alpha = model.robinAlpha( it, quadrature, pt );
          const double factor = coefficient * alpha * quadrature.weight( pt ) * volume;
                                            
          for( int i = 0; i < numBaseFunctions; ++i ) 
          {
            RangeType phi_i;
            baseFunctionSet.evaluate( i, quadrature, pt, phi_i );
            for( int j = 0; j < numBaseFunctions; ++j ) 
			      {
              RangeType phi_j;
              baseFunctionSet.evaluate( j, quadrature, pt, phi_j );
              matrix.add( i, j, factor * (phi_i[ 0 ] * phi_j[ 0 ]) );
            }
          }
        } // end loop over quadraturepoints
      } // end loop over intersections
    } // end method addRobinElementMatrix
    
  protected:
    //! Reference to the data model
    ModelType& model_;
    //! the discrete function space
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    //! verbosity flag
    int verbose_;
    //! number of basis functions
    int numBaseFunctions_;  
    //! storage for basis-function gradients
    JacobianRangeType *gradPhiPtr_;
  }; // end of DefaultElementMatrixIntegrator class



  /*! \class ElementRhsIntegratorInterface 
   *  \brief The ElementRhsIntegratorInterface specifies the methods, which
   *         must be implemented by an ElementRhsIntegrator
   *
   *  This is the Interface class for a local RhsIntegrator. An instance of
   *  a derived class can be used as template parameter in the RhsAssembler
   *  class for assembling the Rhs of a general problem.
   *
   *  In derived classes, the method addElementRhs must be implemented.
   *
   *  No default implementations are done here, an example of a derived class is
   *  the DefaultElementRhsIntegrator
   *
   *  The class uses two template parameters, which is a Traitsclass and 
   *  the final derived class by Barton-Nackman
   */
  template< class TraitsImp, class ModelImp, class ElementRhsIntegratorImp >
  class ElementRhsIntegratorInterface;

  template< class TraitsImp, class ModelImp, class ElementRhsIntegratorImp >
  class ElementRhsIntegratorInterface
    < TraitsImp,
      LinearEllipticModelInterface< typename ModelImp :: FunctionSpaceType,
                                    ModelImp,
                                    typename ModelImp :: Properties >,
      ElementRhsIntegratorImp
    >
  {
  public:
    typedef TraitsImp TraitsType;
    typedef LinearEllipticModelInterface< typename ModelImp :: FunctionSpaceType,
                                          ModelImp,
                                          typename ModelImp :: Properties >
      ModelType;
    
    /*! 
     *   addElementRhs: Interface method that adds a multiple of a local rhs 
     *                on an entity. 
     *
     *   This method must be implemented in derived classes. 
     *
     *   The RhsAssembler class, in which the ElementRhsIntegrator can be used,
     *   uses this class, by a grid-walkthrough and collecting all local 
     *   contributions 
     * 
     *   \param entity the entity over which integration is performed
     *
     *   \param elRhs storage, which is to be increased, writable access by 
     *          operator[] should be possible, i.e. a "LocalFunction" access
     *          is assumed. 
     *
     *   \param coef an optonal coefficient, which is multiplied to the 
     *               increment value before addition
     */
    template< class EntityType, class ElementRhsType >
    inline void addElementRhs ( EntityType &entity, 
                                ElementRhsType &elRhs, 
                                double coefficient = 1 ) // const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().addElementRhs( entity, elRhs, coefficient ) );
    }

    inline const ModelType& model () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().model() );
      return asImp().model();
    }
    
  protected:
    inline const ElementRhsIntegratorImp &asImp () const 
          { 
            return static_cast< const ElementRhsIntegratorImp& >( *this ); 
          }

    inline ElementRhsIntegratorImp &asImp ()
    { 
       return static_cast< ElementRhsIntegratorImp& >( *this ); 
    }
  };



  /** \class DefaultElementRhsIntegrator 
   *  \brief The DefaultElementRhsIntegrator provides implementation of some 
   *         default methods 
   *
   *  For using a derived class, a method addElementRhs must be 
   *  implemented there.
   *  This can be done explicitly, or simply by combining existing methods,
   *  an example use of the class in application as given in 
   *  dune/fem/examples/elliptic/elliptic.cc. 
   *  Then it can be used in the RhsAssembler class.
   *
   *  The class uses three template parameters, a Traitsclass, a 
   *  ModelImp and a ElementRhsIntegratorImp 
   *
   *  The TraitsType typedef collects various typdefs and 
   *
   *  The ModelImp class is a model class providing various data functions, 
   *  which define the model for the elliptic problem. An instance of the 
   *  Model is required in the constructor for the DefaultRhsIntegrator.
   *
   *  The ElementRhsIntegratorImp template parameter is the final derived 
   *  class for Barton-Nackman.
   */
  template< class TraitsImp, class ModelImp, class ElementRhsIntegratorImp >
  class DefaultElementRhsIntegrator
  : public ElementRhsIntegratorInterface
    < TraitsImp, 
      typename ModelImp :: LinearEllipticModelInterfaceType,
      ElementRhsIntegratorImp
    >
  {
  private:
    typedef DefaultElementRhsIntegrator< TraitsImp, ModelImp, ElementRhsIntegratorImp >
      ThisType;
    typedef ElementRhsIntegratorInterface
      < TraitsImp,
        typename ModelImp :: LinearEllipticModelInterfaceType,
        ElementRhsIntegratorImp >
      BaseType;
    
  public:
    typedef typename BaseType :: TraitsType TraitsType;
    typedef typename BaseType :: ModelType ModelType;
   
    typedef typename TraitsType::ElementQuadratureType ElementQuadratureType;
    typedef typename TraitsType::IntersectionQuadratureType 
                     IntersectionQuadratureType;
    typedef typename TraitsType::DiscreteFunctionSpaceType 
                     DiscreteFunctionSpaceType;
    typedef typename TraitsType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename TraitsType::DiscreteFunctionType::LocalFunctionType 
                     LocalFunctionType;
    typedef typename TraitsType::RangeType RangeType;
    typedef typename TraitsType::GridPartType GridPartType;
    typedef typename TraitsType::IntersectionIteratorType 
                     IntersectionIteratorType;

  // the following does not seem to work properly...
  //  enum {quadDegree = TraitsType::quadDegree};
    
  /*======================================================================*/
  /*! 
   *   constructor: store local reference to the underlying model
   *
   *   The local reference model_ is initialized
   *
   *   \param model the model which provides the data functions
   *
   *   \param verbose an optional verbosity flag
   */
  /*======================================================================*/  
    DefaultElementRhsIntegrator( ModelType &model, const DiscreteFunctionSpaceType &dfSpace, int verbose = 0):
            model_(model), discreteFunctionSpace_( dfSpace ), verbose_(verbose)
    {
      if (verbose_)
          std::cout << "entered constructor of " 
                    << "DefaultElementRhsIntegrator\n";
    };

    /** \breif give access to the model
     *
     *   the reference to the locally stored model is returned
     */
    inline ModelType &model() 
    {
      return model_;
    }

    inline const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    template< class EntityType, class ElementRhsType >
    inline void addElementRhs ( EntityType &entity, 
                                ElementRhsType &elRhs, 
                                double coefficient = 1 ) // const
    {
      if( ModelType :: Properties :: hasSource )
        addSourceElementRhs( entity, elRhs, coefficient );

      if( entity.hasBoundaryIntersection() )
      {
        if( ModelType :: Properties :: hasNeumannValues )
          addNeumannElementRhs( entity, elRhs, coefficient );
        if( ModelType :: Properties :: hasRobinValues )
          addRobinElementRhs( entity, elRhs, coefficient );
        if( ModelType :: Properties :: hasGeneralizedNeumannValues )
          addGeneralizedNeumannElementRhs( entity, elRhs, coefficient );
      }
    }
 
    /** \brief adds a multiple of a source term contribution to the rhs 
     *
     *   This method can be used as a building block for problem specific 
     *   ElementRhsIntegrators.
     * 
     *   The elRhs vector is increased by the following values:
     *   \f[
     *      b_i + = coef * \int_{entity}  source * phi_i 
     *   \f]
     *   By accessing only the local basis functions, only few 
     *   values b_i are updated.
     *
     *   If this function is used, the model must provide a 
     *   source() method
     * 
     *   by keeping the method templatized, the model requirements really only 
     *   are demanding, if these building blocks are used. So also simple 
     *   models can be realized.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param elRhs reference ot the local function to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before addition
     */
    template< class EntityType, class ElementRhsType >
    void addSourceElementRhs( EntityType &entity, 
                              ElementRhsType &elRhs,
                              double coefficient = 1 ) // const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;

      assert( ModelType :: Properties :: hasSource );
      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();

      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
      const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
      
      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree);
      const int numQuadraturePoints = quadrature.nop();
            
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
      {
        const CoordinateType &x = quadrature.point( pt );

        const double volume = geometry.integrationElement( x );
        const double factor = coefficient * quadrature.weight( pt ) * volume;
              
        for( int i = 0; i < numBaseFunctions; ++i ) 
        {
          RangeType phi;
          baseFunctionSet.evaluate( i, quadrature, pt, phi );
          
          RangeType source;
          model.source( entity, quadrature, pt, source );

          elRhs[ i ] += factor * (source[ 0 ] * phi[ 0 ]);
        }
      }
    } // end method addSourceElementRhs


    /**  \brief adds a multiple of a Neuman term contribution to the rhs 
     *
     *   This method can be used as a building block for problem specific 
     *   ElementRhsIntegrators.
     * 
     *   The elRhs vector is increased by the following values:
     *   \f[
     *      b_i + = coef * \int_{neumann_boundary of entity} neumannValues * phi_i 
     *   \f]
     *   By accessing only the local basis functions, only few 
     *   values b_i are updated.
     *
     *   If this function is used, the model must provide a 
     *   neumannValues and a boundaryType method
     *
     *   by keeping the method templatized, the model requirements really only 
     *   are demanding, if these building blocks are used. So also simple 
     *   models can be realized.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param elRhs reference ot the local function to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before addition
     */
    template< class EntityType, class ElementRhsType >
    void addNeumannElementRhs ( EntityType &entity, 
                                ElementRhsType &elRhs, 
                                double coefficient = 1 ) // const
    {
      assert( ModelType :: Properties :: hasNeumannValues );
      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();
      
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        // check for boundary
        if( !it.boundary() )
          continue;
        // check for Neumann boundary
        if( model.boundaryType( it ) != ModelType :: Neumann )
          continue;

        const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
 
        // integrate over intersection
        IntersectionQuadratureType
          quadrature( it, TraitsType :: quadDegree, IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop();
        for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
        {
          // the following seems wrong...
          const double volume = it.intersectionGlobal().integrationElement(quadrature.localPoint(pt));
          const double factor = coefficient * quadrature.weight( pt ) * volume;
        
          for( int i = 0; i < numBaseFunctions; ++i )
          {
            RangeType phi;
            baseFunctionSet.evaluate( i, quadrature, pt, phi );  
            
            RangeType value;
            model.neumannValues( it, quadrature, pt, value );
            
            elRhs[ i ] += factor * (value[ 0 ] * phi[ 0 ]);
          }
        } // end loop over quadrature points
      } // end loop over intersections
    }

    /**  adds a multiple of a Robin term contribution to the rhs 
     *
     *   This method can be used as a building block for problem specific 
     *   ElementRhsIntegrators.
     * 
     *   The elRhs vector is increased by the following values:
     *   \f[
     *      b_i + = coef * \int_{robin_boundary of entity} robinValues * phi_i 
     *   \f]
     *   By accessing only the local basis functions, only few 
     *   values b_i are updated.
     *
     *   If this function is used, the model must provide a 
     *   robinValues and a boundaryType method
     *
     *   by keeping the method templatized, the model requirements really only 
     *   are demanding, if these building blocks are used. So also simple 
     *   models can be realized.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param elRhs reference ot the local function to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before addition
     */
    template< class EntityType, class ElementRhsType >
    void addRobinElementRhs ( EntityType &entity, 
                              ElementRhsType &elRhs, 
                              double coefficient = 1 ) // const
    {
      assert( ModelType :: Properties :: hasRobinValues );
      const ModelType &model = this->model();
      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();
     
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        if( !it.boundary() )
          continue;
        if( model.boundaryType( it ) != ModelType :: Robin )
          continue;

        const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions = baseFunctionSet.numBaseFunctions();
 
        // integrate over intersection
        IntersectionQuadratureType
          quadrature( it, TraitsType :: quadDegree, IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop();
        for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
        {  
          // the following seems wrong...
          const double volume = it.intersectionGlobal().integrationElement(quadrature.localPoint(pt));                 
          const double factor = coefficient * quadrature.weight( pt ) * volume;
        
          for( int i = 0; i < numBaseFunctions; ++i ) 
          {
            RangeType phi;
            baseFunctionSet.evaluate( i, quadrature, pt, phi );  
            
            RangeType value;
            model.robinValues( entity, quadrature, pt, value );
            
            elRhs[ i ] += factor * (value[ i ] * phi[ i ]);
          }
        } // end loop over quadrature points
      } // end loop over intersections
    }

    /**  adds a multiple of a Robin term contribution to the rhs 
     *
     *   This method can be used as a building block for problem specific 
     *   ElementRhsIntegrators.
     * 
     *   The elRhs vector is increased by the following values:
     *   \f[
     *      b_i + = coef * \int_{neumann_boundary of entity} generalizedNeumannValues * phi_i 
     *   \f]
     *   By accessing only the local basis functions, only few 
     *   values b_i are updated.
     *
     *   If this function is used, the model must provide a 
     *   robinValues and a boundaryType method
     *
     *   by keeping the method templatized, the model requirements really only 
     *   are demanding, if these building blocks are used. So also simple 
     *   models can be realized.
     *
     *  \param entity the entity over which the intrgration is performed
     *
     *  \param elRhs reference ot the local function to be increased
     *
     *  \param coef an optional weighting coefficient, which is multiplied to the
     *         increase before addition
     */
    template< class EntityType, class ElementRhsType >
    void addGeneralizedNeumannElementRhs( EntityType &entity,
                                          ElementRhsType &elRhs,
                                          double coefficient = 1 ) // const
    {
      assert( ModelType :: Properties :: hasGeneralizedNeumannValues );
      const ModelType &model = this->model();

      const DiscreteFunctionSpaceType &dfSpace = this->discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();
      
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        if( !it.boundary() )
          continue;
        if( model.boundaryValues( it ) != ModelType :: GeneralizedNeumann )
          continue;
        
        const BaseFunctionSetType &baseFunctionSet = dfSpace.baseFunctionSet( entity );
        const int numBaseFunctions =  baseFunctionSet.numBaseFunctions();
 
        // integrate over intersection
        IntersectionQuadratureType
          quadrature( it, TraitsType :: quadDegree, IntersectionQuadratureType :: INSIDE );
        const int numQuadraturePoints = quadrature.nop(); 
        for( int pt = 0; pt < numQuadraturePoints; ++pt ) 
        {  
          // the following seems wrong...
          const double volume = it.intersectionGlobal().integrationElement(quadrature.localPoint(pt));
          const double factor = coefficient * quadrature.weight( pt ) * volume;
        
          for( int i=0; i < numBaseFunctions; ++i ) 
          {
            RangeType phi;
            baseFunctionSet.evaluate( i, quadrature, pt, phi );
            
            RangeType value;
            model.generalizedNeumannValues( it, quadrature, pt, value );
            elRhs[ i ] += factor * (value[ 0 ] * phi[ 0 ]);
          }
        } // end loop over quadrature points
      } // end loop over intersections
    } //end of addGeneralizedNeumannElementRhs
  
  private:
    //! reference to the underlying model specified during construction
    ModelType& model_;
    //! the discrete function space
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    //! verbosity flag;
    int verbose_;
  }; // end class DefaultElementRhsIntegrator



  /** \class RhsAssembler
   *  \brief The RhsAssembler class assembles the right hand side of a 
   *         general finite element problem including boundary correction
   *
   *  The right hand side for a Lagrange-basis of order 1 is assembled, i.e.
   *  \f{eqnarray*}
   *    b_i := \left\{\begin{array}{ll}
   *           g_D(x_i)   &  \mbox{if}\; x_i \;\mbox{is Dirichlet-Bnd.Point}
   *           \\
   *           XX_i       &  \mbox{otherwise}
   *           \end{array}\right.
   *  \f}
   *
   *  where \f$ XX_i \f$ is determined by grid-walkthorugh and calling of the 
   *  ElementRhsIntegrator on each Element.
   *  For instance for a full general elliptic problem
   *  \f{eqnarray*}
   *                 - div(a*grad(u) - b*u) + c*u &=& f   \quad\mbox{in}\quad \Omega  \\
   *                                            u &=& g_D \quad\mbox{in}\quad \Gamma_D\\
   *                           (a*grad(u) -b*u) n &=& g_N \quad\mbox{in}\quad \Gamma_N\\
   *                 (a*grad(u) -b*u) n + alpha*u &=& g_R \quad\mbox{in}\quad \Gamma_R
   *            (a*grad(u) -b*u) n + alphaFunction*u &=& g_GN \quad\mbox{in}\quad \Gamma_neu
   *  \f}
   *
   *  where \f$ a,b,c,g_D,g_N,g_R \f$ are space dependent, alpha a constant and the 
   *  quantities denote
   *              "stiffness"        a(x) 
   *              "velocity"         b(x) 
   *              "mass"             c(x) 
   *              "source"           f(x) 
   *              "dirichletValues"  g_D
   *              "neumannValues"    g_N
   *              "robinValues"      g_R
   *              "alpha"            alpha
   *              "generalizedNeumannValues"   g_GN
   *              "alphaFunction"    alphaFunction                   
                     
   *
   *  the EllipticElementRhsIntegrator in fem/examples/elliptic/elliptic.cc 
   *  results in
   *  \f[
   *    XX_i =    \int_\Omega   f   phi_i
   *            + \int_\Gamma_N g_N phi_i
   *            + \int_\Gamma_R g_R phi_i                   
   *  \f]
   *  The corresponding Matrix is generated by FEOp. Additionally that class 
   *  allows an additional dirichlet-boundary treatment of both the matrix and 
   *  the Rhs
   *
   *  The class requires a Template Argument ElementRhsIntegrator, which must 
   *  provide a TraitsType typedef collecting various typenames and a ModelType
   *
   *  An Instance of the ElementRhsIntegrator is required during 
   *  initialization, and must provide a method addElementRhs(entity, rhs, coef) 
   *  The RhsAssembler class assembles the right hand side by performing a 
   *  grid-walkthrough and repeatedly calling this allElementRhs method. 
   *  finally, a dirichlet treatment is performed by setting dirichlet-DOFs
   *  to the Dirichlet-values as specified by the model (known to the 
   *  ElementRhsIntegrator) 
   */
template <class ElementRhsIntegratorImp>
class RhsAssembler
{
public:  
  typedef ElementRhsIntegratorImp ElementRhsIntegratorType;
  typedef typename ElementRhsIntegratorType::TraitsType TraitsType;
  typedef typename ElementRhsIntegratorType::ModelType ModelType;
  typedef typename TraitsType::DiscreteFunctionType 
          DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridType GridType; 
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType 
          IntersectionIteratorType;
  typedef typename TraitsType::IntersectionQuadratureType 
          IntersectionQuadratureType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename TraitsType::RangeType RangeType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename EntityType :: ctype coordType;
    
/*======================================================================*/
/*! 
 *   constructor: store a reference to the elementRhsIntegrator
 *
 *   \param elRhsInt the ElementRhsIntegrator to be used
 *
 *   \param verbose an optional verbosity flag
 */
/*======================================================================*/
  
  RhsAssembler(ElementRhsIntegratorType& elRhsInt, int verbose = 0)
          : elRhsInt_(elRhsInt), verbose_(verbose)
        {
          if (verbose_)
              std::cout << "entered model() of " 
                        << "RhsAssembler\n";
        };

/*======================================================================*/
/*! 
 *   assemble: assemble the rhs for a general elliptic problem
 *
 *   \param rhs a reference to the rhs, which is to be filled
 */
/*======================================================================*/
  
  void assemble(DiscreteFunctionType& rhs)
        {
          if (verbose_)
              std::cout << "entered assemble() of " 
                        << "RhsAssembler\n";
          rhs.clear();
          
          // grid walkthrough for accumulating rhs
          const DiscreteFunctionSpaceType &dfSpace 
            =  elRhsInt_.discreteFunctionSpace();
          
          IteratorType it    = dfSpace.begin(); 
          IteratorType endit = dfSpace.end(); 
          
          for (;it!=endit;++it)
          {
            // initialize ElementRhsStorage : get LocalFunction
            LocalFunctionType elRhs = rhs.localFunction(*it);            
            elRhsInt_.addElementRhs(*it,elRhs);
          }
          
          // final dirichlet treatment
          this->bndCorrect(rhs);
          if (verbose_)
              std::cout << "finished bndCorrect(rhs) \n";
        };
  
private:

/*======================================================================*/
/*! 
 *   bndCOrrect: boundary correction for the rhs of a general elliptic problem
 *
 *   method is called from assemble() after the grid-walkthrough-accumulation
 *
 *   \param rhs a reference to the rhs, which is to be corrected
 */
/*======================================================================*/

    void bndCorrect( DiscreteFunctionType &rhs )
    {
      if( verbose_ >= 1 )
        std::cout << "entered bndCorrect() of RhsAssembler" << std :: endl;

      // codimension of faces
      enum { faceCodim = 1 };
      
      // type of Lagrange point set
      // note: DiscreteFunctionSpaceType must be a Lagrange discrete function
      //       space
      typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
        LagrangePointSetType;
      // type of iterator over all DoFs on a face
      typedef typename LagrangePointSetType :: template Codim< 1 >
                                            :: SubEntityIteratorType
        FaceDofIteratorType;

      // grid walkthrough for setting Dirichlet-DOFs
      ModelType &model = elRhsInt_.model();
      const DiscreteFunctionSpaceType &fspace = elRhsInt_.discreteFunctionSpace();
      
      const GridPartType & gridPart = fspace.gridPart();
                    
      IteratorType it = fspace.begin();
      const IteratorType endit = fspace.end();
      for( ; it != endit; ++it ) {
        const EntityType &entity = *it;

        // const GeometryType t = entity.geometry().type();
        IntersectionIteratorType nit = gridPart.ibegin( entity );
        const IntersectionIteratorType endnit = gridPart.iend( entity );
        for( ; nit != endnit; ++nit ) {
          if( !nit.boundary() )
            continue;
          if( model.boundaryType( nit ) != ModelType :: Dirichlet )
            continue;
          
          const int faceNumber = nit.numberInSelf();
          LocalFunctionType lf = rhs.localFunction( entity );
          
          const LagrangePointSetType &lagrangePointSet
            = fspace.lagrangePointSet( entity );
          FaceDofIteratorType faceIt
            = lagrangePointSet.template beginSubEntity< faceCodim >( faceNumber );
          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.template endSubEntity< faceCodim >( faceNumber );
          for( ; faceIt != faceEndIt; ++faceIt ) {
            const unsigned int entityDofNumber = *faceIt;
            RangeType phi;
            model.dirichletValues( nit, lagrangePointSet, entityDofNumber, phi );
            lf[ entityDofNumber ] = phi[ 0 ];
          }
         
          #if 0
          // determine all local DOF numbers of intersection vertices
          // enum { dim = EntityType :: dimension };
                  
          // const ReferenceElement<coordType,dim>& refElem; 
          // ReferenceElement<coordType,dim> refElem; 
          // ReferenceElementContainer<coordType,dim> refElemCont;
          // const ReferenceElement<coordType,dim>& 
          // refElem = refElemCont(t);
                  
          // t is geometrytype

          // class, which gives access to local vertex coordinates by a
          // point() routine:
          // VertexPointProvider<EntityType> fakequad(en);
          LagrangeDofHandler< DiscreteFunctionSpaceType >
            dofHandler( fspace, entity );

          // int novx = refElem.size( faceNumber, faceCodim, dim );
          int novx = dofHandler.numDofsOnFace( faceNumber, faceCodim );
          //  assert( novx == dim );
                  
          // we only can treat lagrange-functions with correspondence 
          // between local vertices and basis-functions: 
          // assert(lf.numDofs() == refElem.size(dim));

          for( int j = 0; j < novx; ++j ) {
            // get all local DOF number located on the face 
                    
            // int vx  = refElem.subEntity(face, faceCodim , j , dim );
            int vx = dofHandler.entityDofNum( faceNumber, faceCodim, j );
                    
            // determine Dirichlet-value in this vertex
            RangeType dirval;
            // std::cout << "j =" << j << ", vx = " << vx << " \n";   
            model.dirichletValues( entity, dofHandler, vx, dirval ); 
            // std::cout << "dirval[0] =" << dirval[0] << "\n";
            // set local DOF to determined value
            // here the assumption of identical local numbering of 
            // DOFS and vertex numbering enters!!!!
            lf[ vx ] = dirval;                    
            // std::cout << "lf[vx] =" << lf[vx] << "\n";
          }
          // std::cout << "finished loop of lf-update \n";
          #endif
        }
        // std::cout << "finished loop over intersections \n";
      }
      // std::cout << "finished loop over entities \n";
    }
  
// member variables:
private:
  //! a reference to the elementRhsIntegrator, which is to be used
  ElementRhsIntegratorType& elRhsInt_;
  //! verbosity flag
  int verbose_;
}; // end class RhsAssembler
//! @}
} // end namespace dune

#endif

