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
  template< class TraitsImp, class ElementMatrixIntegratorImp >
  class ElementMatrixIntegratorInterface 
  {
  public:
    typedef TraitsImp                                 TraitsType;
    typedef typename TraitsType::EntityType           EntityType;
    typedef typename TraitsType::ElementMatrixType    ElementMatrixType;


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
     *   \param[in] entity a reference to an 
     *   \param mat instance of the local matrix implementation 
     *   \param[in] coef the scaling coefficient
     */
    void addElementMatrix ( EntityType &entity,
                            ElementMatrixType& mat, 
                            double coef = 1 ) // const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().addElementMatrix( entity, mat, coef ) );
      // don't call the implementation twice!
      //asImp().addElementMatrix( entity, mat, coef);
    }

  protected:
    //! Barton-Nackman method forwarding
    ElementMatrixIntegratorImp& asImp ()
    { 
      return static_cast< ElementMatrixIntegratorImp& >( *this ); 
    }

    //! Barton-Nackman method forwarding
    const ElementMatrixIntegratorImp &asImp( ) const 
    { 
      return static_cast< const ElementMatrixIntegratorImp& >( *this );
    }
  }; // end of ElementMatrixIntegratorInterface class



/*======================================================================*/
/*!
 *  \class DefaultElementMatrixIntegrator
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
/*======================================================================*/
  
  template <class TraitsImp, class ModelImp, class ElementMatrixIntegratorImp>
  class DefaultElementMatrixIntegrator:
        public ElementMatrixIntegratorInterface<TraitsImp, 
                                                ElementMatrixIntegratorImp>
  {
  public:
//! basic typedefs 
    typedef TraitsImp                                 TraitsType;
    typedef ModelImp                                  ModelType;
    typedef typename TraitsType::ElementMatrixType    ElementMatrixType;
    typedef typename TraitsType::EntityType           EntityType;
    typedef typename TraitsType::EntityPointerType    EntityPointerType;

//! derived typedefs for ommitting long type-dereferencing in the methods
    typedef typename TraitsType::ElementQuadratureType ElementQuadratureType;
    typedef typename TraitsType::DiscreteFunctionSpaceType 
            DiscreteFunctionSpaceType;
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

/*======================================================================*/
/*! 
 *   constructor: initialize ElementMatrixIntegrator with given model
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
 *   \param model an instance of the model providing the data functions
 *
 *   \param verbose an optional verbosity flag
 */
/*======================================================================*/

    DefaultElementMatrixIntegrator(ModelType& model, int verbose = 0): 
            model_(model), verbose_(verbose) 
          {
            // determine number of basis functions on first entity 
            if (verbose_)
                std::cout << "entered constructor of " 
                          << " DefaultElementMatrixIntegrator\n";
            
            DiscreteFunctionSpaceType& fspace = 
                this->model().discreteFunctionSpace();    

            if (verbose_)
                std::cout << "got discrete functionspace\n";

            EntityPointerType ep = fspace.begin();
            EntityType& entity = *ep;     

            if (verbose_)
                std::cout << "got first entity\n";
            
//            const typename EntityType::Geometry& geo = entity.geometry();

            if (verbose_)
                std::cout << "successful got geometry\n";
            
            const BaseFunctionSetType& baseSet = 
              fspace.baseFunctionSet( entity );

            if (verbose_)
                std::cout << "got base function set\n";
            
            numBaseFunctions_ = baseSet.numBaseFunctions();
            gradPhiPtr_ = new JacobianRangeType[ numBaseFunctions_ ];

            if (verbose_)
                std::cout << "allocated temporary storage for gradients\n";
          };
    
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
                                        double coefficient = 1.0 ) const
    {
      typedef typename ElementQuadratureType :: CoordinateType CoordinateType;
        
      // get quadrature and function space
      ElementQuadratureType quadrature( entity, TraitsType :: quadDegree );
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = model_.discreteFunctionSpace();

      // get local basis
      const BaseFunctionSetType& baseSet
        =  discreteFunctionSpace.baseFunctionSet( entity );
          
      // assert that allocated space for gradPhiPtr is sufficient!!
      int numBaseFunctions = baseSet.numBaseFunctions();
      assert( numBaseFunctions <= numBaseFunctions_ );
          
      // assert that matrix allocation is sufficient
      assert( matrix.rows() >= numBaseFunctions );
      assert( matrix.cols() >= numBaseFunctions );
          
      const int numQuadraturePoints = quadrature.nop();
      for ( int pt = 0; pt < numQuadraturePoints; ++pt ) {
        // const CoordinateType &x = quadrature.point( pt );

        // get the geometry of the entity
        // const typename EntityImp :: Geometry &geometry = entity.geometry();
 
        // compute gradients of all base functions in point pt
        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimworld, dimworld > &inv
          = entity.geometry().jacobianInverseTransposed( quadrature.point( pt ) );
        const double volume
          = entity.geometry().integrationElement( quadrature.point( pt ) );
            
        for( int i = 0; i < numBaseFunctions; ++i ) {
          JacobianRangeType &gradPhi = gradPhiPtr_[ i ];
          baseSet.jacobian( i, quadrature, pt, gradPhi );
          // multiply with transposed of the jacobian inverse 
          gradPhi[ 0 ] = FMatrixHelp :: mult( inv, gradPhi[ 0 ] );
        }
            
        // evaluate diffusiveFlux for all gradients of basis functions
        const double factor = quadrature.weight( pt ) * volume * coefficient;
        for( int j = 0; j < numBaseFunctions; ++j ) {
          JacobianRangeType psi;
          model_.diffusiveFlux( entity, quadrature, pt, gradPhiPtr_[ j ], psi );
          for( int i = 0; i < numBaseFunctions; ++i ) {
            const double incr = factor * (psi[ 0 ] * gradPhiPtr_[ i ][ 0 ]);
            matrix.add( i, j, incr );
          }
        }            
      }
    } // end addDiffusiveFluxElementMatrix

/*======================================================================*/
/*!
 *  addConvectiveFluxElementMatrix: accumulate convective contributions
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
 *  \param entity the entity over which the intrgration is performed
 *
 *  \param mat reference ot the local element matrix storage to be increased
 *
 *  \param coef an optional weighting coefficient, which is multiplied to the
 *         increase before matrix addition
 */
/*======================================================================*/
    
    template <class EntityImp, class ElementMatrixImp>
    void addConvectiveFluxElementMatrix(EntityImp& entity, 
                                        ElementMatrixImp& mat, 
                                        double coef = 1.0) // const
          {
            // get quadrature and function space
            ElementQuadratureType 
                quad(entity,TraitsType::quadDegree);
            DiscreteFunctionSpaceType& fspace = 
                this->model().discreteFunctionSpace();
            
            // get local basis
            const BaseFunctionSetType &baseSet = 
              fspace.baseFunctionSet( entity );
            
            int numBaseFunctions =  baseSet.numBaseFunctions();             
            // assert that allocated space for gradPhiPtr is sufficient!!
            //assert( numBaseFunctions <= numBaseFunctions_);
            
            // assert that matrix allocation is sufficient
            assert( mat.rows() >= numBaseFunctions );
            assert( mat.cols() >= numBaseFunctions );
            
            for ( int pt=0; pt < quad.nop(); pt++ ) 
            {  
              // calc Jacobian inverse before volume is evaluated 
              const FieldMatrix<double,dimworld,dimworld>& inv = 
                  entity.geometry().jacobianInverseTransposed(quad.point(pt));
              const double vol = 
                  entity.geometry().integrationElement(quad.point(pt));
              
              //for(int i=0; i<numBaseFunctions; i++) 
              //{
              //  baseSet.jacobian(i,quad,pt,gradPhiPtr_[i]);  
              //  // multiply with transpose of jacobian inverse 
              //  gradPhiPtr_[i][0] = FMatrixHelp :: mult ( inv,gradPhiPtr_[i][0] );
              //}
              
              JacobianRangeType gradPhi;
              DomainType ret;
              RangeType phi;
              
              double fact = quad.weight( pt ) * vol * coef;
              for(int i=0; i<numBaseFunctions; i++) 
	      {	    
		// evaluate gradient of base function
		baseSet.jacobian(i,quad,pt,gradPhi);  
		// multiply with transpose of jacobian inverse 
		gradPhi[0] = FMatrixHelp :: mult ( inv,gradPhi[0] );
		
		for (int j=0; j<numBaseFunctions; j++ )
                {
                  // evaluate base function
                  baseSet.evaluate(j,quad,pt,phi);  
                  // evaluate convectiveFlux
                  phi *= fact;
                  this->model().convectiveFlux(entity, quad, pt,  
                                               phi, 
                                               ret);
                  // does this vector-vector product work?
                  double incr =  ret * gradPhi[0];
                  mat.add(i,j, incr );
                } // end inner loop over basis functions            
              }// end outer loop over basis functions            
            } // end loop over quadrature points            
          }; // end method addConvectiveFluxElementMatrix


    
/*======================================================================*/
/*!
 *  addMassElementMatrix: accumulate mass contributions
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
/*======================================================================*/
    
    template <class EntityImp, class ElementMatrixImp>
    void addMassElementMatrix(EntityImp& entity, 
                              ElementMatrixImp& mat, 
                              double coef = 1.0) // const
          {
            // get quadrature and function space
            ElementQuadratureType 
                quad(entity,TraitsType::quadDegree);
            DiscreteFunctionSpaceType& fspace = 
                this->model().discreteFunctionSpace();
            
            // get local basis
            const BaseFunctionSetType &baseSet = 
              fspace.baseFunctionSet( entity );
            
            int numBaseFunctions =  baseSet.numBaseFunctions();             
            // assert that allocated space for gradPhiPtr is sufficient!!
            //assert( numBaseFunctions <= numBaseFunctions_);
            
            // assert that matrix allocation is sufficient
            assert( mat.rows() >= numBaseFunctions );
            assert( mat.cols() >= numBaseFunctions );
            
            for ( int pt=0; pt < quad.nop(); pt++ ) 
	    {  
	      // calc Jacobian inverse before volume is evaluated 
	      // const FieldMatrix<double,dimworld,dimworld>& inv = 
	      //   entity.geometry().jacobianInverseTransposed(quad.point(pt));
	      const double vol = 
                  entity.geometry().integrationElement(quad.point(pt));
	      
	      //for(int i=0; i<numBaseFunctions; i++) 
	      //{
	      //  baseSet.jacobian(i,quad,pt,gradPhiPtr_[i]);  
	      //  // multiply with transpose of jacobian inverse 
	      //  gradPhiPtr_[i][0] = FMatrixHelp :: mult ( inv,gradPhiPtr_[i][0] );
	      //}
	      
	      RangeType phi_i;
	      RangeType ret;
	      RangeType phi_j;
	      
	      double fact = quad.weight( pt ) * vol * coef;
	      for(int i=0; i<numBaseFunctions; i++) 
              {	    
                baseSet.evaluate(i,quad,pt,phi_i);  
                // evaluate gradient of base function
                //baseSet.jacobian(i,quad,pt,gradPhi);  
                // multiply with transpose of jacobian inverse 
                //gradPhi = FMatrixHelp :: mult ( inv,gradPhi );
                
                for (int j=0; j<numBaseFunctions; j++ )
                {
                  // evaluate basis function
                  baseSet.evaluate(j,quad,pt,phi_j); 
                  phi_j *= fact;
                  // evaluate mass function
                  this->model().mass(entity, quad, pt,  
                                     ret);
                  double incr =  ret[0] * phi_i[0] * phi_j[0];
                  mat.add(i,j, incr );
                }            
              }
	    }
          }; // end method addMassElementMatrix
/*======================================================================*/
/*!
 *  addGeneralizedNeumannElementMatrix: accumulate generalized Neumann boundary 
 *  contributions
 *
 *  The method is used for the Neumann boundary of general elliptic problems.
 *  The following matrix is computed, where i,j run over the local dofs
 *  of base functions, which have support on an entity.
 *  \f[
 *     L_ij :=  +  \int_\Gamma_Neu alphaFunction (entity, iquad, pt,it,ret )  phi_i        phi_j     
 *  \f]
 *  The model class is assumed to have an alphaFunction(entity, iquad, pt,it,ret) and a
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
/*======================================================================*/

    template <class EntityImp, class ElementMatrixImp>
    void addGeneralizedNeumannElementMatrix(EntityImp& entity, 
                               ElementMatrixImp& mat, 
                               double coef = 1.0) //  const
          {
	    // for all intersections check whether boundary
            DiscreteFunctionSpaceType& fspace = 
                this->model().discreteFunctionSpace();          
            GridPartType & gridPart = fspace.gridPart();
            
            IntersectionIteratorType 
                endit = gridPart.iend(entity);
            for(IntersectionIteratorType it = 
                    gridPart.ibegin(entity); 
                it != endit; ++it)
                
	
                if(it.boundary())
                {
                  // if boundary, get cog of intersection and check for Robin
                  IntersectionQuadratureType 
                      iquad(it,0,
                            IntersectionQuadratureType::INSIDE);

                  // check, whether real cog integration is performed
                  assert(iquad.nop()==1);
                  if (this->model().boundaryType(entity,iquad,0) == 
                      TraitsType::GeneralizedNeumann)
                  {
                    // if Robin-boundary, then integrate over intersect:
                    // get quadrature and function space
		 
                    IntersectionQuadratureType 
                        iquad(it,TraitsType::quadDegree, 
                              IntersectionQuadratureType::INSIDE);
                    
                    // get local basis
                    const BaseFunctionSetType & baseSet = 
                        fspace.baseFunctionSet(entity);
                    
                    int numBaseFunctions =  baseSet.numBaseFunctions();     
                    // LocalFunctionType lf = rhs.localFunction(entity);
                    
                    for ( int pt=0; pt < iquad.nop(); pt++ ) 
                    {  
		      // the following seems wrong...
                      const double vol = 
                          it.intersectionGlobal().
                          integrationElement(iquad.localPoint(pt));    
RangeType ret;
                      double fact = iquad.weight( pt ) * vol * coef ;
		      this->model().alphaFunction(entity, iquad, pt,  it,
						  ret );
			fact *= ret[0];	     
                                            
                      for(int i=0; i<numBaseFunctions; i++) 
                      {
		                           RangeType phi_i;
                        baseSet.evaluate(i,iquad,pt,phi_i);  
                        for(int j=0; j<numBaseFunctions; j++) 
			{
			
			  RangeType phi_j;
			  baseSet.evaluate(j,iquad,pt,phi_j);  	    
			  double incr =  fact * phi_i[0] * phi_j[0];
			  mat.add(i,j,incr);
			} // end inner loop over basefunctions           
                      } // end outer loop over basefunctions           
                    } // end loop over quadraturepoints
                  } // end of if isrobin-intersection
                } // end of isboundary 
          }; // end method addGeneralizedNeumannElementMatrix

/*======================================================================*/
/*!
 *  addRobinElementMatrix: accumulate Robin boundary contributions
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
/*======================================================================*/

    template <class EntityImp, class ElementMatrixImp>
    void addRobinElementMatrix(EntityImp& entity, 
                               ElementMatrixImp& mat, 
                               double coef = 1.0) //  const
          {
            // for all intersections check whether boundary
            DiscreteFunctionSpaceType& fspace = 
                this->model().discreteFunctionSpace();          
            GridPartType & gridPart = fspace.gridPart();
            
            IntersectionIteratorType 
                endit = gridPart.iend(entity);
            for(IntersectionIteratorType it = 
                    gridPart.ibegin(entity); 
                it != endit; ++it)
                
                if(it.boundary())
                {
                  // if boundary, get cog of intersection and check for Robin
                  IntersectionQuadratureType 
                      iquad(it,0,
                            IntersectionQuadratureType::INSIDE);

                  // check, whether real cog integration is performed
                  assert(iquad.nop()==1);
                  if (this->model().boundaryType(entity,iquad,0) == 
                      TraitsType::Robin)
                  {
                    // if Robin-boundary, then integrate over intersect:
                    // get quadrature and function space
                    IntersectionQuadratureType 
                        iquad(it,TraitsType::quadDegree, 
                              IntersectionQuadratureType::INSIDE);
                    
                    // get local basis
                    const BaseFunctionSetType &baseSet = 
                      fspace.baseFunctionSet( entity );
                    
                    int numBaseFunctions =  baseSet.numBaseFunctions();     
                    // LocalFunctionType lf = rhs.localFunction(entity);
                    
                    for ( int pt=0; pt < iquad.nop(); pt++ ) 
                    {  
                      // the following seems wrong...
                      const double vol = 
                          it.intersectionGlobal().
                          integrationElement(iquad.localPoint(pt));    

                      double fact = iquad.weight( pt ) * vol * coef * 
                          (this->model().alpha());
                                            
                      for(int i=0; i<numBaseFunctions; i++) 
                      {
                        RangeType phi_i;
                        baseSet.evaluate(i,iquad,pt,phi_i);  
                        for(int j=0; j<numBaseFunctions; j++) 
			{
			  RangeType phi_j;
			  baseSet.evaluate(j,iquad,pt,phi_j);  	    
			  double incr =  fact * phi_i[0] * phi_j[0];
			  mat.add(i,j,incr);
			} // end inner loop over basefunctions           
                      } // end outer loop over basefunctions           
                    } // end loop over quadraturepoints
                  } // end of if isrobin-intersection
                } // end of isboundary 
          }; // end method addRobinElementMatrix
    
  protected:
    //! Reference to the data model
    ModelType& model_;
    //! verbosity flag
    int verbose_;
    //! number of basis functions
    int numBaseFunctions_;  
    //! storage for basis-function gradients
    JacobianRangeType *gradPhiPtr_;
  }; // end of DefaultElementMatrixIntegrator class
  

/*======================================================================*/
/*!
 *  \class ElementRhsIntegratorInterface 
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
/*======================================================================*/

  template <class TraitsImp, class ElementRhsIntegratorImp>
  class ElementRhsIntegratorInterface
  {
  public:
    typedef TraitsImp                                 TraitsType;
    
/*======================================================================*/
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
/*======================================================================*/
    
    template <class EntityType, class ElementRhsType>
    void addElementRhs(EntityType &entity, 
                       ElementRhsType &elRhs, 
                       double coef = 1.0) // const
          {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(  (asImp().template 
                                              addElementRhs<EntityType, 
                                              ElementRhsType>
                                              ( entity, elRhs, coef) ) );
            asImp().template addElementRhs<EntityType,ElementRhsType>
                ( entity, elRhs, coef);
          }
    
  protected:
    // Barton-Nackman
    ElementRhsIntegratorImp &asImp() 
          { 
            return static_cast<ElementRhsIntegratorImp&>( *this ); 
          }
    
    const ElementRhsIntegratorImp &asImp( ) const 
          { 
            return static_cast<const ElementRhsIntegratorImp&>( *this ); 
          }
  };
  
/*======================================================================*/
/*!
 *  \class DefaultElementRhsIntegrator 
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
/*======================================================================*/

  template <class TraitsImp, class ModelImp, class ElementRhsIntegratorImp>
  class DefaultElementRhsIntegrator:
        public ElementRhsIntegratorInterface<TraitsImp,  
                                             ElementRhsIntegratorImp>
  {
public:
  typedef ModelImp ModelType;
  typedef typename ModelImp::TraitsType TraitsType;
 
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

  DefaultElementRhsIntegrator(ModelType& model, int verbose = 0):
          model_(model), verbose_(verbose)
  {
    if (verbose_)
        std::cout << "entered constructor of " 
                  << "DefaultElementRhsIntegrator\n";
  };

/*======================================================================*/
/*! 
 *   model: give access to the model for higher classes
 *
 *   the reference to the locally stored model is returned
 */
/*======================================================================*/  

  inline ModelType& model() 
  {
//  this method is simply called too often
//    if (verbose_)
//        std::cout << "entered model() of " 
//                  << "DefaultElementRhsIntegrator";
    return model_;
  };

/*======================================================================*/
/*! 
 *   addSourceElementRhs: adds a multiple of a source term contribution
 *                      to the rhs 
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
/*======================================================================*/  
    
    template <class EntityImp, class ElementRhsImp>
    void addSourceElementRhs(EntityImp &entity, 
                             ElementRhsImp &elRhs, 
                             double coef = 1.0) // const
        {
          // get quadrature and function space
          ElementQuadratureType 
              quad(entity,TraitsType::quadDegree);
          DiscreteFunctionSpaceType& fspace = 
              model_.discreteFunctionSpace();
          
          // get local basis
          const BaseFunctionSetType &baseSet = 
            fspace.baseFunctionSet( entity );
          
          // assert that allocated space for gradPhiPtr is sufficient!!
          int numBaseFunctions =  baseSet.numBaseFunctions();             
          // assert( numBaseFunctions <= numBaseFunctions_);
          
          // assert that matrix allocation is sufficient
          // assert( mat.rows() >= numBaseFunctions );
          // assert( mat.cols() >= numBaseFunctions );
          
          for ( int pt=0; pt < quad.nop(); pt++ ) 
          {  
            const double vol = 
                entity.geometry().integrationElement(quad.point(pt));
            double fact = quad.weight( pt ) * vol * coef;
            
            RangeType ret;
            RangeType phi;
            for(int i=0; i<numBaseFunctions; i++) 
            {
              baseSet.evaluate(i,quad,pt,phi);
              this->model().source(entity, quad, pt, ret);
              double incr =  ret[0] * phi[0] * fact;
              elRhs[i]+= incr;
            }  // end loop DOFs          
          } // end loop quadrature points
        }; // end method addSourceElementRhs


/*======================================================================*/
/*! 
 *   addNeumannElementRhs: adds a multiple of a Neuman term contribution
 *                       to the rhs 
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
/*======================================================================*/  

    template <class EntityImp, class ElementRhsImp>
    void addNeumannElementRhs(EntityImp &entity, 
                             ElementRhsImp &elRhs, 
			  double coef = 1.0) // const
  {
    // for all intersections check whether boundary
    DiscreteFunctionSpaceType& fspace = 
      this->model().discreteFunctionSpace();          
    GridPartType & gridPart = fspace.gridPart();
    
    IntersectionIteratorType 
      endit = gridPart.iend(entity);
    for(IntersectionIteratorType it = 
	  gridPart.ibegin(entity); 
	it != endit; ++it)
      
      if(it.boundary())
	{
	  // if boundary, get cog of intersection and check for Neumann
	  IntersectionQuadratureType 
	    iquad(it,0,IntersectionQuadratureType::INSIDE);
	  // check, whether real cog integration is performed
	  assert(iquad.nop()==1);
	  if (this->model().boundaryType(entity,iquad,0) 
              == TraitsType::Neumann)
	    {
	      // if Neumann-boundary, then integrate over intersect:
	      // get quadrature and function space
	      IntersectionQuadratureType 
		iquad(it,TraitsType::quadDegree, 
                      IntersectionQuadratureType::INSIDE);
	      
	      // get local basis
	      const BaseFunctionSetType &baseSet = 
		    fspace.baseFunctionSet( entity );
	      
	      int numBaseFunctions =  baseSet.numBaseFunctions();             
	      	      
	      for ( int pt=0; pt < iquad.nop(); pt++ ) 
		{  
		  // the following seems wrong...
		  const double vol = 
                      it.intersectionGlobal().
                      integrationElement(iquad.localPoint(pt)); 
//		    iquad.geometry().integrationElement(
//		      iquad.localPoint(pt));
		  // entity.geometry().integrationElement(iquad.point(pt));
		  
		  double fact = iquad.weight( pt ) * vol * coef;
		  
		  // get normal
		  // TraitsType::DomainType outerNormal = 
		  // it.unitOuterNormal(iquad.localPoint(pt));
		  
		  for(int i=0; i<numBaseFunctions; i++) 
                    {
		      RangeType phi;
		      baseSet.evaluate(i,iquad,pt,phi);  
		      
   		      RangeType ret;
		      this->model().neumannValues(entity, iquad, pt, ret);
                      ret[0] *= fact;
		      double incr =  ret[0] * phi[0];
		      elRhs[i] += incr;
                    } // end loop over basefunctions           
		} // end loop over quadraturepoints
	    } // end of if isdirichlet-intersection
	} // end of isboundary
  };

/*======================================================================*/
/*! 
 *   addRobinElementRhs: adds a multiple of a Robin term contribution
 *                     to the rhs 
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
/*======================================================================*/  
    
    template <class EntityImp, class ElementRhsImp>
    void addRobinElementRhs(EntityImp &entity, 
                            ElementRhsImp &elRhs, 
                            double coef = 1.0) // const
          {
    // for all intersections check whether boundary
    DiscreteFunctionSpaceType& fspace = 
      this->model().discreteFunctionSpace();          
    GridPartType & gridPart = fspace.gridPart();
    
    IntersectionIteratorType 
      endit = gridPart.iend(entity);
    for(IntersectionIteratorType it = 
	  gridPart.ibegin(entity); 
	it != endit; ++it)
      
      if(it.boundary())
	{
	  // if boundary, get cog of intersection and check for Neumann
	  IntersectionQuadratureType 
	    iquad(it,0,IntersectionQuadratureType::INSIDE);
	  // check, whether real cog integration is performed
	  assert(iquad.nop()==1);
	  if (this->model().boundaryType(entity,iquad,0) == TraitsType::Robin)
	    {
	      // if Neumann-boundary, then integrate over intersect:
	      // get quadrature and function space
	      IntersectionQuadratureType 
		iquad(it,TraitsType::quadDegree, 
                      IntersectionQuadratureType::INSIDE);
	      
	      // get local basis
	      const BaseFunctionSetType &baseSet = 
            fspace.baseFunctionSet( entity );
	      
	      int numBaseFunctions =  baseSet.numBaseFunctions();             
	      	      
	      for ( int pt=0; pt < iquad.nop(); pt++ ) 
		{  
		  // the following seems wrong...
		  const double vol = 
                      it.intersectionGlobal().
                      integrationElement(iquad.localPoint(pt));                 
//		  const double vol = 
//		    iquad.geometry().integrationElement(
//		      iquad.localPoint(pt));
		  // entity.geometry().integrationElement(iquad.point(pt));
		  
		  double fact = iquad.weight( pt ) * vol * coef;
		  
		  // get normal
		  // TraitsType::DomainType outerNormal = 
		  // it.unitOuterNormal(iquad.localPoint(pt));
		  
		  for(int i=0; i<numBaseFunctions; i++) 
                    {
		      RangeType phi;
		      baseSet.evaluate(i,iquad,pt,phi);  
		      
   		      RangeType ret;
		      this->model().robinValues(entity, iquad, pt, ret);
                      ret[0] *= fact;
		      double incr =  ret[0] * phi[0];
		      elRhs[i] += incr;
                    } // end loop over basefunctions           
		} // end loop over quadraturepoints
	    } // end of if isdirichlet-intersection
	} // end of isboundary
  };

  /*======================================================================*/
/*! 
 *   addGeneralizedNeumannElementRhs: adds a multiple of a Robin term contribution
 *                     to the rhs 
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
/*======================================================================*/  
    
    template <class EntityImp, class ElementRhsImp>
    void addGeneralizedNeumannElementRhs(EntityImp &entity, 
                            ElementRhsImp &elRhs, 
                            double coef = 1.0) // const
          {
    // for all intersections check whether boundary
    DiscreteFunctionSpaceType& fspace = 
      this->model().discreteFunctionSpace();          
    GridPartType & gridPart = fspace.gridPart();
    
    IntersectionIteratorType 
      endit = gridPart.iend(entity);
    for(IntersectionIteratorType it = 
	  gridPart.ibegin(entity); 
	it != endit; ++it)
      
      if(it.boundary())
	{
	  // if boundary, get cog of intersection and check for Neumann
	  IntersectionQuadratureType 
	    iquad(it,0,IntersectionQuadratureType::INSIDE);
	  // check, whether real cog integration is performed
	  assert(iquad.nop()==1);
	  if (this->model().boundaryType(entity,iquad,0) == TraitsType::GeneralizedNeumann)
	    {
	      // if Neumann-boundary, then integrate over intersect:
	      // get quadrature and function space
	      IntersectionQuadratureType 
		iquad(it,TraitsType::quadDegree, 
                      IntersectionQuadratureType::INSIDE);
	      
	      // get local basis
	      const BaseFunctionSetType & baseSet = 
		fspace.baseFunctionSet(entity);
	      
	      int numBaseFunctions =  baseSet.numBaseFunctions();             
	      	      
	      for ( int pt=0; pt < iquad.nop(); pt++ ) 
		{  
		  // the following seems wrong...
		  const double vol = 
                      it.intersectionGlobal().
                      integrationElement(iquad.localPoint(pt));                 
//		  const double vol = 
//		    iquad.geometry().integrationElement(
//		      iquad.localPoint(pt));
		  // entity.geometry().integrationElement(iquad.point(pt));
		  
		  double fact = iquad.weight( pt ) * vol * coef;
		  
		  // get normal
		  // TraitsType::DomainType outerNormal = 
		  // it.unitOuterNormal(iquad.localPoint(pt));
		  
		  for(int i=0; i<numBaseFunctions; i++) 
                    {
		      RangeType phi;
		      baseSet.eval(i,iquad,pt,phi);  
		      
   		      RangeType ret;
		      this->model().generalizedNeumannValues(entity, iquad, pt, it,ret);
                      ret[0] *= fact;
		      double incr =  ret[0] * phi[0];
		      elRhs[i] += incr;
                    } // end loop over basefunctions           
		} // end loop over quadraturepoints
	    } // end of if isdirichlet-intersection
	} // end of isboundary
  }; //end of addGeneralizedNeumannElementRhs

  
  private:
    //! reference to the underlying model specified during construction
    ModelType& model_;
    //! verbosity flag;
    int verbose_;
  }; // end class DefaultElementRhsIntegrator
     

/*======================================================================*/
/*!
 *  \class RhsAssembler
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
/*======================================================================*/
     
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
          DiscreteFunctionSpaceType& fspace 
              =  elRhsInt_.model().discreteFunctionSpace();
          
          IteratorType it    = fspace.begin(); 
          IteratorType endit = fspace.end(); 
          
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
      DiscreteFunctionSpaceType &fspace = model.discreteFunctionSpace();
      
      GridPartType & gridPart = fspace.gridPart();
                    
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
          
          // get center of gravity of intersection and check for dirichlet
          IntersectionQuadratureType
            iquad( nit, 0, IntersectionQuadratureType :: INSIDE );
          assert( iquad.nop() == 1 );
          
          if( model.boundaryType( entity, iquad, 0 ) != TraitsType :: Dirichlet )
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
            model.dirichletValues( entity, lagrangePointSet, entityDofNumber, phi );
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

