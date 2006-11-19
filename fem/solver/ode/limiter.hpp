#ifndef LIMITER_HPP
#define LIMITER_HPP



namespace pardg
{

  // forward declaration
  class MeshBase;
  class Triang2d;
  class DG;
  class DG2d;
  
  
  class Limiter
  {
  public:
    virtual ~Limiter () {}
    virtual void operator()(double *u) = 0;
  };
  


  // limit to piecwise constant data if solution is not smooth
  // on i-th cell, mesh independent
  // smoothness[i] <= 1: smooth, smoothness[i] > 1: not smooth
  class PCLimiter : public Limiter
  {
  public:
    PCLimiter(int dim_system, int num_base_polys, 
	      MeshBase &mesh, DG &dg);
    
    virtual ~PCLimiter () {}
    // from Limiter
    virtual void operator()(double *U);

  private:
    const int dim_system, num_base_polys;
    MeshBase &mesh;
    DG &dg;
  };



  class SlopeLimiter2d : public Limiter
  {
  public:
    SlopeLimiter2d(int dim_system, int poly_order,
		   Triang2d &tr2d, DG2d &dg2d);
    virtual ~SlopeLimiter2d();
    
    // from Limiter
    virtual void operator()(double *U);
    
  protected:
    virtual void pre_op(double *U);
    virtual void post_op(double *U);
    
    const int dim_system, poly_order, num_base_polys;
    Triang2d &tr2d;
    DG2d &dg2d;  

  private:
    // num "integration points", 3 vertices of a triangle
    static const int num_pts = 3; 

    // number_of_base_polys for linear polynomials
    static const int num_lin_polys = 3; 
  
    double *op, *u, *m, *M;
  };


  
} // namespace pardg

#endif
