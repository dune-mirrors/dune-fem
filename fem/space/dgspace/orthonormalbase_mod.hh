#ifndef DUNE_ORTHONORMALBASE_HH
#define DUNE_ORTHONORMALBASE_HH

namespace OrthonormalBase_1D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_line(int i,DomainType xi );
  void grad_line(int i, DomainType xi,
                 JacobianRangeType grad );

}
  
namespace OrthonormalBase_2D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_triangle_2d (int i, DomainType xi );
  double eval_quadrilateral_2d (int i, DomainType xi );
  void grad_triangle_2d (int i, DomainType xi,
                          JacobianRangeType grad );
  void grad_quadrilateral_2d (int i, DomainType xi,
                               JacobianRangeType grad );
}


namespace OrthonormalBase_3D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_tetrahedron_3d (int i, DomainType xi );
  double eval_pyramid_3d (int i, DomainType xi );
  double eval_prism_3d (int i, DomainType xi );
  double eval_hexahedron_3d (int i, DomainType xi );

  void grad_tetrahedron_3d (int i, DomainType xi,
                             JacobianRangeType grad );
  void grad_pyramid_3d (int i, DomainType xi,
                         JacobianRangeType grad );
  void grad_prism_3d (int i, DomainType xi,
                       JacobianRangeType grad );
  void grad_hexahedron_3d (int i, DomainType xi,
                            JacobianRangeType grad );
}
#endif
