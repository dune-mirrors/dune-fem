#ifndef DUNE_ORTHONORMALBASE_HH
#define DUNE_ORTHONORMALBASE_HH

namespace OrthonormalBase_1D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_line(const int baseNum,DomainType xi );
  void grad_line(const int baseNum, DomainType xi,
                 JacobianRangeType grad );

}
  
namespace OrthonormalBase_2D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_triangle_2d (const int baseNum, DomainType xi );
  double eval_quadrilateral_2d (const int baseNum, DomainType xi );
  void grad_triangle_2d (const int baseNum, DomainType xi,
                          JacobianRangeType grad );
  void grad_quadrilateral_2d (const int baseNum, DomainType xi,
                               JacobianRangeType grad );
}


namespace OrthonormalBase_3D 
{
  typedef const double* DomainType;
  typedef double* JacobianRangeType;

  double eval_tetrahedron_3d (const int baseNum, DomainType xi );
  double eval_pyramid_3d (const int baseNum, DomainType xi );
  double eval_prism_3d (const int baseNum, DomainType xi );
  double eval_hexahedron_3d (const int baseNum, DomainType xi );

  void grad_tetrahedron_3d (const int baseNum, DomainType xi,
                             JacobianRangeType grad );
  void grad_pyramid_3d (const int baseNum, DomainType xi,
                         JacobianRangeType grad );
  void grad_prism_3d (const int baseNum, DomainType xi,
                       JacobianRangeType grad );
  void grad_hexahedron_3d (const int baseNum, DomainType xi,
                            JacobianRangeType grad );
}
#endif
