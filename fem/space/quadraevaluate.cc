#include"legendrepoly.hh"

template <class FunctionSpaceType>
double DGBaseFunctionWrapper<FunctionSpaceType>::
eval_quadrilateral_2d_l ( int PolOrd,int i, const DomainType & xi ) const
{
  double x,y;
  LegendrePoly px , py;
  x=xi[0]; y=xi[1];

  assert(i <=(PolOrd+1)*(PolOrd+1)-1);
  px=LegendrePoly((i-i%(PolOrd+1))/(PolOrd+1));
  py=LegendrePoly(i%(PolOrd+1));

  return(px.eval(x)*py.eval(y));
}


/*template <class FunctionSpaceType>
double DGBaseFunctionWrapper<FunctionSpaceType>::
eval_quadrilateral_2d_l ( int i, const DomainType & xi ) const
{
  return( eval_quadrilateral_2d(1,i, xi ));
  }*/



template <class FunctionSpaceType>
void DGBaseFunctionWrapper<FunctionSpaceType>::
grad_quadrilateral_2d_l (int PolOrd, int i, const DomainType & xi, JacobianRangeType & grad ) const
{
	double x, y;
	LegendrePoly px , py;
	x=xi[0]; y=xi[1];

	assert(i <=(PolOrd+1)*(PolOrd+1)-1);
	px=LegendrePoly((i-i%(PolOrd+1))/(PolOrd+1));
	py=LegendrePoly(i%(PolOrd+1));

	grad[0][0]=px.eval1(x)*py.eval(y);
	grad[0][1]=px.eval(x)*py.eval1(y);
	return;
}  
 





 template <class FunctionSpaceType>
double DGBaseFunctionWrapper<FunctionSpaceType>::
 eval_hexahedron_3d_l (int PolOrd, int i, const DomainType & xi ) const
{
  double x, y, z;
  LegendrePoly px, py, pz;
  x=xi[0]; y=xi[1]; z=xi[2];
  int& baseNum = i;
  assert(baseNum <=(PolOrd+1)*(PolOrd+1)*(PolOrd+1)-1);
  px=LegendrePoly((baseNum-baseNum%11-baseNum%((PolOrd+1)*(PolOrd+1)))/((PolOrd+1)*(PolOrd+1)));
  py=LegendrePoly((baseNum-baseNum%11)%121);
  pz=LegendrePoly(baseNum%(PolOrd+1));
       

  return(px.eval(x)*py.eval(y)*pz.eval(z));
}

template <class FunctionSpaceType>
void DGBaseFunctionWrapper<FunctionSpaceType>::
grad_hexahedron_3d_l (int PolOrd, int i, const DomainType & xi, JacobianRangeType & grad ) const
{
  double x, y, z;
  LegendrePoly px , py, pz;
  x=xi[0]; y=xi[1]; z=xi[2];
  int& baseNum = i;
  assert(baseNum <=(PolOrd+1)*(PolOrd+1)*(PolOrd+1)-1);
  px=LegendrePoly((baseNum-baseNum%11-baseNum%((PolOrd+1)*(PolOrd+1)))/((PolOrd+1)*(PolOrd+1)));
  py=LegendrePoly((baseNum-baseNum%11)%121);
  pz=LegendrePoly(baseNum%(PolOrd+1));

  grad[0][0]=px.eval1(x)*py.eval(y)*pz.eval(z);
  grad[0][1]=px.eval(x)*py.eval1(y)*pz.eval(z);
  grad[0][2]=px.eval(x)*py.eval(y)*pz.eval1(z);
  return;
}  

