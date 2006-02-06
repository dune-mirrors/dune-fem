#ifndef DUNE_LEGENDREPOLY_HH
#define DUNE_LEGENDREPOLY_HH

namespace Dune {

class LegendrePoly
{
 public: 
  
  LegendrePoly():
    num_(0)
  { factor_[0]= 1.0;
      factor_[1]= 0;
      factor_[2]= 0;
      factor_[3]= 0;
      factor_[4]= 0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;}

   LegendrePoly(int k) :
    num_(k)
      
  {
    switch(k){
    case 0://p=1 
      factor_[0]= 1.0;
      factor_[1]= 0;
      factor_[2]= 0;
      factor_[3]= 0;
      factor_[4]= 0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 1://p=2x-1 
      factor_[0]= -1.0;
      factor_[1]= 2.0;
      factor_[2]= 0;
      factor_[3]= 0;
      factor_[4]= 0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 2://p=6x^2-6x+1 
      factor_[0]= 1.0;
      factor_[1]=- 6.0;
      factor_[2]= 6.0;
      factor_[3]= 0;
      factor_[4]= 0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;		 
    case 3://p=20*x^3-30*x^2+12*x-1  
      factor_[0]= -1.0;
      factor_[1]= 12.0;
      factor_[2]= -30.0;
      factor_[3]= 20.0;
      factor_[4]= 0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;		 
    case 4://p=70*x^4-140*x^3+90*x^2-20*x+1
      factor_[0]= 1.0;
      factor_[1]= -20.0;
      factor_[2]= 90.0;
      factor_[3]= -140.0;
      factor_[4]= 70.0;
      factor_[5]= 0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 5://p=252*x^5-630*x^4+560*x^3-210*x^2+30*x-1 
      factor_[0]= -1.0;
      factor_[1]= 30.0;
      factor_[2]= -210.0;
      factor_[3]= 560.0;
      factor_[4]= -630.0;
      factor_[5]= 252.0;
      factor_[6]= 0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
       
      break;
    case 6://p=924*x^6-2772*x^5+3150*x^4-1680*x^3+420*x^2-42*x+1
      factor_[0]= 1.0;
      factor_[1]= -42.0;
      factor_[2]= 420.0;
      factor_[3]= -1680.0;
      factor_[4]= 3150.0;
      factor_[5]= -2772.0;
      factor_[6]= 924.0;
      factor_[7]= 0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 7://p=3432*x^7-12012*x^6+16632*x^5-11550*x^4+4200*x^3-756*x^2+56*x-1
      factor_[0]= -1.0;
      factor_[1]= 56.0;
      factor_[2]= -756.0;
      factor_[3]= 4200.0;
      factor_[4]= -11550.0;
      factor_[5]= 16632.0;
      factor_[6]= -12012.0;
      factor_[7]= 3432.0;
      factor_[8]= 0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 8://p=1-72*x+1260*x^2-9240*x^3-72072*x^5+34650*x^4+84084*x^6-51480*x^7+12870*x^8
      factor_[0]= 1.0;
      factor_[1]= -72.00;
      factor_[2]= 1260.0;
      factor_[3]= -9240.0;
      factor_[4]= 34650.0;
      factor_[5]= -72072.0;
      factor_[6]= 84084.0;
      factor_[7]= -51480.0;
      factor_[8]= 12870.0;
      factor_[9]= 0;
      factor_[10]= 0;
      break;
    case 9://p=-1+90*x-1980*x^2+18480*x^3+252252*x^5-90090*x^4-420420*x^6+411840*x^7-218790*x^8+48620*x^9
      factor_[0]= -1.0;
      factor_[1]= 90.0;
      factor_[2]= -1980.0;
      factor_[3]= 18480.0;
      factor_[4]= -90090.0;
      factor_[5]= 252252.0;
      factor_[6]= -420420.0;
      factor_[7]= 411840.0;
      factor_[8]= -218790.0;
      factor_[9]=48620.0;
      factor_[10]= 0;
      break;
    case 10://p= 1-110*x+2970*x^2-34320*x^3-756756*x^5+210210*x^4+1681680*x^6-2333760*x^7+1969110*x^8-923780*x^9+184756*x^10
      factor_[0]= 1.0;
      factor_[1]= -110.0;
      factor_[2]= 2970.0;
      factor_[3]= -34920.0;
      factor_[4]= 210210.0;
      factor_[5]= -756756.0;
      factor_[6]= 1681680.0;
      factor_[7]= -2333760.0;
      factor_[8]= 1969110.0;
      factor_[9]= -923780.0;
      factor_[10]= 184756.0;
      break;
    }
  }
  
  double eval(double x) const
  {  
    double phi=factor_[num_];
    for(int i=num_-1;i>=0;i--)
      {phi=phi*x+factor_[i];}
    return phi;
  }
  
  double eval1(double  x) const
  { 
    double phi=0.0;
    if (num_>=1){
      phi=factor_[num_]*num_;
      for(int i=num_-1;i>=1;i--)
	{phi=phi*x+factor_[i]*i;}
    }
      
    return phi;
  }
     
  double eval2(double x) const
  { 
    double phi=0.0;    
    if (num_>=2){
      phi=factor_[num_]*num_*(num_-1);
      for(int i=num_-1;i>=2;i--){
	phi=phi*x+factor_[i]*i*(i-1);
      }
    }
    return phi;
  }
  
 private:
  double factor_[11];
  int num_;
};

} // end namespace Dune 

#endif
