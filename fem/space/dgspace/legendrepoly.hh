#ifndef DUNE_LEGENDREPOLY_HH
#define DUNE_LEGENDREPOLY_HH

#include <cassert>

namespace Dune
{

class LegendrePoly
{
 public: 
  LegendrePoly()
  { 
  init(); }
  ~LegendrePoly()
  { 
  }
  void init()
  {
    factor_[0][0]= 1.0;
    factor_[0][1]= 0;
    factor_[0][2]= 0;
    factor_[0][3]= 0;
    factor_[0][4]= 0;
    factor_[0][5]= 0;
    factor_[0][6]= 0;
    factor_[0][7]= 0;
    factor_[0][8]= 0;
    factor_[0][9]= 0;
    factor_[0][10]= 0;
    weight[0] = 1.0;

    factor_[1][0]= -1.0;
    factor_[1][1]= 2.0;
    factor_[1][2]= 0;
    factor_[1][3]= 0;
    factor_[1][4]= 0;
    factor_[1][5]= 0;
    factor_[1][6]= 0;
    factor_[1][7]= 0;
    factor_[1][8]= 0;
    factor_[1][9]= 0;
    factor_[1][10]= 0;
    weight[1]=1.73205080756887729352744634151;

    //p=6x^2-6x+1 
    factor_[2][0]= 1.0;
    factor_[2][1]=- 6.0;
    factor_[2][2]= 6.0;
    factor_[2][3]= 0;
    factor_[2][4]= 0;
    factor_[2][5]= 0;
    factor_[2][6]= 0;
    factor_[2][7]= 0;
    factor_[2][8]= 0;
    factor_[2][9]= 0;
    factor_[2][10]= 0;
    weight[2]=2.23606797749978969640917366873;

    //p=20*x^3-30*x^2+12*x-1  
    factor_[3][0]= -1.0;
    factor_[3][1]= 12.0;
    factor_[3][2]= -30.0;
    factor_[3][3]= 20.0;
    factor_[3][4]= 0;
    factor_[3][5]= 0;
    factor_[3][6]= 0;
    factor_[3][7]= 0;
    factor_[3][8]= 0;
    factor_[3][9]= 0;
    factor_[3][10]= 0;
    weight[3]=2.64575131106459059050161575364;

    //p=70*x^4-140*x^3+90*x^2-20*x+1
    factor_[4][0]= 1.0;
    factor_[4][1]= -20.0;
    factor_[4][2]= 90.0;
    factor_[4][3]= -140.0;
    factor_[4][4]= 70.0;
    factor_[4][5]= 0;
    factor_[4][6]= 0;
    factor_[4][7]= 0;
    factor_[4][8]= 0;
    factor_[4][9]= 0;
    factor_[4][10]= 0;
    weight[4]=3.0;
    
    
    //p=252*x^5-630*x^4+560*x^3-210*x^2+30*x-1 
    factor_[5][0]= -1.0;
    factor_[5][1]= 30.0;
    factor_[5][2]= -210.0;
    factor_[5][3]= 560.0;
    factor_[5][4]= -630.0;
    factor_[5][5]= 252.0;
    factor_[5][6]= 0;
    factor_[5][7]= 0;
    factor_[5][8]= 0;
    factor_[5][9]= 0;
    factor_[5][10]= 0;
    weight[5]=3.31662479035539984911493273667;
      
    //p=924*x^6-2772*x^5+3150*x^4-1680*x^3+420*x^2-42*x+1
    factor_[6][0]= 1.0;
    factor_[6][1]= -42.0;
    factor_[6][2]= 420.0;
    factor_[6][3]= -1680.0;
    factor_[6][4]= 3150.0;
    factor_[6][5]= -2772.0;
    factor_[6][6]= 924.0;
    factor_[6][7]= 0;
    factor_[6][8]= 0;
    factor_[6][9]= 0;
    factor_[6][10]= 0;
    weight[6]=3.60555127546398929311922126747;
      
    //p=3432*x^7-12012*x^6+16632*x^5-11550*x^4+4200*x^3-756*x^2+56*x-1
    factor_[7][0]= -1.0;
    factor_[7][1]= 56.0;
    factor_[7][2]= -756.0;
    factor_[7][3]= 4200.0;
    factor_[7][4]= -11550.0;
    factor_[7][5]= 16632.0;
    factor_[7][6]= -12012.0;
    factor_[7][7]= 3432.0;
    factor_[7][8]= 0;
    factor_[7][9]= 0;
    factor_[7][10]= 0;
    weight[7]=3.87298334620741688517926539978;
    
    //p=1-72*x+1260*x^2-9240*x^3-72072*x^5+34650*x^4+84084*x^6-51480*x^7+12870*x^8
    factor_[8][0]= 1.0;
    factor_[8][1]= -72.00;
    factor_[8][2]= 1260.0;
    factor_[8][3]= -9240.0;
    factor_[8][4]= 34650.0;
    factor_[8][5]= -72072.0;
    factor_[8][6]= 84084.0;
    factor_[8][7]= -51480.0;
    factor_[8][8]= 12870.0;
    factor_[8][9]= 0;
    factor_[8][10]= 0;
    weight[8]=4.12310562561766054982140985597;
    
    //p=-1+90*x-1980*x^2+18480*x^3+252252*x^5-90090*x^4-420420*x^6+411840*x^7-218790*x^8+48620*x^9
    factor_[9][0]= -1.0;
    factor_[9][1]= 90.0;
    factor_[9][2]= -1980.0;
    factor_[9][3]= 18480.0;
    factor_[9][4]= -90090.0;
    factor_[9][5]= 252252.0;
    factor_[9][6]= -420420.0;
    factor_[9][7]= 411840.0;
    factor_[9][8]= -218790.0;
    factor_[9][9]=48620.0;
    factor_[9][10]= 0;
    weight[9]=4.35889894354067355223698198386;
    
    //p= 1-110*x+2970*x^2-34320*x^3-756756*x^5+210210*x^4+1681680*x^6-2333760*x^7+1969110*x^8-923780*x^9+184756*x^10
    factor_[10][0]= 1.0;
    factor_[10][1]= -110.0;
    factor_[10][2]= 2970.0;
    factor_[10][3]= -34920.0;
    factor_[10][4]= 210210.0;
    factor_[10][5]= -756756.0;
    factor_[10][6]= 1681680.0;
    factor_[10][7]= -2333760.0;
    factor_[10][8]= 1969110.0;
    factor_[10][9]= -923780.0;
    factor_[10][10]= 184756.0;
    weight[10]= 4.58257569495584000658804719373;

  }
  
  double eval(int num,double x) const
  {  
    assert(0<=num && num<maxPol);
    double phi=factor_[num][num];
    for(int i=num-1;i>=0;i--)
      {phi=phi*x+factor_[num][i];}
    return weight[num]*phi;
  }
  
  double eval1(int num,double x) const
  { 
    assert(0<=num && num<maxPol);
    double phi=0.0;
    if (num>=1){
      phi=factor_[num][num]*num;
      for(int i=num-1;i>=1;i--)
	{phi=phi*x+factor_[num][i]*i;}
    }
      
    return weight[num]*phi;
  }
     
  double eval2(int num,double x) const
  { 
    assert(0<=num && num<maxPol);
    double phi=0.0;    
    if (num>=2){
      phi=factor_[num][num]*num*(num-1);
      for(int i=num-1;i>=2;i--){
	phi=phi*x+factor_[num][i]*i*(i-1);
      }
    }
    return weight[num]*phi;
  }
  
 private:
  enum {maxPol=11};
  double factor_[maxPol][maxPol];
  double weight[maxPol];
};

} // end namespace Dune 

#endif
