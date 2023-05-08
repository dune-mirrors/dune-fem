#ifndef DUNE_FEM_GAUSSPOINTS_HH
#define DUNE_FEM_GAUSSPOINTS_HH

#include <cassert>
#include <vector>

#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Fem
  {

    /*! \class QuadPtsBase
     *  \ingroup Quadrature
     *  \brief one-dimensional quadrature points and their weights
     *
     *  QuadPtsBase is an array of one-dimensional quadratures for the
     *  interval [0,1]. The index of a quadrature equals its number of quadrature
     *  points (so there is no 0-th quadrature).
     *
     *  \note This class implements the Singleton pattern
     */
    class QuadPtsBase
    {
    protected:
      std::vector< std::vector< double > > G; //[MAXP+1][MAXP]; // positions of Gauss points
      std::vector< std::vector< double > > W; //[MAXP+1][MAXP]; // weights associated with points
      std::vector< int > O;                   //[MAXP+1];       // order of the rule

      /*! \brief constructor initializing the points for all orders
       */
      QuadPtsBase ( const int maxp )
       : G( maxp+1, std::vector<double>(maxp,0.0)),
         W( maxp+1, std::vector<double>(maxp,0.0)),
         O( maxp+1, 0 )
      {}

    public:
      /*! \brief obtain the i-th point of the m-th quadrature
       *
       *  \param[in]  m  index of the quadrature
       *  \param[in]  i  number of the point within the quadrature (0 <= i < m)
       *
       *  \returns a double in [0,1] representing the i-th Gauss point
       */
      double point ( int m, int i ) const
      {
        assert(m > 0 && i < m);
        return G[m][i];
      }

      /*! \brief obtain the i-th weight of the m-th quadrature
       *
       *  \param[in]  m  index of the quadrature
       *  \param[in]  i  number of the weight within the quadrature (0 <= i < m)
       *
       *  \returns a double representing the weight i-th Gauss point
       */
      double weight ( int m, int i ) const
      {
        assert(m > 0 && i < m);
        return W[m][i];
      }

      /*! \brief obtain the order of the m-th quadrature
       *
       *  \param[in]  m  index of the quadrature
       *
       *  \returns a double representing the weight i-th Gauss point
       */
      int order ( int m ) const
      {
        return O[m];
      }

      /*! \brief a simple power method
       *
       *  \note This method does not use a template meta program
       *
       *  \param[in] y  base \f$y\f$ of the power
       *  \param[in] d  exponent \f$d\f$ of the power
       *
       *  \returns \f$y^d\f$
       */
      int power ( int y, int d ) const
      {
        int m = 1;
        for( int i = 0; i < d; ++i )
          m *= y;
        return m;
      }
    };


    /*! \class GaussPts
     *  \ingroup Quadrature
     *  \brief one-dimensional Gauss points and their weights
     *
     *  GaussPtr is an array of one-dimensional Gauss quadratures for the
     *  interval [0,1]. The index of a quadrature equals its number of quadrature
     *  points (so there is no 0-th quadrature).
     *
     *  \note This class implements the Singleton pattern
     */
    class GaussPts : public QuadPtsBase
    {
    public:
      //! number of available quadratures
      static const int MAXP=10;

      //! highest quadrature order within the array
      static const int highestOrder=19;

    protected:
      using QuadPtsBase :: G;
      using QuadPtsBase :: W;
      using QuadPtsBase :: O;

    public:
      /*! \brief constructor initializing the Gauss points for all orders
       */
      GaussPts ();
    };

    /*! \class Modified Newton-Cotes (for Lobatto <--> FV identification)
     *  \ingroup Quadrature
     *  \brief one-dimensional modified Newton-Cotes points and their weights
     *         The difference to Newton-Cotes is that the
     *         first and last interval is only h/2 instead of h in the original Newton-Cotes rule.
     *         This allows to compute integrals of DG functions on quadrature
     *         points that correspond to a finite volume submesh.
     *
     *  ModifiedNewtonCotes is an array of one-dimensional quadratures for the
     *  interval [0,1]. The index of a quadrature equals its number of quadrature
     *  points (so there is no 0-th quadrature).
     *
     *  \note This class implements the Singleton pattern
     */
    class ModifiedNewtonCotes : public QuadPtsBase
    {
    public:
      //! number of available quadratures
      static const int MAXP=10;

      //! highest quadrature order within the array
      static const int highestOrder=10;

    protected:
      using QuadPtsBase :: G;
      using QuadPtsBase :: W;
      using QuadPtsBase :: O;

    public:
      /*! \brief constructor initializing the points for all orders
       */
      ModifiedNewtonCotes ();
    };
  } // namespace Fem

} // namespace Dune

#include "gausspoints_implementation.hh"
#include "modifiednewtoncotes_implementation.hh"
#endif // #ifndef DUNE_FEM_GAUSSPOINTS_HH
