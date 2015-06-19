#include <iostream>
#include <fstream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


struct double_per
{
  typedef double_per self;
  static constexpr double period = 2.0*M_PI;
  
  double_per(double val = 0.0) : val(mod(val)) {}
  
  operator double() const { return val; }
  
  inline double value() const { return val; }
  
  self& operator=(double const& rhs) { val = mod(rhs); return *this; }
  self& operator*=(double const& rhs) { val = mod(val*rhs); return *this; }
  self& operator/=(double const& rhs) { val = mod(val/rhs); return *this; }
  self& operator+=(double const& rhs) { val = mod(val+rhs); return *this; }
  self& operator-=(double const& rhs) { val = mod(val-rhs); return *this; }
  
  self& operator=(double_per const& rhs) { val = mod(rhs.value()); return *this; }
  self& operator*=(double_per const& rhs) { val = mod(val*rhs.value()); return *this; }
  self& operator/=(double_per const& rhs) { val = mod(val/rhs.value()); return *this; }
  self& operator+=(double_per const& rhs) { val = mod(val+rhs.value()); return *this; }
  self& operator-=(double_per const& rhs) { val = mod(val-rhs.value()); return *this; }
  
private:
  inline double mod(double val) const
  {
    double ang = fmod(val, period);
    if (ang > M_PI) return ang - period;
    if (-M_PI > ang) return ang + period;
    return ang;
  }
  
  double val;
};

inline double_per operator*(double_per lhs, double rhs) { return lhs*=rhs; }
inline double_per operator/(double_per lhs, double rhs) { return lhs/=rhs; }
inline double_per operator+(double_per lhs, double rhs) { return lhs+=rhs; }
inline double_per operator-(double_per lhs, double rhs) { return lhs-=rhs; }
  
using namespace std;
using namespace mtl;
using namespace itl;

int main() 
{
    ofstream csvout;
    csvout.open("output.csv", ios::out);
    csvout << "t,rhs,phi,phiExact" << endl;
    csvout << "0.0,0.0,0.0,0.0" << endl;

    const int N = 10000;
    double hInv = (1.0 + (double)(N));
    double h = 1. / hInv;
    double c = hInv * hInv;
    typedef compressed2D<double>  matrix_type;

    // Set up a matrix: (1/h^2)*[1 , -2 , 1] stencel
    matrix_type                   A(N, N);
    {
      mat::inserter<matrix_type> ins(A);
      for (int i = 0; i < N; ++i) {
        ins[i][i] << -2.0 * c;
        if (i > 0) ins[i][i-1] << c;
        if (i < N-1) ins[i][i+1] << c;
      }
    }

    // Create an ILU(0) preconditioner
    pc::ilu_0<matrix_type>        L(A);
    //pc::identity<matrix_type>        L(A);
    pc::identity<matrix_type>        R(A);
    
    dense_vector<double_per>          x(N, 0.0), b(N), xE(N), time(N);
    // solution: 2Pi(4t^3 - 6t^2 + 3t) mod (-Pi,Pi]
    // rhs(second derivation): 24Pi(2t-1) mod (-Pi,Pi]
    // t in [0,1]
    for (int i = 0; i < N; ++i) {
      double t = (double)(i+1) * h;
      time[i] = t;
      b[i] = 24.0 * M_PI*(2.0 * t - 1.0);
      xE[i] = 2.0 * M_PI *(4*t*t*t - 6*t*t+ 3*t);
    }

    
    // Termination criterion: r < 1e-6 * b or N iterations
    noisy_iteration<double_per>       iter(b, 500, 1.e-6);
    
    // Solve Ax == b with left preconditioner P
    bicgstab(A, x, b, L, iter);

    
    for (int i = 0; i < N; ++i) {
      csvout << time[i] << "," << b[i] << "," << x[i] << "," << xE[i] << endl;
    }
    csvout << "1.0,0.0,0.0,0.0" << endl;

    //x = xE;
    cout << "||A*x - b||_2 = " << two_norm(A*x - b) << endl;

    return 0;
}
