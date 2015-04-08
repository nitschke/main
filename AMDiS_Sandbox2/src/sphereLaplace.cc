#include "AMDiS.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"
#include "io/VtkVectorWriter.h"
#include "io/ElementFileWriter.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// 1-form -> Rot(z)
class Alpha : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Alpha() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = -x[1];
    conBasis2[1] = x[0];
    conBasis2[2] = 0.0;
    return (con2 * conBasis2) * vec;
  }
};

// <alpha,[p,q]>
class Alpha_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Alpha_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return 2.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};


// Laplace-Beltrami of alpha
class LbAlpha : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbAlpha() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = -x[1];
    conBasis2[1] = x[0];
    conBasis2[2] = 0.0;
    return -2.0 * (con2 * conBasis2) * vec;
  }
};

// <lbalpha,[p,q]>
class LbAlpha_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbAlpha_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return -4.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};


// exact 1-form : gamma = d(xyz)
class DXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> gradF;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      gradF[i] = x[ii] * x[iii] * (1.0 - 3.0*x[i]*x[i]);
    }
    return  gradF * vec;
  }
};

// <dxyz,[p,q]>
class DXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

// <Lcbdxyz,[p,q]>
class LcbDXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LcbDXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  -12.*(q[0]*q[1]*q[2] - p[0]*p[1]*p[2]);
  }
};

// rot(x*y*z)
class RotXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> rot;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      rot[i] = x[i] * (x[ii]*x[ii] - x[iii]*x[iii]);
    }
    return  rot * vec;
  }
};

// <Rot(xyz),[p,q]>
class RotXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return 0.25*pow(1. + p1*q1 + p2*q2 + p3*q3,-2.)*(12.*(-1.*p3*q2 + p2*q3)*
      atan(pow((-1. + p1*q1 + p2*q2 + p3*q3)/(-1. - 1.*p1*q1 - 1.*p2*q2 - 1.*p3*q3),0.5))*
      (q2*q3*pow(p2,4.)*(-3. + 4.*pow(q2,2.) + 3.*pow(q3,2.)) + 
        q2*(p1*p3*q1*(-1.*pow(q2,2.) - 3.*pow(q3,2.) + pow(p3,2.)*(-1. + pow(q2,2.) + 4.*pow(q3,2.))) + 
           q3*(-1. + pow(p3,2.))*(-1.*pow(q2,2.) - 1.*pow(q3,2.) + pow(p3,2.)*(-3. + 3.*pow(q2,2.) + 4.*pow(q3,2.)))) + 
        q2*pow(p2,2.)*(p1*p3*q1*(-3. + 4.*pow(q2,2.) + 7.*pow(q3,2.)) + 
           q3*(3. - 5.*pow(q2,2.) - 4.*pow(q3,2.) + pow(p3,2.)*(-8. + 9.*pow(q2,2.) + 9.*pow(q3,2.)))) + 
        pow(p2,3.)*(p1*q1*q3*(-1. + 4.*pow(q2,2.) + pow(q3,2.)) + 
           p3*(1. + 4.*pow(q2,4.) - 4.*pow(q3,2.) + pow(q2,2.)*(-5. + 9.*pow(q3,2.)) + 3.*pow(q3,4.))) + 
        p2*(-1.*p1*q1*q3*(3.*pow(q2,2.) + pow(q3,2.)) + p1*q1*q3*pow(p3,2.)*(-3. + 7.*pow(q2,2.) + 4.*pow(q3,2.)) + 
           p3*(-3.*pow(q2,4.) + pow(q2,2.)*(3. - 8.*pow(q3,2.)) + 3.*pow(q3,2.) - 3.*pow(q3,4.)) + 
           pow(p3,3.)*(1. + 3.*pow(q2,4.) - 5.*pow(q3,2.) + pow(q2,2.)*(-4. + 9.*pow(q3,2.)) + 4.*pow(q3,4.))))*
      pow((-1. - 1.*p1*q1 - 1.*p2*q2 - 1.*p3*q3)*pow(-1. + p1*q1 + p2*q2 + p3*q3,5),-0.5) + 
     2.*pow(-1. + p1*q1 + p2*q2 + p3*q3,-2.)*(p1*p3*q1*(q2 - 2.*q3)*q3*(q2 + 2.*q3)*(pow(q2,2.) + pow(q3,2.)) - 
        4.*p1*q1*q3*pow(p3,5.)*(-1. + pow(q2,2.) + 2.*pow(q3,2.)) + 
        q2*pow(p2,5.)*(4.*p1*q1*(-1. + 2.*pow(q2,2.) + pow(q3,2.)) + p3*q3*(-3. + 4.*pow(q2,2.) + 3.*pow(q3,2.))) + 
        p1*q1*q3*pow(p3,3.)*(-1.*pow(q2,4.) + pow(q2,2.)*(-1. + 8.*pow(q3,2.)) + 8.*pow(q3,4.)) + 
        q2*pow(p2,3.)*(p1*q1*(-8.*pow(q2,4.) + pow(p3,2.)*(-3. + 8.*pow(q2,2.) - 3.*pow(q3,2.)) + pow(q3,2.) - 8.*pow(q2,2.)*pow(q3,2.) + 
              pow(q3,4.)) + p3*q3*(2. + 5.*pow(p3,2.)*pow(q2,2.) - 4.*pow(q2,4.) - 5.*(pow(p3,2.) + pow(q2,2.))*pow(q3,2.) + 3.*pow(q3,4.))) + 
        p2*q2*(-1.*p1*q1*(-4.*pow(q2,4.) - 3.*pow(q2,2.)*pow(q3,2.) + pow(p3,4.)*(-1. + pow(q2,2.) + 4.*pow(q3,2.)) + 
              pow(p3,2.)*(pow(q2,2.) + 4.*pow(q2,4.) - 3.*pow(q2,2.)*pow(q3,2.) - 4.*pow(q3,4.)) + pow(q3,4.)) + 
           p3*q3*(pow(p3,4.)*(3. - 3.*pow(q2,2.) - 4.*pow(q3,2.)) + (q2 - 1.*q3)*(q2 + q3)*(-2. + 3.*pow(q2,2.) + 3.*pow(q3,2.)) + 
              pow(p3,2.)*(-2. - 3.*pow(q2,4.) + 5.*pow(q2,2.)*pow(q3,2.) + 4.*pow(q3,4.)))) + 
        pow(p3,2.)*(2.*pow(q2,6.) - 1.*(-3. + pow(q2,2.))*pow(q2,2.)*pow(q3,2.) + (5. - 11.*pow(q2,2.))*pow(q3,4.) - 8.*pow(q3,6.)) + 
        pow(p2,2.)*(8.*pow(q2,6.) + p1*q1*q3*pow(p3,3.)*(3. + 3.*pow(q2,2.) - 8.*pow(q3,2.)) + pow(q2,2.)*(-3. + pow(q3,2.))*pow(q3,2.) + 
           pow(q2,4.)*(-5. + 11.*pow(q3,2.)) + pow(p3,4.)*
            (-1. + 2.*pow(q2,4.) + 11.*pow(q3,2.) - 1.*pow(q2,2.)*(1. + 5.*pow(q3,2.)) - 12.*pow(q3,4.)) + 
           p1*p3*q1*q3*(-4.*pow(q2,4.) + pow(q3,2.) - 3.*pow(q2,2.)*pow(q3,2.) + 4.*pow(q3,4.)) - 
           1.*(q2 - 1.*q3)*(q2 + q3)*pow(p3,2.)*(-3. + 8.*pow(q2,4.) - 1.*pow(q3,2.) + pow(q2,2.)*(-1. + 13.*pow(q3,2.)) + 8.*pow(q3,4.)) - 
           2.*pow(q3,6.)) + pow(p2,4.)*(-8.*pow(q2,6.) - 12.*pow(q2,4.)*pow(q3,2.) + p1*p3*q1*q3*(-1. + 4.*pow(q2,2.) + pow(q3,2.)) + 
           pow(p3,2.)*(1. + 12.*pow(q2,4.) + pow(q3,2.) + pow(q2,2.)*(-11. + 5.*pow(q3,2.)) - 2.*pow(q3,4.)) - 1.*pow(q3,4.) - 
           1.*pow(q2,2.)*(-5. + pow(q3,2.) + 2.*pow(q3,4.)) + pow(q3,6.)) + 
        pow(p3,4.)*(pow(q2,4.) - 1.*pow(q2,6.) + (-5. + pow(q2,2.) + 2.*pow(q2,4.))*pow(q3,2.) + 12.*pow(q2,2.)*pow(q3,4.) + 8.*pow(q3,6.)) - 
        1.*pow(p3,6.)*(8.*(-1. + pow(q2,2.))*pow(q3,2.) + 8.*pow(q3,4.) + pow(-1. + pow(q2,2.),2.)) + 
        pow(p2,6.)*(8.*pow(q2,4.) + 8.*pow(q2,2.)*(-1. + pow(q3,2.)) + pow(-1. + pow(q3,2.),2.)) - 
        1.*(q2 - 1.*q3)*(q2 + q3)*pow(pow(q2,2.) + pow(q3,2.),2.)));
  }
};


// Laplace-Beltrami of rot(x*y*z) -> -12*rot(x*y*z)
class LbRotXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbRotXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> rot;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      rot[i] = x[i] * (x[ii]*x[ii] - x[iii]*x[iii]);
    }
    return  -12.0 * (rot * vec);
  }
};


// <LbRot(xyz),[p,q]>
class LbRotXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbRotXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return -3.*pow(1. + p1*q1 + p2*q2 + p3*q3,-2.)*(12.*(-1.*p3*q2 + p2*q3)*
      atan(pow((-1. + p1*q1 + p2*q2 + p3*q3)/(-1. - 1.*p1*q1 - 1.*p2*q2 - 1.*p3*q3),0.5))*
      (q2*q3*pow(p2,4.)*(-3. + 4.*pow(q2,2.) + 3.*pow(q3,2.)) + 
        q2*(p1*p3*q1*(-1.*pow(q2,2.) - 3.*pow(q3,2.) + pow(p3,2.)*(-1. + pow(q2,2.) + 4.*pow(q3,2.))) + 
           q3*(-1. + pow(p3,2.))*(-1.*pow(q2,2.) - 1.*pow(q3,2.) + pow(p3,2.)*(-3. + 3.*pow(q2,2.) + 4.*pow(q3,2.)))) + 
        q2*pow(p2,2.)*(p1*p3*q1*(-3. + 4.*pow(q2,2.) + 7.*pow(q3,2.)) + 
           q3*(3. - 5.*pow(q2,2.) - 4.*pow(q3,2.) + pow(p3,2.)*(-8. + 9.*pow(q2,2.) + 9.*pow(q3,2.)))) + 
        pow(p2,3.)*(p1*q1*q3*(-1. + 4.*pow(q2,2.) + pow(q3,2.)) + 
           p3*(1. + 4.*pow(q2,4.) - 4.*pow(q3,2.) + pow(q2,2.)*(-5. + 9.*pow(q3,2.)) + 3.*pow(q3,4.))) + 
        p2*(-1.*p1*q1*q3*(3.*pow(q2,2.) + pow(q3,2.)) + p1*q1*q3*pow(p3,2.)*(-3. + 7.*pow(q2,2.) + 4.*pow(q3,2.)) + 
           p3*(-3.*pow(q2,4.) + pow(q2,2.)*(3. - 8.*pow(q3,2.)) + 3.*pow(q3,2.) - 3.*pow(q3,4.)) + 
           pow(p3,3.)*(1. + 3.*pow(q2,4.) - 5.*pow(q3,2.) + pow(q2,2.)*(-4. + 9.*pow(q3,2.)) + 4.*pow(q3,4.))))*
      pow((-1. - 1.*p1*q1 - 1.*p2*q2 - 1.*p3*q3)*pow(-1. + p1*q1 + p2*q2 + p3*q3,5),-0.5) + 
     2.*pow(-1. + p1*q1 + p2*q2 + p3*q3,-2.)*(p1*p3*q1*(q2 - 2.*q3)*q3*(q2 + 2.*q3)*(pow(q2,2.) + pow(q3,2.)) - 
        4.*p1*q1*q3*pow(p3,5.)*(-1. + pow(q2,2.) + 2.*pow(q3,2.)) + 
        q2*pow(p2,5.)*(4.*p1*q1*(-1. + 2.*pow(q2,2.) + pow(q3,2.)) + p3*q3*(-3. + 4.*pow(q2,2.) + 3.*pow(q3,2.))) + 
        p1*q1*q3*pow(p3,3.)*(-1.*pow(q2,4.) + pow(q2,2.)*(-1. + 8.*pow(q3,2.)) + 8.*pow(q3,4.)) + 
        q2*pow(p2,3.)*(p1*q1*(-8.*pow(q2,4.) + pow(p3,2.)*(-3. + 8.*pow(q2,2.) - 3.*pow(q3,2.)) + pow(q3,2.) - 8.*pow(q2,2.)*pow(q3,2.) + 
              pow(q3,4.)) + p3*q3*(2. + 5.*pow(p3,2.)*pow(q2,2.) - 4.*pow(q2,4.) - 5.*(pow(p3,2.) + pow(q2,2.))*pow(q3,2.) + 3.*pow(q3,4.))) + 
        p2*q2*(-1.*p1*q1*(-4.*pow(q2,4.) - 3.*pow(q2,2.)*pow(q3,2.) + pow(p3,4.)*(-1. + pow(q2,2.) + 4.*pow(q3,2.)) + 
              pow(p3,2.)*(pow(q2,2.) + 4.*pow(q2,4.) - 3.*pow(q2,2.)*pow(q3,2.) - 4.*pow(q3,4.)) + pow(q3,4.)) + 
           p3*q3*(pow(p3,4.)*(3. - 3.*pow(q2,2.) - 4.*pow(q3,2.)) + (q2 - 1.*q3)*(q2 + q3)*(-2. + 3.*pow(q2,2.) + 3.*pow(q3,2.)) + 
              pow(p3,2.)*(-2. - 3.*pow(q2,4.) + 5.*pow(q2,2.)*pow(q3,2.) + 4.*pow(q3,4.)))) + 
        pow(p3,2.)*(2.*pow(q2,6.) - 1.*(-3. + pow(q2,2.))*pow(q2,2.)*pow(q3,2.) + (5. - 11.*pow(q2,2.))*pow(q3,4.) - 8.*pow(q3,6.)) + 
        pow(p2,2.)*(8.*pow(q2,6.) + p1*q1*q3*pow(p3,3.)*(3. + 3.*pow(q2,2.) - 8.*pow(q3,2.)) + pow(q2,2.)*(-3. + pow(q3,2.))*pow(q3,2.) + 
           pow(q2,4.)*(-5. + 11.*pow(q3,2.)) + pow(p3,4.)*
            (-1. + 2.*pow(q2,4.) + 11.*pow(q3,2.) - 1.*pow(q2,2.)*(1. + 5.*pow(q3,2.)) - 12.*pow(q3,4.)) + 
           p1*p3*q1*q3*(-4.*pow(q2,4.) + pow(q3,2.) - 3.*pow(q2,2.)*pow(q3,2.) + 4.*pow(q3,4.)) - 
           1.*(q2 - 1.*q3)*(q2 + q3)*pow(p3,2.)*(-3. + 8.*pow(q2,4.) - 1.*pow(q3,2.) + pow(q2,2.)*(-1. + 13.*pow(q3,2.)) + 8.*pow(q3,4.)) - 
           2.*pow(q3,6.)) + pow(p2,4.)*(-8.*pow(q2,6.) - 12.*pow(q2,4.)*pow(q3,2.) + p1*p3*q1*q3*(-1. + 4.*pow(q2,2.) + pow(q3,2.)) + 
           pow(p3,2.)*(1. + 12.*pow(q2,4.) + pow(q3,2.) + pow(q2,2.)*(-11. + 5.*pow(q3,2.)) - 2.*pow(q3,4.)) - 1.*pow(q3,4.) - 
           1.*pow(q2,2.)*(-5. + pow(q3,2.) + 2.*pow(q3,4.)) + pow(q3,6.)) + 
        pow(p3,4.)*(pow(q2,4.) - 1.*pow(q2,6.) + (-5. + pow(q2,2.) + 2.*pow(q2,4.))*pow(q3,2.) + 12.*pow(q2,2.)*pow(q3,4.) + 8.*pow(q3,6.)) - 
        1.*pow(p3,6.)*(8.*(-1. + pow(q2,2.))*pow(q3,2.) + 8.*pow(q3,4.) + pow(-1. + pow(q2,2.),2.)) + 
        pow(p2,6.)*(8.*pow(q2,4.) + 8.*pow(q2,2.)*(-1. + pow(q3,2.)) + pow(-1. + pow(q3,2.),2.)) - 
        1.*(q2 - 1.*q3)*(q2 + q3)*pow(pow(q2,2.) + pow(q3,2.),2.)));
  }
};


// e.g. Laplace-Beltrami of exact forms
class Zero : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Zero() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    return 0.0;
  }
};


//projection
class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

  public:
  Proj() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    return x / sqrt(x * x);
  }
};

//jacobi matrix of projection
class JProj : public AbstractFunction<WorldMatrix<double>, WorldVector<double> > {

  public:
  JProj() : AbstractFunction<WorldMatrix<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldMatrix<double> operator()(const WorldVector<double>& x) const 
  {
    double normx = sqrt(x * x); // ||x||
    WorldMatrix<double> J;
    J.vecProduct(x, x); // xXx
    J *= - 1.0 / normx / normx; // - xXx/||x||^2
    for (int i = 0; i < 3; i++) J[i][i] += 1.0; // I - xXx/||x||^2
    return J / normx ;  // (I - xXx/||x||^2) / ||x||
  }


};




// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  const EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());
  //cout << *edgeMesh << endl;

 
  //cout << endl << "********** Alpha *************" << endl;
  //cout <<         "*** Laplace-Beltrami ***" << endl;
  //DofEdgeVector alphad(edgeMesh, "alpha_GL4");
  ////alphad.interpolGL4(new Alpha(), new Proj(), new JProj());
  //alphad.interpolLinTrapz(new Alpha());
  ////alphad.set(new Alpha_d());
  //alphad.writeFile("output/alpha.vtu");

  //DOFVector< WorldVector<double> > alphaSharp = alphad.getSharpFaceAverage();
  //AMDiS::io::VtkVectorWriter::writeFile(alphaSharp, string("output/alphaSharp.vtu"));

  //DofEdgeVector lbAlpha = alphad.laplaceBeltrami();
  //DOFVector< WorldVector<double> > lbAlphaSharp = lbAlpha.getSharpFaceAverage();
  //AMDiS::io::VtkVectorWriter::writeFile(lbAlphaSharp, string("output/lbAlphaSharp.vtu"));

  //DofEdgeVector lbAlphaSol(edgeMesh, "lbAlpha");
  ////lbAlphaSol.interpolGL4(new LbAlpha(), new Proj(), new JProj());
  //lbAlphaSol.set(new LbAlpha_d());
  //double error = lbAlpha.error(lbAlphaSol);
  //cout << "Error l2: " << (error * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << error << endl;
  //cout << "ErrorMax: " << lbAlpha.errorMax(lbAlphaSol) << endl;

  //DofEdgeVector lbAlphaError = lbAlphaSol - lbAlpha;
  //lbAlphaError.writeFile("output/lbAlphaErrorEdges.vtu");
  //DOFVector< WorldVector<double> > lbAlphaErrorSharp = lbAlphaError.getSharpFaceAverage();
  //AMDiS::io::VtkVectorWriter::writeFile(lbAlphaErrorSharp, string("output/lbAlphaSharpError.vtu"));

  //cout <<         "*** Laplace-CoBeltrami ***" << endl;
  //DofEdgeVector lcbAlpha = alphad.laplaceCoBeltrami();
  //lcbAlpha.writeFile("output/lcbalpha.vtu");

  //DofEdgeVector lcbAlphaSol(edgeMesh, "lcbAlpha");
  //lcbAlphaSol.set(0.0);

  //error = lcbAlpha.error(lcbAlphaSol);
  //cout << "Error l2: " << (error * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << error << endl;
  //cout << "ErrorMax: " << lcbAlpha.errorMax(lcbAlphaSol) << endl;

  //cout <<         "*** Laplace-deRham ***" << endl;
  //DofEdgeVector ldrAlpha(edgeMesh, "ldrAlpha");
  //ldrAlpha.set(0.0);
  //ldrAlpha -= lbAlpha;
  //ldrAlpha -= lcbAlpha;
  //ldrAlpha.writeFile("output/ldralpha.vtu");
  //
  //DofEdgeVector ldrAlphaSol(edgeMesh, "ldrAlpha");
  //ldrAlphaSol.set(0.0);
  //ldrAlphaSol -= lbAlphaSol;

  //error = ldrAlpha.error(ldrAlphaSol);
  //cout << "Error l2: " << (error * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << error << endl;
  //cout << "ErrorMax: " << ldrAlpha.errorMax(ldrAlphaSol) << endl;






  //cout << endl << "********** Rot(x*y*z) *************" << endl;
  //DofEdgeVector rotxyz(edgeMesh, "rotxyz");
  //rotxyz.interpolGL4(new RotXYZ(), new Proj(), new JProj());
  ////rotxyz.interpolLinTrapz(new RotXYZ());
  ////rotxyz.set(new RotXYZ_d());
  //rotxyz.writeFile("output/rotxyz.vtu");
  //DOFVector< WorldVector<double> > rotxyzSharp = rotxyz.getSharpFaceAverage();
  //AMDiS::io::VtkVectorWriter::writeFile(rotxyzSharp, string("output/rotxyzSharp.vtu"));

  //cout <<         "*** Laplace-Beltrami ***" << endl;
  //DofEdgeVector lbrotxyz = rotxyz.laplaceBeltrami();
  //lbrotxyz.writeFile("output/lbrotxyz.vtu");

  //DofEdgeVector lbrotxyzSol(edgeMesh, "lbRotXYZSol");
  //lbrotxyzSol.interpolGL4(new LbRotXYZ(), new Proj(), new JProj());
  ////lbrotxyz.set(new LbRotXYZ_d());
  //lbrotxyzSol.writeFile("output/lbrotxyzSol.vtu");
  //double errorrotxyz = lbrotxyz.error(lbrotxyzSol);
  //cout << "Error l2: " << (errorrotxyz * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << errorrotxyz << endl;
  //cout << "ErrorMax: " << lbrotxyz.errorMax(lbrotxyzSol) << endl;

  //cout <<         "*** Laplace-CoBeltrami ***" << endl;
  //DofEdgeVector lcbrotxyz = rotxyz.laplaceCoBeltrami();
  //lcbrotxyz.writeFile("output/lcbrotxyz.vtu");

  //DofEdgeVector lcbrotxyzSol(edgeMesh, "lcbrotxyz");
  //lcbrotxyzSol.set(0.0);

  //errorrotxyz = lcbrotxyz.error(lcbrotxyzSol);
  //cout << "Error l2: " << (errorrotxyz * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << errorrotxyz << endl;
  //cout << "ErrorMax: " << lcbrotxyz.errorMax(lcbrotxyzSol) << endl;

  //cout <<         "*** Laplace-deRham ***" << endl;
  //DofEdgeVector ldrrotxyz(edgeMesh, "ldrrotxyz");
  //ldrrotxyz.set(0.0);
  //ldrrotxyz -= lbrotxyz;
  //ldrrotxyz -= lcbrotxyz;
  //ldrrotxyz.writeFile("output/ldrrotxyz.vtu");
  //
  //DofEdgeVector ldrrotxyzSol(edgeMesh, "ldrrotxyz");
  //ldrrotxyzSol.set(0.0);
  //ldrrotxyzSol -= lbrotxyzSol;

  //errorrotxyz = ldrrotxyz.error(ldrrotxyzSol);
  //cout << "Error l2: " << (errorrotxyz * edgeMesh->getVol()) << endl;
  //cout << "Error   : " << errorrotxyz << endl;
  //cout << "ErrorMax: " << ldrrotxyz.errorMax(ldrrotxyzSol) << endl;
  
  cout << endl << "********** d(x*y*z) *************" << endl;
  DofEdgeVector dxyz(edgeMesh, "dxyz");
  //dxyz.interpolGL4(new DXYZ(), new Proj(), new JProj());
  dxyz.interpolNC(new DXYZ(), 7, new Proj());
  //dxyz.interpolLinTrapz(new DXYZ());
  //dxyz.set(new DXYZ_d());
  dxyz.writeFile("output/dxyz.vtu");
  DOFVector< WorldVector<double> > dxyzSharp = dxyz.getSharpFaceAverage();
  AMDiS::io::VtkVectorWriter::writeFile(dxyzSharp, string("output/dxyzSharp.vtu"));

  cout <<         "*** Laplace-CoBeltrami ***" << endl;
  DofEdgeVector lbdxyz = dxyz.laplaceBeltrami();
  lbdxyz.writeFile("output/lbdxyz.vtu");

  DofEdgeVector lbdxyzSol(edgeMesh, "lbDXYZSol");
  lbdxyzSol.set(0.0);
  lbdxyzSol.writeFile("output/lbdxyzSol.vtu");
  double errordxyz = lbdxyz.error(lbdxyzSol);
  cout << "Error l2: " << (errordxyz * edgeMesh->getVol()) << endl;
  cout << "Error   : " << errordxyz << endl;
  cout << "ErrorMax: " << lbdxyz.errorMax(lbdxyzSol) << endl;

  cout <<         "*** Laplace-CoBeltrami ***" << endl;
  DofEdgeVector lcbdxyz = dxyz.laplaceCoBeltrami();
  lcbdxyz.writeFile("output/lcbdxyz.vtu");

  DofEdgeVector lcbdxyzSol(edgeMesh, "lcbdxyz");
  lcbdxyzSol.set(new LcbDXYZ_d());

  errordxyz = lcbdxyz.error(lcbdxyzSol);
  cout << "Error l2: " << (errordxyz * edgeMesh->getVol()) << endl;
  cout << "Error   : " << errordxyz << endl;
  cout << "ErrorMax: " << lcbdxyz.errorMax(lcbdxyzSol) << endl;

  cout <<         "*** Laplace-deRham ***" << endl;
  DofEdgeVector ldrdxyz(edgeMesh, "ldrdxyz");
  ldrdxyz.set(0.0);
  ldrdxyz -= lbdxyz;
  ldrdxyz -= lcbdxyz;
  ldrdxyz.writeFile("output/ldrdxyz.vtu");
  
  DofEdgeVector ldrdxyzSol(edgeMesh, "ldrdxyz");
  ldrdxyzSol.set(0.0);
  ldrdxyzSol -= lcbdxyzSol;

  errordxyz = ldrdxyz.error(ldrdxyzSol);
  cout << "Error l2: " << (errordxyz * edgeMesh->getVol()) << endl;
  cout << "Error   : " << errordxyz << endl;
  cout << "ErrorMax: " << ldrdxyz.errorMax(ldrdxyzSol) << endl;

  //map<int, std::vector<double> > alphaFaceSharp = alphadGL4.getSharpOnFaces();
  //AMDiS::io::ElementFileWriter::writeFile(alphaFaceSharp, sphere.getFeSpace()->getMesh(), "output/alphaFaceSharp");

  //// === create adapt info ===
  //AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  //// === create adapt ===
  //AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
	//				       &sphere,
	//				       adaptInfo);
  

  // ===== start adaption loop =====
  //adapt->adapt();

  //sphere.writeFiles(adaptInfo, true);
  


  //DOFVector<WorldVector<double> > WDV(sphere.getFeSpace(),"vector");
  //WDV.interpol(new Alpha());
  //AMDiS::io::VtkVectorWriter::writeFile(WDV, string("output/vector.vtu"));

  AMDiS::finalize();
}


