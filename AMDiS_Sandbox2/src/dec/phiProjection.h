#include "AMDiS.h"

namespace AMDiS {

class PhiProject : public Projection {

public:

  PhiProject(int id, 
             ProjectionType type, 
             AbstractFunction<double, WorldVector<double> > *phi_,
             AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi_,
             double eps_ = 1.0e-6, int nMax_ = 100) : Projection(id, type), phi(phi_), gradPhi(gradPhi_), eps(eps_), nMax(nMax_), proj(this), jproj(this) {}

  
  void project(WorldVector<double> &x);

  void project(const WorldVector<double> &x, WorldVector<double> &v);

  WorldVector<double> getNormal(const WorldVector<double> &x);


  AbstractFunction<WorldVector<double>, WorldVector<double> >* getProjection() {return &proj;}

  AbstractFunction<WorldMatrix<double>, WorldVector<double> >* getJProjection(double dh_ = 1.E-6) {
    jproj.setdh(dh_);
    return &jproj;
  }


protected:

  class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {
    
      public:
      Proj(PhiProject *phiProject) : AbstractFunction<WorldVector<double>, WorldVector<double> >(), pp(phiProject) {}
    
      /// Implementation of AbstractFunction::operator().
      WorldVector<double> operator()(const WorldVector<double>& x) const 
      {
        WorldVector<double> xx = x;
        pp->project(xx);
        return xx;
      }

      private:
        PhiProject *pp;
   };

  //jacobi matrix of projection
  class JProj : public AbstractFunction<WorldMatrix<double>, WorldVector<double> > {
  
    public:
    JProj(PhiProject *phiProject) : AbstractFunction<WorldMatrix<double>, WorldVector<double> >(), pp(phiProject), dh(1.E-6) {}

    void setdh(double dh_) {
      dh = dh_;
    }
  
    /// Implementation of AbstractFunction::operator().
    WorldMatrix<double> operator()(const WorldVector<double>& x) const 
    {
      WorldMatrix<double> J;
      J.set(0.0);

      //WorldVector<double> x0 = x;
      //pp->project(x0);

      for (int i = 0; i < 3; ++i) {
        WorldVector<double> x1 = x;
        WorldVector<double> xM1 = x;
        x1[i] += dh;
        xM1[i] -= dh;
        //x.print();
        //std::cout << std::endl;
        //x1.print();
        //xM1.print();
        //std::cout << std::endl;
        pp->project(x1);
        pp->project(xM1);
        //x1.print();
        //xM1.print();
        //std::cout << std::endl;
        WorldVector<double> Ji = 0.5 * (1./dh) * (x1 - xM1);
        for (int j = 0; j < 3; ++j) {
          J[i][j] = Ji[j];
        }
      }

      //J.print();
      //WAIT;

      return J;
    }
    private:
    
    PhiProject *pp;

    double dh;
  };

  Proj proj;
  JProj jproj;

  double eps;
  int nMax;
  
  AbstractFunction<double, WorldVector<double> > *phi;
  AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi;

};

}
