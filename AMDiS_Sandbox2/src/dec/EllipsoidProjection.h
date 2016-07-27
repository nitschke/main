#ifndef ELLIPSOIDPROJECT_H
#define ELLIPSOIDPROJECT_H

namespace AMDiS {

  class EllipsoidProject : public Projection
  {
  public:
    /// Constructor.
    EllipsoidProject(int id, 
		ProjectionType type,
    double cx,
    double cy,
    double cz) 
      : Projection(id, type),
	proj(cx, cy, cz), jproj(cx, cy, cz)
    {}

    /// Destructor.
    virtual ~EllipsoidProject() {}

    /// Implementation of Projection::project();
    void project(WorldVector<double> &x) 
    {
      x = proj(x);
    }

     AbstractFunction<WorldVector<double>, WorldVector<double> >* getProjection() {return &proj;}

     AbstractFunction<WorldMatrix<double>, WorldVector<double> >* getJProjection() {return &jproj;}

  protected:
  
  class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {
    
      public:
      Proj(double cx_, double cy_, double cz_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(), cx(cx_), cy(cy_), cz(cz_) {}
    
      /// Implementation of AbstractFunction::operator().
      WorldVector<double> operator()(const WorldVector<double>& x) const 
      {
        WorldVector<double> xS(x);
        //scale:
        xS[0] /= cx;
        xS[1] /= cy;
        xS[2] /= cz;
        return x / std::sqrt(xS * xS);
      }

      private:
        // axis scales
        double cx;
        double cy;
        double cz;
   };

  //jacobi matrix of projection
  class JProj : public AbstractFunction<WorldMatrix<double>, WorldVector<double> > {
  
    public:
    JProj(double cx_, double cy_, double cz_) : AbstractFunction<WorldMatrix<double>, WorldVector<double> >(), cx(cx_), cy(cy_), cz(cz_) {}
  
    /// Implementation of AbstractFunction::operator().
    WorldMatrix<double> operator()(const WorldVector<double>& x) const 
    {
      WorldVector<double> xS(x); // LE*x
      //scale:
      xS[0] /= cx;
      xS[1] /= cy;
      xS[2] /= cz;
      double normx = std::sqrt(xS * xS); // ||x||_{LE^2}
      WorldMatrix<double> J;
      xS[0] /= cx;
      xS[1] /= cy;
      xS[2] /= cz; // => xS = LE^2*x
      J.vecProduct(x, xS); // xX(LE^2*x)
      J *= - 1.0 / normx / normx; // - xXx/||x||^2
      for (int i = 0; i < 3; i++) J[i][i] += 1.0; // I - xXx/||x||^2
      return J / normx ;  // (I - xXx/||x||^2) / ||x||
    }

     private:
        // axis scales
        double cx;
        double cy;
        double cz;
  };


    Proj proj;
    JProj jproj;
  
  };

}

#endif

