#ifndef SPHEREPROJECT_H
#define SPHEREPROJECT_H

namespace AMDiS {

  class SphereProject : public Projection
  {
  public:
    /// Constructor.
    SphereProject(int id, 
		ProjectionType type,
		double radius = 1) 
      : Projection(id, type),
	radius_(radius), proj(), jproj()
    {}

    /// Destructor.
    virtual ~SphereProject() {}

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


    Proj proj;
    JProj jproj;

    /// Radius of the ball.
    double radius_;
  };

}

#endif

