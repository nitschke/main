#include "AMDiS.h"

using namespace std;
using namespace AMDiS;
class Quadratic {

  public:
    Quadratic(const FiniteElemSpace *finiteElemSpace);

    DOFVector<WorldVector<double> > getNormals();

    void writeVTK(string name);

  private: 

    FiniteElemSpace *feSpace;
    mtl::dense_vector<double> coeffs;

    class Phi : public AbstractFunction<double, WorldVector<double> >
    {
    public:
      Phi(mtl::dense_vector<double> coeffs_) : AbstractFunction<double, WorldVector<double> >(1), coeffs(coeffs_) {}
    
      double operator()(const WorldVector<double>& coords) const 
      {
        double x = coords[0];
        double y = coords[1];
        double z = coords[2];
        return 0.5 * (coeffs[0] * x * x + coeffs[1] * y * y + coeffs[2] * z * z)
                + coeffs[3] * y * z + coeffs[4] * x * z + coeffs[5] * x * y
                + coeffs[6] * x + coeffs[7] * y + coeffs[8] * z - 1.0;
      }

      private:
        mtl::dense_vector<double> coeffs;
    };

    class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
    {
    public:
      GradPhi(mtl::dense_vector<double> coeffs_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(1), coeffs(coeffs_) {}
    
      WorldVector<double> operator()(const WorldVector<double>& x) const 
      {
        WorldVector<double> rval(x);
        for (int k = 0; k < 3; k++){
          int k1 = (k+1)%3;
          int k2 = (k+2)%3;
          rval[k] = coeffs[k] * x[k] + coeffs[6-k-k1] * x[k1] + coeffs[6-k-k2] * x[k2] + coeffs[6+k];
        }
        return rval;
      }

      private:
        mtl::dense_vector<double> coeffs;
    };



};


