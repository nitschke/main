#include "DOFCondition.h"

void ProblemStatExt::setDOFCondition(int iFESpace, 
                         int iDOF, 
                         double val,
                         bool verbose) 
  {
    FUNCNAME("ProblemVecExt::setDOFCondition");

    MSG("upgrade Systemmatrix ...\n");

    int n= getSystemMatrix(0,0)->getUsedSize();
    int nM= getSystemMatrix()->getSize();

    mtl::compressed2D<double> additionalRow(n,n);
    {
      mtl::matrix::inserter< mtl::compressed2D<double> > additionalRowIns(additionalRow);
      additionalRowIns[iDOF][iDOF]= 1.0;
    }

    Matrix< mtl::compressed2D<double> > A(nM,nM);
    if (verbose) cout << "before upgrade:\n";
    for (int i= 0; i < nM ; i++)
      for (int j= 0; j < nM ; j++)
      {
        A[i][j]= getSystemMatrix(i,j)->getBaseMatrix();
        if (verbose)
        {
          getSystemMatrix(i,j)->calculateNnz();
          cout << "nnz per row (" << i << "," << j << "): " 
               << getSystemMatrix(i,j)->getNnz() << endl;
        }
      }

    Matrix< mtl::compressed2D<double> > newA(nM,nM);
    for (int i= 0; i < nM ; i++)
      for (int j= i; j < nM ; j++)
      {
        newA[i][j]= mtl::compressed2D<double>(n,n);
        for(int k= 0; k < nM ; k++)
          newA[i][j]+= trans(A[k][i])*A[k][j];
        if ( j > i )
          newA[j][i]= trans(newA[i][j]);
      }
    newA[iFESpace][iFESpace]+= additionalRow;

    if (verbose) cout << "\nafter upgrade:\n";
    for (int i= 0; i < nM ; i++)
      for (int j= 0; j < nM ; j++)
      {
        getSystemMatrix(i,j)->getBaseMatrix()= newA[i][j];
        if (verbose)
        {
          getSystemMatrix(i,j)->calculateNnz();
          cout << "nnz per row (" << i << "," << j << "): " 
               << getSystemMatrix(i,j)->getNnz() << endl;
        }
      }



          
    Vector< mtl::dense_vector<double> > b(nM);
    for (int i= 0; i < nM; i++)
    {
      b[i]= mtl::dense_vector<double>(n,0.0);
      for (int j= 0; j < n; j++)
        (b[i])[j]= (*(getRhs()->getDOFVector(i)))[j];
    }

    Vector< mtl::dense_vector<double> > newb(nM);
    for (int i= 0; i < nM ; i++)
    {
      newb[i]= mtl::dense_vector<double>(n,0.0);
      for (int k= 0; k < nM ; k++)
        newb[i]+= trans(A[k][i])*b[k];
    }
    (newb[iFESpace])[iDOF]+= val;

    for (int i= 0; i < nM; i++)
      for (int j= 0; j < n; j++)
        (*(getRhs()->getDOFVector(i)))[j]= (newb[i])[j];

    MSG("... OK!\n");
  };
