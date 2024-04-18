#include "pch.h"

using std::vector;

struct estimatedArea
{
   vector <double> rW{ }, zW{ };
   vector <vector <int>> MW{ };
};

struct parametersOfAreas
{
   int numArea = 0;
   vector <vector <double>> tensor{ };
   vector <double> kappa{ }, etta{ };
};

struct grid
{
   vector <double> r{ }, z{ };
};

struct matrices
{
   vector <vector <double>> G1 = { { 1, -1 }, { -1, 1 } },
      M1 = { { 2, 1 }, { 1, 2 } },
      M2 = { { 1, 1 }, { 1, 3 } },
      H1 = { { -0.5, 0.5 }, { -0.5, 0.5 } },
      H2 = { { -1, 1 }, { -2, 2 } };
};

struct matrixProfile
{
   vector <int> ig{ }, jg{ };
};

struct SLAE
{
   vector <double> di{ }, ggl{ }, ggu{ }, b{ };
   vector <std::function<double(double, double)>> F
   {
      [](double r, double z) { return sin(z) + cos(z); },
      [](double r, double z) { return 0; },
      [](double r, double z) { return 0; }
   };
};

struct functionsBC
{
   vector <std::function<double(double, double)>> firstBC
   {
      [](double r, double z) { return 1; },
      [](double r, double z) { return cos(z) + sin(z); },
      [](double r, double z) { return cos(32) + sin(32); },
      [](double r, double z) { return 1 + z; }
   };

   vector <std::function<double(double, double)>> secondBC
   {
      [](double r, double z) { return 0; },
      [](double r, double z) { return 7 + z; },
      [](double r, double z) { return r + 5; },
      [](double r, double z) { return 1 + z; }
   };
};

struct LOS
{
   vector <double> r1{ }, rk{ }, z1{ }, p1{ }, Ar{ }, q{ }, mult{ };
};

struct IArrays
{
   vector <int> Ir{ }, Iz{ };
};

void readingArea(estimatedArea &eA);
void readingParamersOfAreas(vector <parametersOfAreas> &par);
void readingBoundaryConditions(vector <vector <int>> &bC);
void createGrid(estimatedArea &eA, grid &g, IArrays &I);
int numberOfEstimatedSubArea(IArrays &I, vector <vector <int>> &MW, int p, int s, int &l);
bool IsFictiousNode(IArrays &I, vector <vector <int>> &MW, int p, int s, int &l);
void generatePortrait(matrixProfile &mp, SLAE &slae, int nR, int nZ);
void addLocalElement(matrixProfile &mp, SLAE &slae, double elem, int i, int j);
double calcLambda(vector <double> &kappa, vector <double> &etta);
int mu(int i);
int nu(int i);
void boundaryConditions(IArrays &I, grid &g, matrixProfile &mp, SLAE &slae, vector <vector <int>> &bC, functionsBC &func);
void calcGlobalMatrixAndVector(matrixProfile &mp, SLAE &slae, grid &g, IArrays &I, estimatedArea &eA, vector <parametersOfAreas> &par, matrices &M);
void calcLU(matrixProfile &mp, SLAE &A, SLAE &LU);
void localOptimalSchemeLU(matrixProfile &mp, SLAE &A, SLAE &LU, LOS &v, vector <double> &p, int maxIter, double eps);
void multOfMatrix(matrixProfile &mp, SLAE &A, vector <double> &x, vector <double> &F);
void calcY(matrixProfile &mp, SLAE &LU, vector <double> &y, vector <double> &b);
void calcX(matrixProfile &mp, SLAE &LU, vector <double> &y, vector <double> &x);
void calcDiscrepancy(matrixProfile &mp, SLAE &A, LOS &v, vector <double> &x, double &normb);
void calcVectorMultCoef(vector <double> &a, double coef, vector <double> &res);
void calcSumVectors(vector <double> &a, vector <double> &b, vector <double> &res);
double scalarMult(vector <double> &a, vector <double> &b);
double valueFuncAtPoint(grid &g, double r, double z, vector <double> &p);
void outputVector(vector <double> &p);