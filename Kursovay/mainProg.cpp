#include "pch.h"

using std::vector;

int main()
{
   estimatedArea eA{ };
   vector <parametersOfAreas> par{ };
   grid g{ };
   IArrays I{ };
   matrices M{ };
   vector <vector <int>> bC{ };
   functionsBC func{ };
   vector <double> p{ };
   SLAE slae{ }, LU{ };
   matrixProfile mp{ };
   LOS v{ };

   readingArea(eA);
   readingParamersOfAreas(par);
   createGrid(eA, g, I);
   generatePortrait(mp, slae, g.r.size(), g.z.size());
   calcGlobalMatrixAndVector(mp, slae, g, I, eA, par, M);
   readingBoundaryConditions(bC);
   boundaryConditions(I, g, mp, slae, bC, func);
   calcLU(mp, slae, LU);
   localOptimalSchemeLU(mp, slae, LU, v, p, 10000, 1e-15);
   double r = 1.1, z = 1.1;
   //double valueFunc = valueFuncAtPoint(g, r, z, p);
   //std::cout << valueFunc << std::endl;
   outputVector(p);
   //printf_s("\n%.15lf\n", valueFunc);
   printf_s("\n");
   //for (double z = 0; z <= 2; z += 1)
      for (double z = 0.5, r = 2; z <= 2.; z += 0.5)
      {
         double valueFunc = valueFuncAtPoint(g, r, z, p);
         printf_s("%.15lf\n", valueFunc);
      }

   return 0;
}

double valueFuncAtPoint(grid &g, double r, double z, vector <double> &p)
{
   auto &gR = g.r, &gZ = g.z;
   const int nR = gR.size(), nZ = gZ.size();
   if (r < gR[0] || r > gR[nR - 1] || z < gZ[0] || z > gZ[nZ - 1])
   {
      std::cout << "The point does not belong to the area." << std::endl;
      return -1;
   }
   
   int begR = 0, endR = nR - 1, begZ = 0, endZ = nZ - 1;

   while (!(gR[begR] <= r && r <= gR[begR + 1]))
   {
      int indR = (begR + endR) / 2;
      if (gR[indR] < r)
         begR = indR;
      else
         endR = indR;
   }

   while (!(gZ[begZ] <= z && z <= gZ[begZ + 1]))
   {
      int indZ = (begZ + endZ) / 2;
      if (gZ[indZ] < z)
         begZ = indZ;
      else
         endZ = indZ;
   }

   vector <int> globalNumbers = { begZ * nR + begR, begZ * nR + begR + 1, (begZ + 1) * nR + begR, (begZ + 1) * nR + begR + 1 };

   double r0 = gR[begR], r1 = gR[begR + 1], z0 = gZ[begZ], z1 = gZ[begZ + 1], hr = r1 - r0, hz = z1 - z0;
   double valueFuncP = p[globalNumbers[0]] * (r1 - r) / hr * (z1 - z) / hz +
                       p[globalNumbers[1]] * (r - r0) / hr * (z1 - z) / hz +
                       p[globalNumbers[2]] * (r1 - r) / hr * (z - z0) / hz +
                       p[globalNumbers[3]] * (r - r0) / hr * (z - z0) / hz;
   return valueFuncP;
}

void outputVector(vector <double> &p)
{
   const int n = p.size();

   printf_s("p:\n");
   for (int i = 0; i < n; i++)
      printf_s("%.15lf\n", p[i]);
}



