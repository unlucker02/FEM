#include "pch.h"

void calcLU(matrixProfile &mp, SLAE &A, SLAE &LU)
{
   auto &ig = mp.ig, &jg = mp.jg;
   auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di, &L = LU.ggl, &U = LU.ggu, &diL = LU.di;
   LU.b = A.b;
   const int n = di.size();
   diL.resize(n);
   const int sizeMat = ig[n];
   L.resize(sizeMat);
   U.resize(sizeMat);

   for (int i = 0; i < n; i++)
   {
      double sumDi = 0;
      int i0 = ig[i];
      int i1 = ig[i + 1];
      for (int k = i0; k < i1; k++)
      {
         double suml = 0, sumu = 0;
         int j = jg[k];
         int j0 = ig[j];
         int j1 = ig[j + 1];
         for (int ik = i0, kj = j0; ik < i1 && kj < j1; )
         {
            if (jg[ik] > jg[kj]) kj++;
            else if (jg[ik] < jg[kj]) ik++;
            else
            {
               suml += L[ik] * U[kj];
               sumu += L[kj] * U[ik];
               ik++;
               kj++;
            }
         }
         L[k] = (ggl[k] - suml);
         U[k] = (ggu[k] - sumu) / diL[j];
         sumDi += L[k] * U[k];
      }
      diL[i] = di[i] - sumDi;
   }
}

void localOptimalSchemeLU(matrixProfile &mp, SLAE &A, SLAE &LU, LOS &v, vector <double> &p, int maxIter, double eps)
{
   auto &r1 = v.r1, &z1 = v.z1, &p1 = v.p1, &mult = v.mult, &rk = v.rk, &Ar = v.Ar, &q = v.q;
   const int n = A.di.size();
   double normb = 0;
   p.resize(n);
   r1.resize(n);
   z1.resize(n);
   p1.resize(n);
   mult.resize(n);
   rk.resize(n);
   Ar.resize(n);
   q.resize(n);
   calcDiscrepancy(mp, A, v, p, normb);
   calcY(mp, LU, r1, r1);
   calcX(mp, LU, r1, z1);
   multOfMatrix(mp, A, z1, p1);
   calcY(mp, LU, p1, p1);
   double scalarr = scalarMult(r1, r1),
      discrepancy = sqrt(scalarr / normb);

   for (int k = 1; k < maxIter && discrepancy > eps; k++)
   {
      double scalarp = scalarMult(p1, p1),
         alpha = scalarMult(p1, r1) / scalarp;
      calcVectorMultCoef(z1, alpha, mult);
      calcSumVectors(p, v.mult, p);
      calcVectorMultCoef(p1, -alpha, mult);
      calcSumVectors(r1, mult, r1);

      calcX(mp, LU, r1, rk);
      multOfMatrix(mp, A, rk, Ar);
      calcY(mp, LU, Ar, q);
      double betta = -scalarMult(p1, q) / scalarp;
      calcVectorMultCoef(z1, betta, mult);
      calcSumVectors(rk, mult, z1);
      calcVectorMultCoef(p1, betta, mult);
      calcSumVectors(q, mult, p1);
      discrepancy = sqrt(scalarMult(r1, r1) / scalarr);
      std::cout << k << " " << discrepancy << std::endl;
   }
   normb = 0;
   calcDiscrepancy(mp, A, v, p, normb);
   discrepancy = sqrt(scalarMult(r1, r1) / normb);
   std::cout << "Final discrepancy: " << discrepancy << std::endl;
}

void multOfMatrix(matrixProfile &mp, SLAE &A, vector <double> &x, vector <double> &F)
{
   auto &ig = mp.ig, &jg = mp.jg;
   auto &di = A.di, &ggl = A.ggl, &ggu = A.ggu;
   const int n = di.size();

   for (int i = 0; i < n; i++)
   {
      F[i] = di[i] * x[i];
      int i0 = ig[i], i1 = ig[i + 1];
      for (i0; i0 < i1; i0++)
      {
         int j = jg[i0];
         F[i] += ggl[i0] * x[j];
         F[j] += ggu[i0] * x[i];
      }
   }
}

void calcY(matrixProfile &mp, SLAE &LU, vector <double> &b, vector <double> &y)
{
   auto &ig = mp.ig, &jg = mp.jg;
   auto &di = LU.di, &L = LU.ggl;
   const int n = di.size();

   for (int i = 0; i < n; i++)
   {
      double sum = 0;
      int i0 = ig[i], i1 = ig[i + 1];
      for (i0; i0 < i1; i0++)
      {
         int j = jg[i0];
         sum += L[i0] * y[j];
      }
      y[i] = (b[i] - sum) / di[i];
   }
}

void calcX(matrixProfile &mp, SLAE &LU, vector <double> &y, vector <double> &x)
{
   auto &ig = mp.ig, &jg = mp.jg;
   auto &U = LU.ggu;
   const int n = ig.size() - 1;
   vector <double> v = y;

   for (int i = n - 1; i >= 0; i--)
   {
      x[i] = v[i];
      int i0 = ig[i], i1 = ig[i + 1];
      for (i0; i0 < i1; i0++)
      {
         int j = jg[i0];
         v[j] -= x[i] * U[i0];
      }
   }
}

void calcDiscrepancy(matrixProfile &mp, SLAE &A, LOS &v, vector <double> &x, double &normb)
{
   auto &ig = mp.ig, &jg = mp.jg;
   auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di, &b = A.b, &r1 = v.r1;
   const int n = di.size();

   for (int i = 0; i < n; i++)
   {
      normb += b[i] * b[i];
      r1[i] = b[i] - di[i] * x[i];
      int i0 = ig[i], i1 = ig[i + 1];
      for (i0; i0 < i1; i0++)
      {
         int j = mp.jg[i0];
         r1[i] -= ggl[i0] * x[j];
         r1[j] -= ggu[i0] * x[i];
      }
   }
}

void calcVectorMultCoef(vector <double> &a, double coef, vector <double> &res)
{
   const int n = a.size();
   for (int i = 0; i < n; i++)
      res[i] = a[i] * coef;
}

void calcSumVectors(vector <double> &a, vector <double> &b, vector <double> &res)
{
   const int n = a.size();
   for (int i = 0; i < n; i++)
      res[i] = a[i] + b[i];
}

double scalarMult(vector <double> &a, vector <double> &b)
{
   int n = a.size();
   double res = 0;
   for (int i = 0; i < n; i++)
      res += a[i] * b[i];
   return res;
}