#include "pch.h"

void generatePortrait(matrixProfile &mp, SLAE &slae, int nR, int nZ)
{
   auto &ig = mp.ig, &jg = mp.jg;
   const int n = nR * nZ;
   vector <std::set<int>> list{ }; // Список, который содержит глобальные номера j баз-ых ф-ий,
   list.resize(n);                 // связанных с баз-ой ф-ей с глобальным номером i
   slae.di.resize(n);

   for (int s = 0; s < nZ - 1; s++)
      for (int p = 0; p < nR - 1; p++)
      {
         vector <int> globalNum{ (s + 1) * nR + p + 1, (s + 1) * nR + p,
                                   s * nR + p + 1, s * nR + p };
         for (int i = 0; i < 4; i++)
         {
            int ind1 = globalNum[i];
            for (int j = i + 1; j < 4; j++)
            {
               int ind2 = globalNum[j];
               list[ind1].insert(ind2);
            }
         }
      }

   // Теперь по списку строим портрет
   ig.resize(n + 1);
   ig[0] = 0;
   ig[1] = 0;
   for (int i = 0; i < n; i++)
      ig[i + 1] = ig[i] + list[i].size();
   int ign = ig[n];
   jg.resize(ign);
   slae.ggl.resize(ign);
   slae.ggu.resize(ign);
   for (int i = 0, k = 0; i < n; i++)
      for (auto j : list[i])
      {
         jg[k] = j;
         k++;
      }
}

void addLocalElement(matrixProfile &mp, SLAE &slae, double elem, int i, int j)
{
   auto &di = slae.di, &ggl = slae.ggl, &ggu = slae.ggu;
   auto &ig = mp.ig, &jg = mp.jg;

   if (i == j)
      di[i] += elem;
   else
   {
      if (i > j)
      {
         int beg = ig[i], end = ig[i + 1] - 1;
         while (jg[beg] != j)
         {
            int ind = (beg + end) / 2;
            if (jg[ind] < j)
               beg = ind + 1;
            else
               end = ind;
         }
         ggl[beg] += elem;
      }
      else
      {
         int beg = ig[j], end = ig[j + 1] - 1;
         while (jg[beg] != i)
         {
            int ind = (beg + end) / 2;
            if (jg[ind] < i)
               beg = ind + 1;
            else
               end = ind;
         }
         ggu[beg] += elem;
      }
   }
}

double calcLambda(vector <double> &kappa, vector <double> &etta)
{
   int mCount = kappa.size();
   double lambda = 0;
   for (int m = 0; m < mCount; m++)
      lambda += kappa[m] / etta[m];
   return lambda;
}

int mu(int i)
{
   return i % 2;
}

int nu(int i)
{
   return i / 2;
}

void boundaryConditions(IArrays &I, grid &g, matrixProfile &mp, SLAE &slae, vector <vector <int>> &bC, functionsBC &func)
{
   auto &ig = mp.ig, &jg = mp.jg, &Ir = I.Ir, &Iz = I.Iz;
   auto &r = g.r, &z = g.z, &di = slae.di, &ggl = slae.ggl, &ggu = slae.ggu, &b = slae.b;
   auto &firstBC = func.firstBC, &secondBC = func.secondBC;
   const int countCond = bC.size(), nR = r.size(), nZ = z.size(), n = ig[di.size()];

   for (int k = 0; k < countCond; k++)
   {
      int typeCond = bC[k][0];
      if (typeCond == 2)
      {
         int numFunc = bC[k][1] - 1, p0 = bC[k][2], p1 = bC[k][3], s0 = bC[k][4], s1 = bC[k][5];
         if (p0 == p1)
         {
            int p = Ir[p0], s = Iz[s0];
            s1 = Iz[s1];
            for (s; s < s1; s++)
            {
               double tetta1 = secondBC[numFunc](r[p], z[s]),
                  tetta2 = secondBC[numFunc](r[p], z[s + 1]), coef = (z[s + 1] - z[s]) / 6.0;
               int globalNum1 = nR * s + p, globalNum2 = nR * (s + 1) + p;
               b[globalNum1] += coef * (2 * tetta1 + tetta2);
               b[globalNum2] += coef * (tetta1 + 2 * tetta2);
            }
         }

         if (s0 == s1)
         {
            int s = Iz[s0], p = Ir[p0];
            p1 = Ir[p1];
            for (p; p < p1; p++)
            {
               double tetta1 = secondBC[numFunc](r[p], z[s]),
                  tetta2 = secondBC[numFunc](r[p + 1], z[s]), hr = r[p + 1] - r[p],
                  coef1 = hr * r[p] / 6.0, coef2 = hr * hr / 12.0;
               int globalNum1 = nR * s + p, globalNum2 = nR * s + p + 1;
               b[globalNum1] += coef1 * (2 * tetta1 + tetta2) + coef2 * (tetta1 + tetta2);
               b[globalNum2] += coef1 * (tetta1 + tetta2) + coef2 * (tetta1 + 3 * tetta2);
            }
         }
      }
   }

   for (int k = 0; k < countCond; k++)
   {
      int typeCond = bC[k][0];
      if (typeCond == 1)
      {
         int numFunc = bC[k][1] - 1, p0 = bC[k][2], p1 = bC[k][3], s0 = bC[k][4], s1 = bC[k][5];
         if (p0 == p1)
         {
            int p = Ir[p0], s = Iz[s0];
            s1 = Iz[s1];
            for (s; s <= s1; s++)
            {
               int globalNum = nR * s + p, i0 = ig[globalNum], i1 = ig[globalNum + 1];
               for (i0; i0 < i1; i0++)
                  ggl[i0] = 0;
               int j0 = ig[globalNum + 1];
               for (j0; j0 < n; j0++)
                  if (jg[j0] == globalNum) ggu[j0] = 0;
               di[globalNum] = 1;
               b[globalNum] = firstBC[numFunc](r[p], z[s]);
            }
         }

         if (s0 == s1)
         {
            int s = Iz[s0], p = Ir[p0];
            p1 = Ir[p1];
            for (p; p <= p1; p++)
            {
               int globalNum = nR * s + p, i0 = ig[globalNum], i1 = ig[globalNum + 1];
               for (i0; i0 < i1; i0++)
                  ggl[i0] = 0;
               int j0 = ig[globalNum + 1];
               for (j0; j0 < n; j0++)
                  if (jg[j0] == globalNum) ggu[j0] = 0;
               di[globalNum] = 1;
               b[globalNum] = firstBC[numFunc](r[p], z[s]);
            }
         }
      }
   }
}


void calcGlobalMatrixAndVector(matrixProfile &mp, SLAE &slae, grid &g, IArrays &I, estimatedArea &eA, vector <parametersOfAreas> &par, matrices &M)
{
   auto &r = g.r, &z = g.z, &b = slae.b;
   auto &G1 = M.G1, &M1 = M.M1, &M2 = M.M2, &H1 = M.H1, &H2 = M.H2;
   auto &Ir = I.Ir, &Iz = I.Iz;
   auto &MW = eA.MW;
   auto &F = slae.F;
   int nR = r.size(), nZ = z.size();
   const int n = nR * nZ;
   b.resize(n);
   int l = 0, countAreas = par.size();

   for (int s = 0; s < nZ - 1; s++)
   {
      for (int p = 0; p < nR - 1; p++)
      {
         int numArea = numberOfEstimatedSubArea(I, MW, p, s, l);
         vector <int> globalNumbers = { nR * s + p, nR * s + p + 1, nR * (s + 1) + p, nR * (s + 1) + p + 1 };

         if (numArea != -1)
         {
            double hz = z[s + 1] - z[s], hr = r[p + 1] - r[p], rp = r[p], coefGr = (2.0 * rp + hr) / (2.0 * hr),
               coefGz = 1.0 / hz, coefMr1 = hr * rp / 6.0, coefMr2 = hr * hr / 12.0,
               coefMz = hz / 6.0, coefHr1 = rp, coefHr2 = hr / 6.0;
            int k = 0;
            for (k = 0; numArea != par[k].numArea; k++);
            auto &kappa = par[k].kappa, &etta = par[k].etta;
            auto &tensor = par[k].tensor;
            double lambda = calcLambda(kappa, etta);
            for (int i = 0; i < 4; i++)
            {
               vector <double> f = { F[k](r[p], z[s]), F[k](r[p + 1], z[s]),
                                     F[k](r[p], z[s + 1]), F[k](r[p + 1], z[s + 1]) };
               double sumbi = 0;
               int mui = mu(i), nui = nu(i);
               for (int j = 0; j < 4; j++)
               {
                  int muj = mu(j), nuj = nu(j);
                  double elemGij =
                     lambda * (tensor[0][0] * coefGr * G1[mui][muj] * coefMz * M1[nui][nuj] +
                        tensor[0][1] * (coefHr1 * H1[mui][muj] + coefHr2 * H2[mui][muj]) * H1[nuj][nui] +
                        tensor[1][0] * (coefHr1 * H1[muj][mui] + coefHr2 * H2[muj][mui]) * H1[nui][nuj] +
                        tensor[1][1] * (coefMr1 * M1[mui][muj] + coefMr2 * M2[mui][muj]) * coefGz * G1[nui][nuj]);
                  addLocalElement(mp, slae, elemGij, globalNumbers[i], globalNumbers[j]);
                  sumbi += f[j] * (coefMr1 * M1[mui][muj] + coefMr2 * M2[mui][muj]) * coefMz * M1[nui][nuj];
               }
               b[globalNumbers[i]] += sumbi;
            }
         }
         else
         {
            auto &di = slae.di;
            vector <vector<int>> localNum = { { p, s }, { p + 1, s }, { p, s + 1 }, { p + 1, s + 1 } };
            int l1 = 0;
            for (int i = 0; i < 4; i++)
               if (IsFictiousNode(I, MW, localNum[i][0], localNum[i][1], l1))
                  di[globalNumbers[i]] = 1;
         }
      }
   }
}