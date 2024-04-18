#include "pch.h"

void createGrid(estimatedArea &eA, grid &g, IArrays &I)
{
   auto &r = g.r, &z = g.z, &rW = eA.rW, &zW = eA.zW;
   auto &Ir = I.Ir, &Iz = I.Iz;
   const int nR = rW.size() - 1, nZ = zW.size() - 1;
   Ir.resize(nR + 1);
   Iz.resize(nZ + 1);
   std::ifstream fin("grid.txt");
   int nRk = 0;
   r.resize(1, rW[0]);
   for (int i = 0, j = 1; i < nR; i++, j++)
   {
      int countInterval = 0;
      double coef = 0, step = 0;
      fin >> countInterval >> coef;
      nRk += countInterval;
      r.resize(nRk + 1);
      if (coef != 1)
      {
         double sumProg = (pow(coef, countInterval) - 1.0) / (coef - 1.0);
         step = (rW[i + 1] - rW[i]) / sumProg;
      }
      else
         step = (rW[i + 1] - rW[i]) / countInterval;
      int jk = 1;
      for (j; j < nRk; j++, jk++)
         if (coef != 1)
            r[j] = rW[i] + step * (pow(coef, jk) - 1.0) / (coef - 1.0);
         else
            r[j] = rW[i] + step * jk;
      r[j] = rW[i + 1];
      Ir[i + 1] = j;
   }

   int nZk = 0;
   z.resize(1, zW[0]);
   for (int i = 0, j = 1; i < nZ; i++, j++)
   {
      int countInterval = 0;
      double coef = 0, step = 0;
      fin >> countInterval >> coef;
      nZk += countInterval;
      z.resize(nZk + 1);
      if (coef != 1)
      {
         double sumProg = (pow(coef, countInterval) - 1.0) / (coef - 1.0);
         step = (zW[i + 1] - zW[i]) / sumProg;
      }
      else
         step = (zW[i + 1] - zW[i]) / countInterval;
      int jk = 1;
      for (j; j < nZk; j++, jk++)
         if (coef != 1)
            z[j] = zW[i] + step * (pow(coef, jk) - 1.0) / (coef - 1.0);
         else
            z[j] = zW[i] + step * jk;
      z[j] = zW[i + 1];
      Iz[i + 1] = j;
   }
   fin.close();
}

int numberOfEstimatedSubArea(IArrays &I, vector <vector <int>> &MW, int p, int s, int &l)
{
   auto &Ir = I.Ir, &Iz = I.Iz;
   const int L = MW.size(), L1 = l;

   for (l; l < L; l++)
   {
      int mr0 = Ir[MW[l][1]], mr1 = Ir[MW[l][2]],
         mz0 = Iz[MW[l][3]], mz1 = Iz[MW[l][4]], m = MW[l][0];
      if (mr0 <= p && p <= mr1 && mr0 <= (p + 1) && (p + 1) <= mr1 &&
         mz0 <= s && s <= mz1 && mz0 <= (s + 1) && (s + 1) <= mz1)
         return m;
   }

   for (l = 0; l < L1; l++)
   {
      int mr0 = Ir[MW[l][1]], mr1 = Ir[MW[l][2]],
         mz0 = Iz[MW[l][3]], mz1 = Iz[MW[l][4]], m = MW[l][0];
      if (mr0 <= p && p <= mr1 && mr0 <= (p + 1) && (p + 1) <= mr1 &&
         mz0 <= s && s <= mz1 && mz0 <= (s + 1) && (s + 1) <= mz1)
         return m;
   }
   return -1;
}

bool IsFictiousNode(IArrays &I, vector <vector <int>> &MW, int p, int s, int &l)
{
   auto &Ir = I.Ir, &Iz = I.Iz;
   const int L = MW.size(), L1 = l;

   for (l; l < L; l++)
   {
      int mr0 = Ir[MW[l][1]], mr1 = Ir[MW[l][2]],
         mz0 = Iz[MW[l][3]], mz1 = Iz[MW[l][4]];
      if (mr0 <= p && p <= mr1 && mz0 <= s && s <= mz1)
         return false;
   }
   for (l = 0; l < L1; l++)
   {
      int mr0 = Ir[MW[l][1]], mr1 = Ir[MW[l][2]],
         mz0 = Iz[MW[l][3]], mz1 = Iz[MW[l][4]];
      if (mr0 <= p && p <= mr1 && mz0 <= s && s <= mz1)
         return false;
   }
   return true;
}
