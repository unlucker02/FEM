#include "pch.h"

void readingArea(estimatedArea &eA)
{
   auto &rW = eA.rW, &zW = eA.zW;
   auto &MW = eA.MW;

   std::ifstream fin("estimatedArea.txt");
   int nrW = 0, nzW = 0, L = 0;
   fin >> nrW;
   rW.resize(nrW);
   for (int i = 0; i < nrW; i++)
      fin >> rW[i];
   fin >> nzW;
   zW.resize(nzW);
   for (int i = 0; i < nzW; i++)
      fin >> zW[i];
   fin >> L;
   MW.resize(L);
   for (int l = 0; l < L; l++)
   {
      MW[l].resize(5);
      fin >> MW[l][0];
      for (int i = 1; i < 5; i++)
      {
         int num = 0;
         fin >> num;
         MW[l][i] = num - 1;
      }
   }
   fin.close();
}

void readingParamersOfAreas(vector <parametersOfAreas> &par)
{
   std::ifstream fin("parametersOfAreas.txt");
   int countAreas = 0;
   fin >> countAreas;
   par.resize(countAreas);

   for (int i = 0; i < countAreas; i++)
   {
      auto &tensor = par[i].tensor;
      auto &kappa = par[i].kappa, &etta = par[i].etta;
      auto &numArea = par[i].numArea;
      int mCount = 0;
      fin >> numArea;
      tensor.resize(2);
      for (int ki = 0; ki < 2; ki++)
      {
         tensor[ki].resize(2);
         for (int kj = 0; kj < 2; kj++)
            fin >> tensor[ki][kj];
      }
      fin >> mCount;
      kappa.resize(mCount);
      etta.resize(mCount);
      for (int m = 0; m < mCount; m++)
         fin >> kappa[m];
      for (int m = 0; m < mCount; m++)
         fin >> etta[m];
   }
}

void readingBoundaryConditions(vector <vector <int>> &bC)
{
   int countConditions = 0;
   std::ifstream fin("boundaryConditions.txt");
   fin >> countConditions;
   bC.resize(countConditions);
   for (int i = 0; i < countConditions; i++)
   {
      bC[i].resize(6);
      fin >> bC[i][0] >> bC[i][1];
      int p0 = 0, p1 = 0, s0 = 0, s1 = 0;
      fin >> p0 >> p1 >> s0 >> s1;
      bC[i][2] = p0 - 1;
      bC[i][3] = p1 - 1;
      bC[i][4] = s0 - 1;
      bC[i][5] = s1 - 1;
   }
}