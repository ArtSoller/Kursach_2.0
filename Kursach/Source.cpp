#include "Header.h"

namespace solution
{
   void Difur::Grid_Generator()
   {
      ifstream in("boundaries.txt");
      in >> x1 >> x2 >> kx >> y1 >> y2 >> ky >> Ne; //вводим границы, количество элементов и коэффициенты релаксации
      in.close();
      int i = 2;
      while (i <= sqrt(Ne))
      {
         if (Ne % i == 0)
         {
            f = true;
            break;
         }
         i++;
      }
      if (f)
      {
         Ny = i;
         Nx = Ne / i;
      }
      else
      {
         Ny = 1;
         Nx = Ne;
      }
      hx.resize(2 * Nx + 1);
      hy.resize(2 * Ny + 1);
      double sum = 0;
      for (int j = 0; j <= Nx - 1; j++)
         sum += pow(kx, j);
      double w = (x2 - x1) / sum;
      hx[0] = x1;
      double nstu;
      for (int i = 1; i <= 2 * Nx; i += 2)
      {
         nstu = hx[i - 1] + w * pow(kx, i - 1);// заполнение разбиения по x
         hx[i] = (nstu + hx[i - 1]) / 2;
         hx[i + 1] = nstu;
         //cout << hx[i] << ' ' << hx[i + 1] << ' ';
      }
      sum = 0;
      for (int j = 0; j <= Ny - 1; j++)
         sum += pow(ky, j);
      w = (y2 - y1) / sum;
      cout << '\n';
      hy[0] = y1;
      for (int i = 1; i <= 2 * Ny; i += 2)
      {
         nstu = hy[i - 1] + w * pow(ky, i - 1);// заполнение разбиения по x
         hy[i] = (nstu + hy[i - 1]) / 2;
         hy[i + 1] = nstu;
         //cout << hy[i] << ' ' << hy[i + 1] << ' ';
      }
      int Nn = (2 * Nx + 1) * (2 * Ny + 1); // кол-во координат узлов
      XY.resize(Nn);
      ofstream out("XY.txt");
      out << Nn << '\n';
      for (int q = 0; q < Nn; q++)
      {
         XY[q].first = hx[q % (2 * Nx + 1)];
         XY[q].second = hy[q / (2 * Nx + 1)];
         out << XY[q].first << ' ' << XY[q].second << '\n';
      }
      out.close();
      elems.resize(Ne);
      int o = 1, u = 0;
      for (int k = 0; k < Ny; k++)
      {
         for (int j = 0; j < Nx; j++)
         {
            if (o % (2 * Nx + 1) == 0)
            {
               o += (2 * Nx + 2);
            }
            elems[u].resize(10);
            elems[u][0] = o;
            elems[u][1] = o + 1;
            elems[u][2] = o + 2;
            elems[u][3] = o + (2 * Nx + 1);
            elems[u][4] = 1 + o + (2 * Nx + 1);
            elems[u][5] = 2 + o + (2 * Nx + 1);
            elems[u][6] = o + (2 * Nx + 1) * 2;
            elems[u][7] = 1 + o + (2 * Nx + 1) * 2;
            elems[u][8] = 2 + o + (2 * Nx + 1) * 2;
            elems[u][9] = 1;
            o += 2;
            u++;
         }
      }
      out.open("elems.txt");
      out << elems.size() << '\n';
      for (int k = 0; k < Ne; k++)
         out << elems[k][0] << ' ' << elems[k][1] << ' ' << elems[k][2] << ' ' << elems[k][3] << ' ' << elems[k][4] << ' ' << elems[k][5] << ' ' << elems[k][6] << ' ' << elems[k][7] << ' ' << elems[k][8] << ' ' << elems[k][9] << ' ' << endl;
      out.close();
   }

   void Difur::PD_Grid_Generator()
   {
      ifstream in("boundaries.txt");
      in >> x1 >> x2 >> y1 >> y2 >> hor >> ver >> kx >> ky >> Perf1 >> Perf2 >> Nperf; //вводим границы, количество элементов и коэффициенты релаксации
      in.close();

      hx.resize(hor);
      hy.resize(ver);
      double sum = 0;
      float h0_perf = (Perf2 - Perf1) / (2 * Nperf);
      Nx = hor * 2 + 1;
      Ny = ver * 2 + 1;
      for (int i = 0; i < (Nx - 1) / 2; i++)
      {
         hx[i] = pow(kx, i);
         sum += hx[i];
      }

      double Hx = (x2 - x1) / (2 * sum); //шаг по r
      sum = 0;
      for (int i = 0; i < (Ny - 1) / 2; i++)
      {
         hy[i] = pow(ky, i);
         sum += hy[i];
      }
      double Hy = (y2 - y1 - (Perf2 - Perf1)) / (2 * sum);
      ver += Nperf;
      NElems = hor * ver;
      elems.resize(NElems);
      XY.resize(Nx * (Ny + 2 * Nperf));
      for (int k = 0; k < hor; k++)
         for (int j = 0; j < ver; j++) // Номер элемента
         {
            elems[k + j * hor].resize(10);
            elems[k + j * hor][0] = 1 + k * 2 + Nx * 2 * j;
            elems[k + j * hor][1] = 2 + k * 2 + Nx * 2 * j;
            elems[k + j * hor][2] = 3 + k * 2 + Nx * 2 * j;
            elems[k + j * hor][3] = 1 + k * 2 + Nx * (2 * j + 1);
            elems[k + j * hor][4] = 2 + k * 2 + Nx * (2 * j + 1);
            elems[k + j * hor][5] = 3 + k * 2 + Nx * (2 * j + 1);
            elems[k + j * hor][6] = 1 + k * 2 + Nx * (2 * j + 2);
            elems[k + j * hor][7] = 2 + k * 2 + Nx * (2 * j + 2);
            elems[k + j * hor][8] = 3 + k * 2 + Nx * (2 * j + 2);
            elems[k + j * hor][9] = 1;
         }
      int a, b; // a - начальный, b - конечный номер точки в зоне перфорации
      a = Ny / 2;
      b = a + 2 * Nperf;
      Ny += 2 * Nperf;
      for (int k = 0; k < Nx * Ny; k++) // Расставляем точки в матрице
      {
         int line = k / Nx; // Номер строки
         int column = k - line * Nx; // Номер столбца
         float rAdd = 0, zAdd = 0;
         if (column == Nx - 1)
            XY[k].first = x2;
         else
         {
            for (int n = 0; n < column && n < Nx / 2; n++)
               rAdd += hx[n] * Hx;
            if (column > Nx / 2)
            {
               for (int n = 0; n < column - Nx / 2; n++)
                  rAdd += Hx * hx[hx.size() - 1 - n]; // берем с конца
            }
            XY[k].first = x1 + rAdd;
         }
         if (line == Ny - 1)
            XY[k].second = y2;
         else
         {
            if (line == a)
               XY[k].second = Perf1;
            else if (line == b)
               XY[k].second = Perf2;
            else if (a < line && line < b)
               XY[k].second = Perf1 + h0_perf * (line - a);
            else if (line < a)
            {
               for (int n = 0; n < line && n < (Ny - Nperf * 2) / 2; n++)
                  zAdd += hy[n] * Hy;
               XY[k].second = y1 + zAdd;
            }

            else if (line > b)
            {
               for (int n = 0; n < line - b; n++)
                  zAdd += Hy * hy[hy.size() - 1 - n]; // берем с конца
               XY[k].second = Perf2 + zAdd;
            }
         }
      }
      ofstream out("XY.txt");
      out << XY.size() << '\n';
      for (int q = 0; q < XY.size(); q++)
         out << XY[q].first << ' ' << XY[q].second << '\n';
      out.close();
      out.open("elems.txt");
      out << elems.size() << '\n';
      for (int k = 0; k < elems.size(); k++)
         out << elems[k][0] << ' ' << elems[k][1] << ' ' << elems[k][2] << ' ' << elems[k][3] << ' ' << elems[k][4] << ' ' << elems[k][5] << ' ' << elems[k][6] << ' ' << elems[k][7] << ' ' << elems[k][8] << ' ' << elems[k][9] << ' ' << endl;
      out.close();
   }

   void Difur::Read()
   {
      ifstream in("XY.txt");// координаты (x, y)
      in >> Nn;
      P.resize(Nn);
      XY.resize(Nn);
      for (int i = 0; i < Nn; i++)
         in >> XY[i].first >> XY[i].second;
      in.close();
      in.open("elems.txt");
      in >> Ne;
      elems.resize(Ne);
      for (int i = 0; i < Ne; i++)
      {
         elems[i].resize(10);
         in >> elems[i][0] >> elems[i][1] >> elems[i][2] >> elems[i][3] >> elems[i][4] >> elems[i][5] >> elems[i][6] >> elems[i][7] >> elems[i][8] >> elems[i][9];
      }
      in.close();
      int N;
      in.open("mat.txt");
      in >> N;
      mat.resize(N);
      for (int i = 0; i < N; i++)
         in >> mat[i].first >> mat[i].second;
      in.close();
      in.open("b2.txt");// краевое условие II рода
      in >> N;
      b2.resize(N);
      for (int i = 0; i < N; i++)
      {
         b2[i].resize(3);
         in >> b2[i][0] >> b2[i][1] >> b2[i][2];
      }
      in.close();
      in.open("b1.txt");// краевое условие I рода
      in >> N;
      b1.resize(N);
      for (int i = 0; i < N; i++)
         in >> b1[i].first >> b1[i].second;
      in.close();
   }

   void Difur::Portrait()
   {
      vector<vector<int>> adjacency;
      adjacency.resize(Nn);
      for (int i = 0; i < elems.size(); i++) //Заполняем матрицу смежности
      {
         int n = adjacency[elems[i][1] - 1].size();
         adjacency[elems[i][1] - 1].resize(n + 1);
         adjacency[elems[i][1] - 1][n] = elems[i][0];

         n = adjacency[elems[i][2] - 1].size();
         adjacency[elems[i][2] - 1].resize(n + 2);
         adjacency[elems[i][2] - 1][n] = elems[i][0];
         adjacency[elems[i][2] - 1][n + 1] = elems[i][1];

         n = adjacency[elems[i][3] - 1].size();
         adjacency[elems[i][3] - 1].resize(n + 3);
         adjacency[elems[i][3] - 1][n] = elems[i][0];
         adjacency[elems[i][3] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][3] - 1][n + 2] = elems[i][2];

         n = adjacency[elems[i][4] - 1].size();
         adjacency[elems[i][4] - 1].resize(n + 4);
         adjacency[elems[i][4] - 1][n] = elems[i][0];
         adjacency[elems[i][4] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][4] - 1][n + 2] = elems[i][2];
         adjacency[elems[i][4] - 1][n + 3] = elems[i][3];

         n = adjacency[elems[i][5] - 1].size();
         adjacency[elems[i][5] - 1].resize(n + 5);
         adjacency[elems[i][5] - 1][n] = elems[i][0];
         adjacency[elems[i][5] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][5] - 1][n + 2] = elems[i][2];
         adjacency[elems[i][5] - 1][n + 3] = elems[i][3];
         adjacency[elems[i][5] - 1][n + 4] = elems[i][4];

         n = adjacency[elems[i][6] - 1].size();
         adjacency[elems[i][6] - 1].resize(n + 6);
         adjacency[elems[i][6] - 1][n] = elems[i][0];
         adjacency[elems[i][6] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][6] - 1][n + 2] = elems[i][2];
         adjacency[elems[i][6] - 1][n + 3] = elems[i][3];
         adjacency[elems[i][6] - 1][n + 4] = elems[i][4];
         adjacency[elems[i][6] - 1][n + 5] = elems[i][5];

         n = adjacency[elems[i][7] - 1].size();
         adjacency[elems[i][7] - 1].resize(n + 7);
         adjacency[elems[i][7] - 1][n] = elems[i][0];
         adjacency[elems[i][7] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][7] - 1][n + 2] = elems[i][2];
         adjacency[elems[i][7] - 1][n + 3] = elems[i][3];
         adjacency[elems[i][7] - 1][n + 4] = elems[i][4];
         adjacency[elems[i][7] - 1][n + 5] = elems[i][5];
         adjacency[elems[i][7] - 1][n + 6] = elems[i][6];

         n = adjacency[elems[i][8] - 1].size();
         adjacency[elems[i][8] - 1].resize(n + 8);
         adjacency[elems[i][8] - 1][n] = elems[i][0];
         adjacency[elems[i][8] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][8] - 1][n + 2] = elems[i][2];
         adjacency[elems[i][8] - 1][n + 3] = elems[i][3];
         adjacency[elems[i][8] - 1][n + 4] = elems[i][4];
         adjacency[elems[i][8] - 1][n + 5] = elems[i][5];
         adjacency[elems[i][8] - 1][n + 6] = elems[i][6];
         adjacency[elems[i][8] - 1][n + 7] = elems[i][7];

      }

      for (int i = 0; i < Nn; i++) // ОТсортируем по возрастанию и удалим повторы, для составления в последствии векторов ig и jg
      {
         sort(adjacency[i].begin(), adjacency[i].end());
         auto last = unique(adjacency[i].begin(), adjacency[i].end());
         adjacency[i].erase(last, adjacency[i].end());
      }

      di.resize(Nn);
      B.resize(Nn);
      ig.resize(Nn + 1);
      ig[0] = ig[1] = 1;
      for (int i = 2; i < Nn + 1; i++) //Заполняем массив ig
         ig[i] = ig[i - 1] + adjacency[i - 1].size();
      for (int i = 1; i < Nn; i++) //Заполняем массив jg
         jg.insert(jg.end(), adjacency[i].begin(), adjacency[i].end());
      ggu.resize(ig.back() - 1);
      ggl.resize(ig.back() - 1);

      //for (vector<int> row : adjacency)
      //{
      //   for (int elem : row)
      //      cout << elem << " ";
      //   cout << endl;
      //}
   }

   void Difur::output()
   {
      for (int i = 0; i < P.size(); i++)
         cout << P[i] << ' ';
   }

   void Difur::GaMbo()
   {
      vector<vector<double>> G, M;// локальная матрица жесткости
      vector<double> b;// локальный вектор правой части
      G.resize(9);
      M.resize(9);
      for (int i = 0; i < 9; i++)
         G[i].resize(9);
      for (int i = 0; i < 9; i++)
         M[i].resize(9);
      b.resize(9);
      double hx, hy, l, g; // шаги по Х и Y, лямбда, гамма
      for (int i = 0; i < Ne; i++)
      {
         hx = (XY[elems[i][8] - 1].first - XY[elems[i][0] - 1].first);
         hy = (XY[elems[i][8] - 1].second - XY[elems[i][0] - 1].second);
         l = mat[elems[i][9] - 1].first;
         g = mat[elems[i][9] - 1].second;

         G[0][0] = G[2][2] = G[6][6] = G[8][8] = l * (hy / hx * 28 / 90 + hx / hy * 28 / 90);
         G[0][1] = G[1][0] = G[1][2] = G[2][1] = G[6][7] = G[7][6] = G[7][8] = G[8][7] = l * (hy / hx * (-32) / 90 + hx / hy * 14 / 90);
         G[0][2] = G[2][0] = G[6][8] = G[8][6] = l * (hy / hx * 4 / 90 + hx / hy * (-7) / 90);
         G[0][3] = G[2][5] = G[3][0] = G[3][6] = G[5][2] = G[5][8] = G[6][3] = G[8][5] = l * (hy / hx * 14 / 90 + hx / hy * (-32) / 90);
         G[0][4] = G[1][3] = G[1][5] = G[2][4] = G[3][1] = G[3][7] = G[4][0] = G[4][2] = G[4][6] = G[4][8] = G[5][1] = G[5][7] = G[6][4] = G[7][3] = G[7][5] = G[8][4] = l * (hy / hx * (-16) / 90 + hx / hy * (-16) / 90);
         G[0][5] = G[2][3] = G[3][2] = G[3][8] = G[5][0] = G[5][6] = G[6][5] = G[8][3] = l * (hy / hx * 2 / 90 + hx / hy * 8 / 90);
         G[0][6] = G[2][8] = G[6][0] = G[8][2] = l * (hy / hx * (-7) / 90 + hx / hy * 4 / 90);
         G[0][7] = G[1][6] = G[1][8] = G[2][7] = G[6][1] = G[7][0] = G[7][2] = G[8][1] = l * (hy / hx * 8 / 90 + hx / hy * 2 / 90);
         G[0][8] = G[2][6] = G[6][2] = G[8][0] = l * (hy / hx * (-1) / 90 + hx / hy * (-1) / 90);
         G[1][1] = G[7][7] = l * (hy / hx * 64 / 90 + hx / hy * 112 / 90);
         G[1][4] = G[4][1] = G[4][7] = G[7][4] = l * (hy / hx * 32 / 90 + hx / hy * (-128) / 90);
         G[1][7] = G[7][1] = l * (hy / hx * (-16) / 90 + hx / hy * 16 / 90);
         G[3][3] = G[5][5] = l * (hy / hx * 112 / 90 + hx / hy * 64 / 90);
         G[3][4] = G[4][3] = l * (hy / hx * (-128) / 90 + hx / hy * 32 / 90);
         G[3][5] = G[5][3] = l * (hy / hx * 16 / 90 + hx / hy * (-16) / 90);
         G[4][4] = l * (hy / hx * 256 / 90 + hx / hy * 256 / 90);
         G[4][5] = G[5][4] = l * (hy / hx * (-128) / 90 + hx / hy * 32 / 90);

         M[0][0] = g * hx * hy / 900 * 16;  M[0][1] = g * hx * hy / 900 * 8;     M[0][2] = g * hx * hy / 900 * (-4);   M[0][3] = g * hx * hy / 900 * 8;
         M[0][4] = hx * hy / 900 * 4;     M[0][5] = g * hx * hy / 900 * (-2);  M[0][6] = g * hx * hy / 900 * (-4);   M[0][7] = g * hx * hy / 900 * (-2);
         M[0][8] = g * hx * hy / 900 * 1;

         M[1][0] = g * hx * hy / 900 * 8;   M[1][1] = g * hx * hy / 900 * 64;  M[1][2] = g * hx * hy / 900 * 8;    M[1][3] = g * hx * hy / 900 * 4;
         M[1][4] = g * hx * hy / 900 * 32;  M[1][5] = g * hx * hy / 900 * 4;   M[1][6] = g * hx * hy / 900 * (-2); M[1][7] = g * hx * hy / 900 * (-16);
         M[1][8] = g * hx * hy / 900 * (-2);

         M[2][0] = g * hx * hy / 900 * (-4);  M[2][1] = g * hx * hy / 900 * 8;  M[2][2] = g * hx * hy / 900 * 16;  M[2][3] = g * hx * hy / 900 * (-2);
         M[2][4] = g * hx * hy / 900 * 4;     M[2][5] = g * hx * hy / 900 * 8;  M[2][6] = g * hx * hy / 900 * 1;   M[2][7] = g * hx * hy / 900 * (-2);
         M[2][8] = g * hx * hy / 900 * (-4);

         M[3][0] = g * hx * hy / 900 * 8;  M[3][1] = g * hx * hy / 900 * 4;     M[3][2] = g * hx * hy / 900 * (-2);  M[3][3] = g * hx * hy / 900 * 64;
         M[3][4] = g * hx * hy / 900 * 32; M[3][5] = g * hx * hy / 900 * (-16); M[3][6] = g * hx * hy / 900 * 8;     M[3][7] = g * hx * hy / 900 * 4;
         M[3][8] = g * hx * hy / 900 * (-2);

         M[4][0] = g * hx * hy / 900 * 4;    M[4][1] = g * hx * hy / 900 * 32;  M[4][2] = g * hx * hy / 900 * 4;  M[4][3] = g * hx * hy / 900 * 32;
         M[4][4] = g * hx * hy / 900 * 256;  M[4][5] = g * hx * hy / 900 * 32;  M[4][6] = g * hx * hy / 900 * 4;  M[4][7] = g * hx * hy / 900 * 32;
         M[4][8] = g * hx * hy / 900 * 4;

         M[5][0] = g * hx * hy / 900 * (-2);  M[5][1] = g * hx * hy / 900 * 4;   M[5][2] = g * hx * hy / 900 * 8;     M[5][3] = g * hx * hy / 900 * (-16);
         M[5][4] = g * hx * hy / 900 * 32;    M[5][5] = g * hx * hy / 900 * 64;  M[5][6] = g * hx * hy / 900 * (-2);  M[5][7] = g * hx * hy / 900 * 4;
         M[5][8] = g * hx * hy / 900 * 8;

         M[6][0] = g * hx * hy / 900 * (-4);  M[6][1] = g * hx * hy / 900 * (-2);  M[6][2] = g * hx * hy / 900 * 1;   M[6][3] = g * hx * hy / 900 * 8;
         M[6][4] = g * hx * hy / 900 * 4;     M[6][5] = g * hx * hy / 900 * (-2);  M[6][6] = g * hx * hy / 900 * 16;  M[6][7] = g * hx * hy / 900 * 8;
         M[6][8] = g * hx * hy / 900 * (-4);

         M[7][0] = g * hx * hy / 900 * (-2);  M[7][1] = g * hx * hy / 900 * (-16);  M[7][2] = g * hx * hy / 900 * (-2);   M[7][3] = g * hx * hy / 900 * 4;
         M[7][4] = g * hx * hy / 900 * 32;    M[7][5] = g * hx * hy / 900 * 4;      M[7][6] = g * hx * hy / 900 * 8;      M[7][7] = g * hx * hy / 900 * 64;
         M[7][8] = g * hx * hy / 900 * 8;

         M[8][0] = g * hx * hy / 900 * 1;    M[8][1] = g * hx * hy / 900 * (-2);  M[8][2] = g * hx * hy / 900 * (-4);   M[8][3] = g * hx * hy / 900 * (-2);
         M[8][4] = g * hx * hy / 900 * 4;    M[8][5] = g * hx * hy / 900 * 8;     M[8][6] = g * hx * hy / 900 * (-4);   M[8][7] = g * hx * hy / 900 * 8;
         M[8][8] = g * hx * hy / 900 * 16;

         for (int k = 0; k < 9; k++) //Заполняем локальный вектор
            for (int j = 0; j < 9; j++)
               b[k] += M[k][j] * Func(XY[elems[i][j] - 1].first, XY[elems[i][j] - 1].second);

         for (int j = 0; j < 9; j++)//Заполняем глобальный вектор и главную диагональ глобальной матрицы
         {
            B[elems[i][j] - 1] += b[j];
            di[elems[i][j] - 1] += G[j][j] +M [j][j];
         }

         int column, line;
         int n = 0;
         for (int j = 0; j < 9; j++) //вносим локальную матрицу в глобальную
            for (int k = 0; k < 9; k++)
            {
               n = 0;//столбец
               if (elems[i][j] > elems[i][k])
               {
                  line = elems[i][j] - 1;
                  column = elems[i][k];
                  while (jg[ig[line] - 1 + n] != column)
                     n++;
                  ggl[ig[line] - 1 + n] += G[j][k] + M[j][k];
               }
               else if (elems[i][j] < elems[i][k])
               {
                  line = elems[i][k] - 1;
                  column = elems[i][j];
                  while (jg[ig[line] - 1 + n] != column)
                     n++;
                  ggu[ig[line] - 1 + n] += G[k][j] + M[k][j];
               }
            }
      }
   }

   /*double Difur::Solution_in_point(double x, double y)
   {
      for (int i = 0; i < elems.size(); i++)
      {
         if (XY[elems[i][0] - 1].first < x && XY[elems[i][3] - 1].first > x && XY[elems[i][0] - 1].second < y && XY[elems[i][3] - 1].second > y)
 
      }*/
      /*  }*/


   double Difur::Tetta(int n, double x, double y) //
   {
      switch (n)
      {
      case -1:
         return -1;
         break;
         break;
      case 1:
         return 1;
         break;
      case 6:
         return 6;
         break;
      case -2:
         return -2;
         break;
      }
   }

   void Difur::Boundary_conditions() // учет краевых условий
   {
      for (int i = 0; i < b2.size(); i++) // краевое условие II рода
      {
         int el = b2[i][0];
         int f, s, t; // первый, второй и третий узел грани
         double h; // шаг
         switch ((int)b2[i][1])
         {
         case 1:
            f = elems[el - 1][0];
            s = elems[el - 1][1];
            t = elems[el - 1][2];
            h = (XY[t - 1].first - XY[f - 1].first);
            break;
         case 2:
            f = elems[el - 1][2];
            s = elems[el - 1][5];
            t = elems[el - 1][8];
            h = (XY[t - 1].second - XY[f - 1].second);
            break;
         case 3: //верхняя грань
            f = elems[el - 1][6];
            s = elems[el - 1][7];
            t = elems[el - 1][8];
            h = (XY[t - 1].first - XY[f - 1].first);
            break;
         case 4: //левая грань
            f = elems[el - 1][0];
            s = elems[el - 1][3];
            t = elems[el - 1][6];
            h = (XY[t - 1].second - XY[f - 1].second);
            break;
         }
         B[f - 1] += h / 30 * (4 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) - Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
         B[s - 1] += h / 30 * (2 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 16 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 2 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
         B[t - 1] += h / 30 * ((-1) * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 4 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
      }

      for (int i = 0; i < b1.size(); i++) // краевое условие I рода
      {
         double n = b1[i].first;// Номер краевого узла
         di[n - 1] = 1;
         for (int j = 0; j < ig[n] - ig[n - 1]; j++)// Обнуляем строчку нижнего треугольника
            ggl[ig[n - 1] + j - 1] = 0;

         for (int i = n; i < Nn; i++)// Зануляем столбцы в верхнем треугольнике, пробегая по ним
         {
            int k = -1;
            int j = 0;
            do
            {
               if (jg[ig[i] - 1 + j] == n)
                  k = j;
               j++;
            } while (j < ig[i + 1] - ig[i]);
            if (k != -1)
               ggu[ig[i] - 1 + k] = 0;
         }
         B[n - 1] = b1[i].second;
      }
   }

   void Difur::LU()
   {
      d.resize(di.size());
      gl.resize(ggl.size());
      gu.resize(ggu.size());
      int i0, i1, j, kk, ll, ll_1, lll, kkk;
      double sd, sl, su;
      for (int i = 0; i < Nn; i++)
      {
         i0 = ig[i] - 1;
         i1 = ig[i + 1] - 1;
         sd = 0;
         for (int jj = i0; jj < i1; jj++)
         {
            j = jg[jj] - 1;
            kk = i0;
            ll = ig[j] - 1;
            ll_1 = ig[j + 1] - 1;
            sl = 0;
            su = 0;
            while (ll < ll_1 && kk < jj)
            {
               lll = jg[ll] - 1;
               kkk = jg[kk] - 1;

               if (lll == kkk)
               {
                  su += gl[ll] * gu[kk];
                  sl += gu[ll] * gl[kk];
                  kk++;
                  ll++;
               }
               else
               {
                  if (kkk < lll)
                     kk++;
                  else
                     ll++;
               }
            }
            gl[jj] = (ggl[jj] - sl);
            gu[jj] = (ggu[jj] - su) / d[j];
            sd += gu[jj] * gl[jj];
         }
         d[i] = di[i] - sd;
      }

      for (int i = 0; i < Nn; i++)// начальное приближение равно 0
         P[i] = 0;
   }

   void Difur::mult_A(vector<double>& vect, vector<double>& res)
   {
      int i0, i1, j;
      for (int i = 0; i < Nn; i++)
      {
         i0 = ig[i] - 1;
         i1 = ig[i + 1] - 1;
         res[i] = di[i] * vect[i];
         for (int k = i0; k < i1; k++)
         {
            j = jg[k];
            res[i] += ggl[k] * vect[j - 1];
            res[j - 1] += ggu[k] * vect[i];
         }
      }
   }

   void Difur::Direct(vector<double>& vec, vector<double>& res)// прямой ход
   {
      for (int i = 0; i < Nn; i++)
      {
         res[i] = vec[i];
         double s = 0;
         int j;
         for (int k = ig[i] - 1; k < ig[i + 1] - 1; k++)
         {
            j = jg[k] - 1;
            s += res[j] * gl[k];
         }
         res[i] -= s;
         res[i] /= d[i];
      }
   }

   void Difur::Reverse(vector<double>& vec, vector<double>& res)//обратный ход
   {
      for (int i = 0; i < Nn; i++)
         res[i] = vec[i];
      for (int i = Nn - 1; i >= 0; i--)
      {
         for (int k = ig[i + 1] - 2; k >= ig[i] - 1; k--)
         {
            int j = jg[k] - 1;
            res[j] -= res[i] * gu[k];
         }
      }
   }

   double norm_vector(int n, vector<double>& vec)
   {
      double sum = 0;
      for (int i = 0; i < n; i++)
         sum += vec[i] * vec[i];
      return sqrt(sum);
   }

   double scal_prod(vector<double>& a, vector<double>& b, int n)
   {
      double sum = 0;
      for (int i = 0; i < n; i++)
         sum += a[i] * b[i];
      return sum;
   }

   void mult_coef(vector<double>& vec, double k, vector<double>& res, int n)
   {
      for (int i = 0; i < n; i++)
         res[i] = k * vec[i];
   }

   void sum_vector(vector<double>& a, vector<double>& b, vector<double>& res, int n)
   {
      for (int i = 0; i < n; i++)
         res[i] = a[i] + b[i];
   }


   void Difur::LOS()
   {
      vector<double> r, p, z, LAUr, L, res;
      r.resize(Nn);
      p.resize(Nn);
      z.resize(Nn);
      LAUr.resize(Nn);
      L.resize(Nn);
      res.resize(Nn);
      double a, b, norm;
      int k;
      mult_A(P, res);
      for (int i = 0; i < Nn; i++)
         r[i] = B[i] - res[i];
      Direct(r, r);
      Reverse(r, z);
      mult_A(z, res);
      Direct(res, p);
      double normV = norm_vector(Nn, B);
      for (k = 0; (norm_vector(Nn, r) / normV) > eps && k < max_iter; k++)
      {
         norm = scal_prod(p, p, Nn);
         a = scal_prod(p, r, Nn) / norm;
         mult_coef(z, a, res, Nn);
         sum_vector(P, res, P, Nn);
         mult_coef(p, -a, res, Nn);
         sum_vector(r, res, r, Nn);
         Reverse(r, res);//res=U_-1*r
         mult_A(res, LAUr);
         Direct(LAUr, LAUr);
         b = -scal_prod(p, LAUr, Nn) / norm;
         mult_coef(z, b, z, Nn);
         sum_vector(res, z, z, Nn);
         mult_coef(p, b, p, Nn);
         sum_vector(p, LAUr, p, Nn);
      }
   }
}