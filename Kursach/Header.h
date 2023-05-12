#pragma once
#include<iostream>
#include<fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

namespace solution
{
   class Difur
   {
   private:
      vector<vector<double>> b2, elems; // второе краевое условие (номер элемента, его локальная грань и значение), элементы разбиения с узлами
      vector<pair<double, double>> XY, mat, b1; // координаты узлов (х,у), материал (лямбда, гамма), первое краевое условие (глобал номер узла, значение функции)
      vector<double> hx, hy, P, d, gl, gu, di, ggl, ggu, B;
      vector<int> ig, jg;
      double x1, x2, y1, y2, kx, ky, eps = 1e-30; // границы по X и Y, коэффициенты релаксакции по X и Y, эпсилон
      int Nx, Ny, Ne, Nn, max_iter = 10000; // кол-во разбиений по X и по Y, кол-во элементов, количество узлов, максим кол-во итераций
      int NElems, ver, hor, Nperf; //для ПД
      double Perf1, Perf2; // зона перфорации
      bool f = false;
   public:
      void Grid_Generator(); // составляем сетку
      void PD_Grid_Generator(); // составляем сетку
      void Read(); // считываем и подготавливаем данные
      void Portrait(); // для создания векторов ig
      void GaMbo(); // построение матрицы
      double Tetta(int n, double x, double y); // пока вообще ХЗ
      double Func(double x, double y) //функция правой части
      {
         return x*x-2;
      }
      void Boundary_conditions(); // учет краевых условий
      void LU();
      void LOS();
      void output();
      double P_in_point(double x, double y);

      // Блок вспомагательных функций для решения СЛАУ через LU+ЛОС
      void mult_A(vector<double>& vect, vector<double>& res);
      void Direct(vector<double>& vec, vector<double>& res);
      void Reverse(vector<double>& vec, vector<double>& res);
   };
}