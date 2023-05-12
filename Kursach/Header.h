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
      vector<vector<double>> b2, elems; // ������ ������� ������� (����� ��������, ��� ��������� ����� � ��������), �������� ��������� � ������
      vector<pair<double, double>> XY, mat, b1; // ���������� ����� (�,�), �������� (������, �����), ������ ������� ������� (������ ����� ����, �������� �������)
      vector<double> hx, hy, P, d, gl, gu, di, ggl, ggu, B;
      vector<int> ig, jg;
      double x1, x2, y1, y2, kx, ky, eps = 1e-30; // ������� �� X � Y, ������������ ����������� �� X � Y, �������
      int Nx, Ny, Ne, Nn, max_iter = 10000; // ���-�� ��������� �� X � �� Y, ���-�� ���������, ���������� �����, ������ ���-�� ��������
      int NElems, ver, hor, Nperf; //��� ��
      double Perf1, Perf2; // ���� ����������
      bool f = false;
   public:
      void Grid_Generator(); // ���������� �����
      void PD_Grid_Generator(); // ���������� �����
      void Read(); // ��������� � �������������� ������
      void Portrait(); // ��� �������� �������� ig
      void GaMbo(); // ���������� �������
      double Tetta(int n, double x, double y); // ���� ������ ��
      double Func(double x, double y) //������� ������ �����
      {
         return x*x-2;
      }
      void Boundary_conditions(); // ���� ������� �������
      void LU();
      void LOS();
      void output();
      double P_in_point(double x, double y);

      // ���� ��������������� ������� ��� ������� ���� ����� LU+���
      void mult_A(vector<double>& vect, vector<double>& res);
      void Direct(vector<double>& vec, vector<double>& res);
      void Reverse(vector<double>& vec, vector<double>& res);
   };
}