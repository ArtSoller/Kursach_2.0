#include"Header.h"
using namespace solution;

void showmenu()
{
   cout << "1 - построение сетки и вывод координат" << endl
      << "2 - подготовка начального приближения" << endl
      << "Не целочисленное значение - выход" << endl;
}

int main()
{
   setlocale(LC_ALL, "Russian");
   Difur object;
   object.Grid_Generator();// построение сетки и вывод координат
   //object.PD_Grid_Generator();// построение сетки и вывод координат
   object.Read();// ввод данных
   object.Portrait(); // список смежности элементов (портрет)?
   object.GaMbo();
   object.Boundary_conditions();
   object.LU();
   object.LOS();
   object.output();
   return 0;
}