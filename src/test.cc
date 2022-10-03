#include "s21_matrix_oop.h"

using namespace std;
void PrintMatrix(S21Matrix mat);

int main() {
  int rows = 3;
  int cols = 3;
  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  given(0, 0) = 1.0;
  given(0, 1) = 2.0;
  given(0, 2) = 3.0;
  given(1, 0) = 4.0;
  given(1, 1) = 5.0;
  given(1, 2) = 6.0;
  given(2, 0) = 7.0;
  given(2, 1) = 8.0;
  given(2, 2) = 9.0;
  given.SetRows(4);
  cout << "Given after SetRows(4)" << endl;
  for (unsigned int i = 0; i < given.GetRows(); i++) {
    for (unsigned int j = 0; j < given.GetCols(); j++) {
      cout << given(i,j) << " ";
    }
      cout << endl;
  }
  given.SetRows(2);
  given.SetCols(2);
  cout << "Given after SetRows(2) SetCols(2)" << endl;
  for (unsigned int i = 0; i < given.GetRows(); i++) {
    for (unsigned int j = 0; j < given.GetCols(); j++) {
      cout << given(i,j) << " ";
    }
      cout << endl;
  }

  expected(0, 0) = 10.0;
  expected(0, 1) = 11.0;
  expected(0, 2) = 12.0;
  expected(1, 0) = 13.0;
  expected(1, 1) = 14.0;
  expected(1, 2) = 15.0;
  expected(2, 0) = 16.0;
  expected(2, 1) = 17.0;
  expected(2, 2) = 18.0;
  expected.SetCols(4);
  cout << "Expected after SetCols(4)" << endl;
  for (unsigned int i = 0; i < expected.GetRows(); i++) {
    for (unsigned int j = 0; j < expected.GetCols(); j++) {
      cout << expected(i,j) << " ";
    }
      cout << endl;
  }
  expected.SetCols(2);
  expected.SetRows(2);
  cout << "Expected after SetRows(2) SetCols(2)" << endl;
  for (unsigned int i = 0; i < expected.GetRows(); i++) {
    for (unsigned int j = 0; j < expected.GetCols(); j++) {
      cout << expected(i,j) << " ";
    }
      cout << endl;
  }
  
  cout << endl;
  cout << "Operators" << endl;
  cout << endl;
  
  cout << "Oper_Sum" << endl;
  S21Matrix oper_Sum = given + expected;
  PrintMatrix(oper_Sum);

  cout << "Oper_Sub" << endl;
  S21Matrix oper_Sub = given - expected;
  PrintMatrix(oper_Sub);

  cout << "Oper_MulNM" << endl;
  S21Matrix oper_MulNM = 10 * given;
  PrintMatrix(oper_MulNM);

  cout << "Oper_MulMN" << endl;
  S21Matrix oper_MulMN = given * 10;
  PrintMatrix(oper_MulMN);

  cout << "Oper_MulM" << endl;
  S21Matrix oper_MulM = given * expected;
  PrintMatrix(oper_MulM);
  
  cout << endl;
  cout << "Operators=" << endl;
  cout << endl;

  cout << "Sum" << endl;
  S21Matrix Sum(2,2); Sum += given;
  PrintMatrix(Sum);

  cout << "Sub" << endl;
  S21Matrix Sub(2,2); Sub -= given;
  PrintMatrix(Sub);

  cout << "MulN" << endl;
  S21Matrix MulN(2,2); MulN *= 10;
  PrintMatrix(MulN);

  cout << "MulM" << endl;
  S21Matrix MulM(2,2); MulM *= given;
  PrintMatrix(MulM);

  cout << endl;
  cout << "Func Matrix" << endl;
  cout << endl;

  cout << "Func Sum" << endl;
  S21Matrix func_sum(2,2); func_sum.SumMatrix(given);
  PrintMatrix(func_sum);

  cout << "Func Sub" << endl;
  S21Matrix func_sub(2,2); func_sub.SubMatrix(given);
  PrintMatrix(func_sub);

  cout << "Func MulNum" << endl;
  S21Matrix func_muln(2,2); func_muln.MulNumber(10);
  PrintMatrix(func_muln);

  cout << "Func MulMatrix" << endl;
  S21Matrix func_mulm(2,2); func_mulm.MulMatrix(given);
  PrintMatrix(func_mulm);

  cout << "Func Transpose" << endl;
  S21Matrix func_trans(2,2); func_trans = given.Transpose();
  PrintMatrix(func_trans);
}

void PrintMatrix(S21Matrix mat) {
  for (unsigned int i = 0; i < mat.GetRows(); i++) {
    for (unsigned int j = 0; j < mat.GetCols(); j++) {
      cout << mat(i,j) << " ";
    }
      cout << endl;
  }
}