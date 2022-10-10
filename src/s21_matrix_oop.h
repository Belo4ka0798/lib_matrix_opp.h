#ifndef __SRC_S21_MATRIX_OOP_H__
#define __SRC_S21_MATRIX_OOP_H__

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
private:
  int rows_, cols_;
  double **matrix_;
  void CreateMatrix(S21Matrix &mat);
  void DeleteMatrix(S21Matrix &mat);
  void CopyMatrix(const S21Matrix &mat);
  S21Matrix DetDel(int row, int col, S21Matrix &A);

public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);

  ~S21Matrix();

  int get_cols() const;
  int get_rows() const;
  void set_rows(int rows);
  void set_cols(int cols);

  bool operator==(const S21Matrix &other);
  double &operator()(int row, int col);
  const double &operator()(int row, int col) const;
  S21Matrix operator+(const S21Matrix &other) const;
  S21Matrix operator-(const S21Matrix &other) const;
  S21Matrix operator*(const S21Matrix &other) const;
  S21Matrix operator*(const double num) const;
  friend S21Matrix operator*(const double num, const S21Matrix &other);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);

  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
};

#endif // __SRC_S21_MATRIX_OOP_H__
