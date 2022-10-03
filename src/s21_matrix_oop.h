#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
private:
  unsigned int rows_, cols_;
  double **matrix_;
  void CreateMatrix(S21Matrix &mat) const;
  void DeleteMatrix(S21Matrix &mat) const;
  void CopyMatrix(const S21Matrix &mat);
  S21Matrix DetDel(unsigned int row, unsigned int col, S21Matrix &A);

public:
  // constructor

  S21Matrix();
  S21Matrix(unsigned int rows, unsigned int cols); // size constructor
  S21Matrix(const S21Matrix &other);               // copy constructor
  S21Matrix(S21Matrix &&other);                    // perenos constructor

  // destructor

  ~S21Matrix();

  // operators

  bool operator==(const S21Matrix &other);
  double &operator()(unsigned int row, unsigned int col);
  const double &operator()(unsigned int row, unsigned int col) const;
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(const double num, const S21Matrix &other);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);

  // get and set

  unsigned int GetRows() { return rows_; }
  unsigned int GetCols() { return cols_; }
  void SetRows(unsigned int rows);
  void SetCols(unsigned int cols);

  // operations

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
};

#endif // SRC_S21_MATRIX_OOP_H_
