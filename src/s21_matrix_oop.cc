#include "./s21_matrix_oop.h"

// Constructors

S21Matrix::S21Matrix() : rows_(3), cols_(3) { CreateMatrix(*this); }

S21Matrix::S21Matrix(unsigned int rows, unsigned int cols)
    : rows_(rows), cols_(cols) {
  if (rows_ < 1 || cols_ < 1) {
    throw std::out_of_range("Incorrect input! (rows or cols < 1)");
  }
  CreateMatrix(*this);
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  // CreateMatrix(*this);
  // for (unsigned int i = 0; i < other.rows_; i++) {
  //     for (unsigned int j = 0; j < other.cols_; j++) {
  //         matrix_[i][j] = other.matrix_[i][j];
  //     }
  // }
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

// Destructor

S21Matrix::~S21Matrix() { DeleteMatrix(*this); }

// mutators

void S21Matrix::SetCols(unsigned int col) {
  if (col < 1) {
    throw std::logic_error("Incorrect input! (col < 1)");
  }
  unsigned int helpcol = col;
  if (col > cols_) {
    helpcol = cols_;
  }
  for (unsigned int i = 0; i < rows_; i++) {
    double *helper = new double[col]{}; 
    for (unsigned int j = 0; j < helpcol; j++) {
      helper[j] = matrix_[i][j];
    }
    delete[] matrix_[i];
    matrix_[i] = helper;
  }
  cols_ = col;
}

void S21Matrix::SetRows(unsigned int row) {
  if (row < 1)
    throw std::logic_error("Incorrect input! (row < 1)");
  
  unsigned int helprow = row;
  if (row > rows_) {
    helprow = rows_;
  }
  double **helper = new double *[row]();
  for (unsigned int i = 0; i < row; i++) {
    helper[i] = new double[cols_]();
  }
   
  for (unsigned int i = 0; i < helprow; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      helper[i][j] = matrix_[i][j];
    }
  }
  DeleteMatrix(*this);
  matrix_ = helper;
  rows_ = row;
}

// Operations

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error(
        "Incorrect input, matrices should have the same size");
  }

  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-07) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error(
        "Incorrect input, matrices should have the same size");
  }
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error(
        "Incorrect input, matrices should have the same size");
  }
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "Incorrect input, matrices should have the same size");
  }
  S21Matrix result{rows_, other.cols_};
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < other.cols_; j++) {
      for (unsigned int m = 0; m < cols_; m++) {
        result(i, j) += matrix_[i][m] * other.matrix_[m][j];
      }
    }
  }
  *this = result;
}

S21Matrix S21Matrix::Transpose() {
  // leaks
  S21Matrix res{cols_, rows_};
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::out_of_range("cols != rows");
  S21Matrix result{rows_, cols_};
  S21Matrix del{rows_,cols_};
  double res = 0;
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      // leaks
      del = DetDel(i, j, *this);
      res = del.Determinant();
      result(i, j) = pow(-1, i + j) * res;
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (cols_ != rows_)
    throw std::out_of_range("cols != rows");
  double result = 0;
  if (cols_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else if (cols_ == 1) {
    result = matrix_[0][0];
  } else {
    S21Matrix del;
    double res = 0;
    for (unsigned int j = 0; j < cols_; j++) {
      del = DetDel(0, j, *this);
      res = del.Determinant();
      result += pow(-1, j) * matrix_[0][j] * res;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double opr = this->Determinant();
  if (fabs(opr) < 1e-7)
    throw std::out_of_range("Determinant = 0");
  if (rows_ < 2 || cols_ < 2)
    throw std::out_of_range("Incorrect matrix");
  S21Matrix result = this->CalcComplements().Transpose();
  result.MulNumber(1.0 / opr);
  return result;
}

// Operators

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

double &S21Matrix::operator()(unsigned int row, unsigned int col) {
  if (row >= rows_ || col >= cols_) {
    throw std::out_of_range("Incorrect input! (row/col >= max row/col)");
  }
  return matrix_[row][col];
}

const double &S21Matrix::operator()(unsigned int row, unsigned int col) const {
  if (row >= rows_ || col >= cols_) {
    throw std::out_of_range("Incorrect input! (row/col >= max row/col)");
  }
  return matrix_[row][col];
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result_mat(*this);
  result_mat.SumMatrix(other);
  return result_mat;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result_mat(*this);
  result_mat.SubMatrix(other);
  return result_mat;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result_mat(*this);
  result_mat.MulMatrix(other);
  return result_mat;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result_mat(*this);
  result_mat.MulNumber(num);
  return result_mat;
}

S21Matrix operator*(const double num, const S21Matrix &other) {
  S21Matrix result_mat(other);
  result_mat.MulNumber(num);
  return result_mat;
}

// DRY
S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    DeleteMatrix(*this);
    rows_ = other.rows_;
    cols_ = other.cols_;
    // CreateMatrix(*this);
    // for (unsigned int i = 0; i < rows_; i++) {
    //     for (unsigned int j = 0; j < cols_; j++) {
    //         matrix_[i][j] = other.matrix_[i][j];
    //     }
    // }
    CopyMatrix(other);
  }
  return *this;
}

// DRY
S21Matrix &S21Matrix::operator=(S21Matrix &&other) {
  if (this != &other) {
    DeleteMatrix(*this);
    rows_ = other.rows_;
    cols_ = other.cols_;
    // CreateMatrix(*this);
    // for (unsigned int i = 0; i < rows_; i++) {
    //     for (unsigned int j = 0; j < cols_; j++) {
    //         matrix_[i][j] = other.matrix_[i][j];
    //     }
    // }
    CopyMatrix(other);
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

// support

void S21Matrix::CreateMatrix(S21Matrix &mat) const {
  mat.matrix_ = new double *[rows_]{};
  for (unsigned int i = 0; i < rows_; i++) {
    mat.matrix_[i] = new double[cols_]{};
  }
}

void S21Matrix::DeleteMatrix(S21Matrix &mat) const {
  for (unsigned int i = 0; i < rows_; i++) {
    delete[] mat.matrix_[i];
  }
  delete[] mat.matrix_;
}

void S21Matrix::CopyMatrix(const S21Matrix &mat) {
  CreateMatrix(*this);
  for (unsigned int i = 0; i < rows_; i++) {
    for (unsigned int j = 0; j < cols_; j++) {
      matrix_[i][j] = mat.matrix_[i][j];
    }
  }
}

// leaks
S21Matrix S21Matrix::DetDel(unsigned int row, unsigned int col, S21Matrix &A) {
  unsigned int n = 0, m = 0;
  S21Matrix result{A.rows_ - 1, A.cols_ - 1};
  for (unsigned int i = 0; i < A.rows_; i++) {
    for (unsigned int j = 0; j < A.cols_; j++) {
      if (i != row && j != col) {
        result(n, m) = A.matrix_[i][j];
        m++;
        if (m == A.cols_ - 1) {
          n++;
          m = 0;
        }
      }
    }
  }
  return result;
}
