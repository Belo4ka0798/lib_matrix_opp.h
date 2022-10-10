#include "s21_matrix_oop.h"
#include <cstdlib>
#include <gtest/gtest.h>

void randm(S21Matrix &m) {
  for (int i = 0; i < m.get_rows(); i++)
    for (int j = 0; j < m.get_cols(); j++)
      m(i, j) = rand() % 10;
}

void printm(S21Matrix &m) {
  for (int i = 0; i < m.get_rows(); i++) {
    for (int j = 0; j < m.get_cols(); j++)
      std::cout << m(i, j) << " ";
    std::cout << std::endl;
  }
}

TEST(test, EqMatrix_1) {
  S21Matrix josky, flex;
  for (int i = 0; i < josky.get_rows(); i++) {
    for (int j = 0; j < josky.get_cols(); j++) {
      josky(i, j) = 1;
      flex(i, j) = 1;
    }
  }
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, EqMatrix_2) {
  S21Matrix josky, flex;
  josky(2, 2) = 1;
  EXPECT_FALSE(josky.EqMatrix(flex));
}

TEST(test, SumMatrix_1) {
  S21Matrix josky(2, 2);
  josky(0, 0) = 2;
  josky(0, 1) = 4;
  josky(1, 0) = 6;
  josky(1, 1) = 8;
  S21Matrix flex = josky * 2;
  josky.SumMatrix(josky);
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, SumMatrix_2) {
  S21Matrix josky(2, 3), flex(2, 2);

  EXPECT_THROW(josky.SumMatrix(flex);, std::logic_error);
}

TEST(test, SubMatrix_1) {
  S21Matrix josky(2, 2), flex(2, 2);
  josky(0, 0) = 2;
  josky(0, 1) = 4;
  josky(1, 0) = 6;
  josky(1, 1) = 8;
  josky.SubMatrix(josky);
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, SubMatrix_2) {
  S21Matrix josky(2, 2), flex(2, 2);
  josky(0, 0) = 2;
  josky(0, 1) = 4;
  josky(1, 0) = 6;
  josky(1, 1) = 8;
  flex.SubMatrix(josky);
  josky *= -1;
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, MulNumber_1) {
  S21Matrix josky(2, 2), flex(2, 2);
  josky(0, 0) = 2;
  josky(0, 1) = 4;
  josky(1, 0) = 6;
  josky(1, 1) = 8;
  flex = josky;
  josky.MulNumber(2);
  flex += flex;
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, MulNumber_2) {
  S21Matrix josky(2, 2), flex(2, 2);
  josky(0, 0) = 2;
  josky(0, 1) = 4;
  josky(1, 0) = 6;
  josky(1, 1) = 8;
  josky.MulNumber(0);
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, MulMatrix_1) {
  S21Matrix josky(2, 3), flex(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      flex(j, i) = count;
      count++;
    }
  }
  josky.MulMatrix(flex);
  S21Matrix res(2, 2);
  res(0, 0) = 14;
  res(0, 1) = 32;
  res(1, 0) = 32;
  res(1, 1) = 77;
  EXPECT_TRUE(josky.EqMatrix(res));
}

TEST(test, MulMatrix_2) {
  S21Matrix josky(2, 3), flex(4, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      flex(j, i) = count;
      count++;
    }
  }
  EXPECT_THROW(josky.MulMatrix(flex), std::logic_error);
}

TEST(test, Transpose_1) {
  S21Matrix josky(2, 3), flex(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      flex(j, i) = count;
      count++;
    }
  }
  josky = josky.Transpose();
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, Transpose_2) {
  S21Matrix josky(2, 3), flex(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      flex(j, i) = count;
      count++;
    }
  }
  flex = flex.Transpose();
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, CalcComplements_1) {
  S21Matrix josky(2, 2), flex(2, 2);
  int count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      josky(i, j) = count;
      count++;
    }
  }
  josky = josky.CalcComplements();
  flex(0, 0) = 4;
  flex(0, 1) = -3;
  flex(1, 0) = -2;
  flex(1, 1) = 1;
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, Determinant_1) {
  S21Matrix josky;
  int count = 1;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      count++;
    }
  }
  EXPECT_DOUBLE_EQ(0, josky.Determinant());
}

TEST(test, Determinant_2) {
  S21Matrix josky;
  int count = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      josky(i, j) = count;
      count += i + j;
    }
    count *= 2;
  }
  EXPECT_DOUBLE_EQ(-12, josky.Determinant());
}

TEST(test, Inverse_1) {
  S21Matrix josky(2, 2), flex(2, 2);
  int count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      josky(i, j) = count;
      count++;
    }
  }
  josky = josky.InverseMatrix();
  flex(0, 0) = -2;
  flex(0, 1) = 1;
  flex(1, 0) = 1.5;
  flex(1, 1) = -0.5;
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(test, Inverse_2) {
  S21Matrix josky(2, 2), flex(2, 2);
  josky(0, 0) = 1;
  josky(0, 1) = 1;
  josky(1, 0) = 3;
  josky(1, 1) = 2;
  josky = josky.InverseMatrix();
  flex(0, 0) = -2;
  flex(0, 1) = 1;
  flex(1, 0) = 3;
  flex(1, 1) = -1;
  EXPECT_TRUE(josky.EqMatrix(flex));
}

TEST(get_set, t1) {
  S21Matrix m1(2, 2);
  m1(1, 1) = 123.456;
  EXPECT_EQ(m1.get_cols(), 2);
  EXPECT_EQ(m1.get_rows(), 2);
  EXPECT_DOUBLE_EQ(m1(1, 1), 123.456);
}

TEST(set_col, t1) {
  S21Matrix m(2, 2);
  m(1, 1) = 123.456;
  m.set_cols(3);
  m(1, 2) = 999.999;
  EXPECT_EQ(m.get_cols(), 3);
  EXPECT_DOUBLE_EQ(m(1, 1), 123.456);
  EXPECT_DOUBLE_EQ(m(1, 2), 999.999);
}

TEST(set_row, t1) {
  S21Matrix m1(2, 2);
  m1(1, 0) = 123.456;
  m1.set_rows(3);
  m1(2, 0) = m1(1, 0);
  EXPECT_EQ(m1.get_rows(), 3);
  EXPECT_EQ(m1.get_cols(), 2);
  EXPECT_DOUBLE_EQ(m1(1, 0), m1(2, 0));
}

TEST(constructors, t1) {
  S21Matrix m1;
  EXPECT_EQ(m1.get_cols(), 3);
  EXPECT_EQ(m1.get_rows(), 3);

  S21Matrix m2(5, 5);
  m2(1, 1) = 123.456;
  EXPECT_EQ(m2.get_cols(), 5);
  EXPECT_EQ(m2.get_rows(), 5);

  S21Matrix m3(m2);
  EXPECT_EQ(m2 == m3, true);

  S21Matrix m4(std::move(m3));
  EXPECT_EQ(m3.get_cols(), 0);
  EXPECT_EQ(m3.get_rows(), 0);
  EXPECT_EQ(m4.get_rows(), 5);
  EXPECT_EQ(m4.get_cols(), 5);
  EXPECT_EQ(m4(1, 1), m2(1, 1));
}

TEST(func1, eq) {
  S21Matrix m1(123, 123);
  S21Matrix m2(m1);

  EXPECT_EQ(m1.EqMatrix(m2), true);
}

TEST(func2, summat) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);

  m1(0, 0) = 1;
  m1(0, 2) = 2;
  m1(2, 0) = 3;

  m2.SumMatrix(m1);

  EXPECT_EQ(m2(0, 0), 1);
  EXPECT_EQ(m2(0, 2), 2);
  EXPECT_EQ(m2(2, 0), 3);
}

TEST(func3, submat) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);

  m1(0, 0) = 1;
  m1(0, 2) = 2;
  m1(2, 0) = 3;

  m2.SubMatrix(m1);

  EXPECT_EQ(m2(0, 0), -1);
  EXPECT_EQ(m2(0, 2), -2);
  EXPECT_EQ(m2(2, 0), -3);
}

TEST(func4, mulnum) {
  S21Matrix m2(3, 3);

  m2(0, 0) = 1.3;
  m2(0, 2) = 2.1;
  m2(2, 0) = 3.3;

  m2.MulNumber(0.412);

  EXPECT_EQ(m2(0, 0), 1.3 * 0.412);
  EXPECT_EQ(m2(0, 2), 2.1 * 0.412);
  EXPECT_EQ(m2(2, 0), 3.3 * 0.412);
}

TEST(func5, mulmat) {
  S21Matrix m1(2, 1);
  S21Matrix m2(1, 2);

  m1(0, 0) = 32;
  m1(1, 0) = 12;

  m2(0, 0) = 2;
  m2(0, 1) = 3;

  m1.MulMatrix(m2);

  EXPECT_DOUBLE_EQ(m1(0, 0), 32 * 2);
  EXPECT_DOUBLE_EQ(m1(1, 1), 12 * 3);
}

TEST(func6, transpose) {
  int rows = 2;
  int cols = 3;

  S21Matrix m(rows, cols);
  S21Matrix res(cols, rows);

  for (int i = 0, c = 1; i < rows; i++) {
    for (int j = 0; j < cols; j++, c++) {
      m(i, j) = c;
      res(j, i) = c;
    }
  }

  m = m.Transpose();

  ASSERT_TRUE(m == res);
}

TEST(func7, calccal) {
  S21Matrix m(3, 12);
  EXPECT_ANY_THROW(m.CalcComplements());
}

TEST(func7, calcalc) {
  int rows = 3;
  int cols = 3;

  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  given(0, 0) = 1.0;
  given(0, 1) = 2.0;
  given(0, 2) = 3.0;
  given(1, 0) = 0.0;
  given(1, 1) = 4.0;
  given(1, 2) = 2.0;
  given(2, 0) = 5.0;
  given(2, 1) = 2.0;
  given(2, 2) = 1.0;

  expected(0, 0) = 0.0;
  expected(0, 1) = 10.0;
  expected(0, 2) = -20.0;
  expected(1, 0) = 4.0;
  expected(1, 1) = -14.0;
  expected(1, 2) = 8.0;
  expected(2, 0) = -8.0;
  expected(2, 1) = -2.0;
  expected(2, 2) = 4.0;

  S21Matrix res = given.CalcComplements();

  ASSERT_TRUE(res == expected);
}

TEST(func8, det) {
  S21Matrix m(5, 5);

  for (uint32_t i = 0; i < 5; i++)
    for (uint32_t j = 0; j < 5; j++)
      m(i, j) = j;

  ASSERT_NEAR(0, m.Determinant(), 1e-06);
}

TEST(func8, det1) {
  int size = 5;
  S21Matrix m(size, size);

  m(0, 1) = 6;
  m(0, 2) = -2;
  m(0, 3) = -1;
  m(0, 4) = 5;
  m(1, 3) = -9;
  m(1, 4) = -7;
  m(2, 1) = 15;
  m(2, 2) = 35;
  m(3, 1) = -1;
  m(3, 2) = -11;
  m(3, 3) = -2;
  m(3, 4) = 1;
  m(4, 0) = -2;
  m(4, 1) = -2;
  m(4, 2) = 3;
  m(4, 4) = -2;

  double res = m.Determinant();
  ASSERT_NEAR(res, 2480, 1e-6);
}

TEST(func8, det2) {
  int size = 3;
  S21Matrix m(size, size);

  m(0, 0) = 2;
  m(0, 1) = 3;
  m(0, 2) = 1;
  m(1, 0) = 7;
  m(1, 1) = 4;
  m(1, 2) = 1;
  m(2, 0) = 9;
  m(2, 1) = -2;
  m(2, 2) = 1;

  double res = m.Determinant();
  ASSERT_NEAR(res, -32, 1e-6);
}

TEST(func9, inverse) {
  int size = 3;
  S21Matrix given(size, size);
  S21Matrix expected(size, size);

  expected(0, 0) = 44300.0 / 367429.0;
  expected(0, 1) = -236300.0 / 367429.0;
  expected(0, 2) = 200360.0 / 367429.0;
  expected(1, 0) = 20600.0 / 367429.0;
  expected(1, 1) = 56000.0 / 367429.0;
  expected(1, 2) = -156483.0 / 367429.0;
  expected(2, 0) = 30900.0 / 367429.0;
  expected(2, 1) = 84000.0 / 367429.0;
  expected(2, 2) = -51010.0 / 367429.0;

  given(0, 0) = 2.8;
  given(0, 1) = 1.3;
  given(0, 2) = 7.01;
  given(1, 0) = -1.03;
  given(1, 1) = -2.3;
  given(1, 2) = 3.01;
  given(2, 0) = 0;
  given(2, 1) = -3;
  given(2, 2) = 2;

  ASSERT_TRUE(given.InverseMatrix() == expected);
}

TEST(func9, iverse1) {
  int size = 3;
  S21Matrix given(size, size);
  S21Matrix expected(size, size);
  expected(0, 0) = 1.0;
  expected(0, 1) = -1.0;
  expected(0, 2) = 1.0;
  expected(1, 0) = -38.0;
  expected(1, 1) = 41.0;
  expected(1, 2) = -34.0;
  expected(2, 0) = 27.0;
  expected(2, 1) = -29.0;
  expected(2, 2) = 24.0;

  given(0, 0) = 2.0;
  given(0, 1) = 5.0;
  given(0, 2) = 7.0;
  given(1, 0) = 6.0;
  given(1, 1) = 3.0;
  given(1, 2) = 4.0;
  given(2, 0) = 5.0;
  given(2, 1) = -2.0;
  given(2, 2) = -3.0;
  ASSERT_TRUE(given.InverseMatrix() == expected);
}

TEST(method1, first) {
  S21Matrix m;
  EXPECT_EQ(m.get_cols(), 3);
  EXPECT_EQ(m.get_rows(), 3);
}

TEST(method2, second) {
  S21Matrix m(11, 11);
  EXPECT_EQ(m.get_cols(), 11);
  EXPECT_EQ(m.get_rows(), 11);
  EXPECT_EQ(m(10, 10), 0);
}

TEST(operator, mul) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);

  for (int i = 0; i < m1.get_rows(); ++i)
    for (int j = 0; j < m1.get_cols(); ++j) {
      m1(i, j) = (i + 1) * j + 1;
      m2(i, j) = (i + 1) * j + 1;
    }

  m1.MulMatrix(m2);
  EXPECT_EQ(3, m1.get_cols());
  EXPECT_EQ(3, m1.get_rows());
  for (int i = 0; i < m1.get_rows(); ++i)
    for (int j = 0; j < m1.get_cols(); ++j)
      if (j == 0)
        EXPECT_DOUBLE_EQ(6 + 3 * i, m1(i, j));
      else if (j == 1)
        EXPECT_DOUBLE_EQ(20 + 11 * i, m1(i, j));
      else
        EXPECT_DOUBLE_EQ(34 + 19 * i, m1(i, j));

  S21Matrix m3(3, 1);
  ASSERT_ANY_THROW(m3.MulMatrix(m2));

  for (int i = 0; i < m3.get_rows(); ++i)
    m3(i, 0) = 2;
  m2.MulMatrix(m3);
  EXPECT_EQ(1, m2.get_cols());
  EXPECT_EQ(3, m2.get_rows());
  for (int i = 0; i < m2.get_rows(); ++i)
    EXPECT_DOUBLE_EQ(12 + 6 * i, m2(i, 0));
}

TEST(method3, third) {
  S21Matrix m1(123, 123);

  S21Matrix m2(m1);

  EXPECT_EQ(m1.get_cols(), m2.get_cols());
  EXPECT_EQ(m1.get_rows(), m2.get_rows());
}

TEST(method4, fourth) {
  S21Matrix m1(123, 123);

  S21Matrix m2(std::move(m1));

  EXPECT_EQ(m1.get_cols(), 0);
  EXPECT_EQ(m1.get_rows(), 0);
  EXPECT_EQ(m2.get_cols(), 123);
  EXPECT_EQ(m2.get_rows(), 123);
}

TEST(valgr, gg) {
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
  given.set_rows(4);

  given.set_rows(2);
  given.set_cols(2);

  expected(0, 0) = 10.0;
  expected(0, 1) = 11.0;
  expected(0, 2) = 12.0;
  expected(1, 0) = 13.0;
  expected(1, 1) = 14.0;
  expected(1, 2) = 15.0;
  expected(2, 0) = 16.0;
  expected(2, 1) = 17.0;
  expected(2, 2) = 18.0;
  expected.set_cols(4);

  expected.set_cols(2);
  expected.set_rows(2);

  S21Matrix oper_Sum = given + expected;
  S21Matrix oper_Sub = given - expected;
  S21Matrix oper_MulNM = 10 * given;
  S21Matrix oper_MulMN = given * 10;
  S21Matrix oper_MulM = given * expected;
  S21Matrix Sum(2, 2);
  Sum += given;
  S21Matrix Sub(2, 2);
  Sub -= given;
  S21Matrix MulN(2, 2);
  MulN *= 10;
  S21Matrix MulM(2, 2);
  MulM *= given;
  S21Matrix func_sum(2, 2);
  func_sum.SumMatrix(given);
  S21Matrix func_sub(2, 2);
  func_sub.SubMatrix(given);
  S21Matrix func_muln(2, 2);
  func_muln.MulNumber(10);
  S21Matrix func_mulm(2, 2);
  func_mulm.MulMatrix(given);
  S21Matrix func_trans(2, 2);
  func_trans = given.Transpose();
  S21Matrix func_inverse(2, 2);
  func_inverse = given.InverseMatrix();
}

int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
  return 0;
}