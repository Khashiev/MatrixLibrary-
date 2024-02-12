#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
}

void default_matrix(S21Matrix& other) {
  int k = 1;
  for (int i = 0; i < other.GetRows(); ++i) {
    for (int j = 0; j < other.GetCols(); ++j) {
      other(i, j) = k;
      ++k;
    }
  }
}

namespace test {

TEST(Constructor, Default) {
  S21Matrix test = S21Matrix();
  EXPECT_EQ(test.GetRows(), 2);
  EXPECT_EQ(test.GetCols(), 2);
}

TEST(Constructor, with_args) {
  S21Matrix test(2, 4);
  EXPECT_EQ(test.GetRows(), 2);
  EXPECT_EQ(test.GetCols(), 4);
}

TEST(Constructor, Copy) {
  S21Matrix test(2, 4);
  test(0, 0) = 3;
  S21Matrix test2(test);
  EXPECT_EQ(test2.GetRows(), 2);
  EXPECT_EQ(test2.GetCols(), 4);
  EXPECT_EQ(test2(0, 0), 3);
}

TEST(Constructor, Move) {
  S21Matrix test(2, 4);
  test(0, 0) = 3;
  S21Matrix test2 = std::move(test);
  EXPECT_EQ(test2.GetRows(), 2);
  EXPECT_EQ(test2.GetCols(), 4);
  EXPECT_EQ(test2(0, 0), 3);
}

TEST(Set, Rows) {
  S21Matrix test(2, 2);
  test(0, 0) = 3;
  test(0, 1) = 1;
  test(1, 0) = 2;
  test(1, 1) = -2;
  EXPECT_EQ(test.GetRows(), 2);
  test.SetRows(3);
  EXPECT_EQ(test.GetRows(), 3);
  EXPECT_EQ(test(1, 1), -2);
  test.SetRows(1);
  test.SetCols(1);
  test.SetRows(2);
  test.SetCols(2);
  EXPECT_EQ(test(1, 1), 0);
}

TEST(Set, Cols) {
  S21Matrix test(2, 1);
  test(0, 0) = 3;
  test(1, 0) = 2;
  EXPECT_EQ(test.GetCols(), 1);
  test.SetCols(10);
  EXPECT_EQ(test.GetCols(), 10);
  EXPECT_EQ(test(1, 9), 0);
  EXPECT_EQ(test(0, 0), 3);
}

TEST(Arithmethic, EqMatrix) {
  S21Matrix test1(3, 3);
  S21Matrix test2(3, 3);
  default_matrix(test1);
  default_matrix(test2);
  EXPECT_EQ(test1.EqMatrix(test2), true);
  EXPECT_EQ(test2.EqMatrix(test1), true);
  EXPECT_EQ(test1 == test2, true);
  EXPECT_EQ(test2 == test1, true);
  test1(0, 0) = 10;
  EXPECT_EQ(test1.EqMatrix(test2), false);
  EXPECT_EQ(test2.EqMatrix(test1), false);
  EXPECT_EQ(test1 == test2, false);
  EXPECT_EQ(test2 == test1, false);
}

TEST(Arithmethic, SumMatrix) {
  S21Matrix test1 = S21Matrix(2, 2);
  S21Matrix test2 = S21Matrix(2, 2);
  default_matrix(test1);
  default_matrix(test2);
  test1 += test2;
  EXPECT_EQ(test1(0, 0), 2);
  EXPECT_EQ(test1(0, 1), 4);
  EXPECT_EQ(test1(1, 0), 6);
  EXPECT_EQ(test1(1, 1), 8);

  S21Matrix test3 = test1 + test2;
  EXPECT_EQ(test3(0, 0), 3);
  EXPECT_EQ(test3(0, 1), 6);
  EXPECT_EQ(test3(1, 0), 9);
  EXPECT_EQ(test3(1, 1), 12);
}

TEST(Arithmethic, SubMatrix) {
  S21Matrix test1 = S21Matrix(2, 2);
  S21Matrix test2 = S21Matrix(2, 2);
  default_matrix(test1);
  default_matrix(test2);
  test1 *= 3;
  test1 -= test2;
  EXPECT_EQ(test1(0, 0), 2);
  EXPECT_EQ(test1(0, 1), 4);
  EXPECT_EQ(test1(1, 0), 6);
  EXPECT_EQ(test1(1, 1), 8);

  S21Matrix test3 = test1 - test2;
  EXPECT_EQ(test3(0, 0), 1);
  EXPECT_EQ(test3(0, 1), 2);
  EXPECT_EQ(test3(1, 0), 3);
  EXPECT_EQ(test3(1, 1), 4);
}

TEST(Arithmethic, MulMatrix) {
  S21Matrix test1 = S21Matrix(3, 2);
  S21Matrix test2 = S21Matrix(2, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 4;
  test1(1, 0) = 2;
  test1(1, 1) = 5;
  test1(2, 0) = 3;
  test1(2, 1) = 6;

  test2(0, 0) = 1;
  test2(0, 1) = -1;
  test2(0, 2) = 1;
  test2(1, 0) = 2;
  test2(1, 1) = 3;
  test2(1, 2) = 4;

  test1 *= test2;

  EXPECT_EQ(test1(0, 0), 9);
  EXPECT_EQ(test1(0, 1), 11);
  EXPECT_EQ(test1(0, 2), 17);
  EXPECT_EQ(test1(1, 0), 12);
  EXPECT_EQ(test1(1, 1), 13);
  EXPECT_EQ(test1(1, 2), 22);
  EXPECT_EQ(test1(2, 0), 15);
  EXPECT_EQ(test1(2, 1), 15);
  EXPECT_EQ(test1(2, 2), 27);

  S21Matrix test3 = S21Matrix(2, 2);
  S21Matrix test4 = S21Matrix(2, 2);
  default_matrix(test3);
  default_matrix(test4);

  S21Matrix test5 = test3 * test4;
  EXPECT_EQ(test5(0, 0), 7);
  EXPECT_EQ(test5(0, 1), 10);
  EXPECT_EQ(test5(1, 0), 15);
  EXPECT_EQ(test5(1, 1), 22);
}

TEST(Arithmethic, MulNum) {
  S21Matrix test1 = S21Matrix(2, 2);
  S21Matrix test2 = S21Matrix(2, 2);
  default_matrix(test1);
  default_matrix(test2);
  S21Matrix test3 = 3 * test1;
  S21Matrix test4 = test2 * 3;
  EXPECT_EQ(test3(0, 0), 3);
  EXPECT_EQ(test3(0, 1), 6);
  EXPECT_EQ(test3(1, 0), 9);
  EXPECT_EQ(test3(1, 1), 12);
  EXPECT_EQ(test4(0, 0), 3);
  EXPECT_EQ(test4(0, 1), 6);
  EXPECT_EQ(test4(1, 0), 9);
  EXPECT_EQ(test4(1, 1), 12);
}

TEST(Calc, Transpose) {
  S21Matrix test1 = S21Matrix(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 20;
  test1(1, 1) = -2;
  test1 = test1.Transpose();
  EXPECT_EQ(test1(0, 0), 3);
  EXPECT_EQ(test1(0, 1), 20);
  EXPECT_EQ(test1(1, 0), 1);
  EXPECT_EQ(test1(1, 1), -2);
}

TEST(Calc, Complements) {
  S21Matrix test1(3, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 0;
  test1(1, 1) = 4;
  test1(1, 2) = 2;
  test1(2, 0) = 5;
  test1(2, 1) = 2;
  test1(2, 2) = 1;

  test1 = test1.CalcComplements();

  EXPECT_EQ(test1(0, 0), 0);
  EXPECT_EQ(test1(0, 1), 10);
  EXPECT_EQ(test1(0, 2), -20);
  EXPECT_EQ(test1(1, 0), 4);
  EXPECT_EQ(test1(1, 1), -14);
  EXPECT_EQ(test1(1, 2), 8);
  EXPECT_EQ(test1(2, 0), -8);
  EXPECT_EQ(test1(2, 1), -2);
  EXPECT_EQ(test1(2, 2), 4);
}

TEST(Calc, Determinant) {
  S21Matrix test1(2, 2);
  test1(0, 0) = 3;
  test1(0, 1) = 1;
  test1(1, 0) = 20;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), -26);
  test1(0, 0) = 2;
  test1(0, 1) = 1;
  test1(1, 0) = 4;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), -8);
  test1(0, 0) = 2;
  test1(0, 1) = -5;
  test1(1, 0) = 1;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), 1);
  test1(0, 0) = 3;
  test1(0, 1) = -1;
  test1(1, 0) = 7;
  test1(1, 1) = -2;
  EXPECT_EQ(test1.Determinant(), 1);
}

TEST(Calc, Inverse) {
  S21Matrix test1(3, 3);
  test1(0, 0) = 2;
  test1(0, 1) = 5;
  test1(0, 2) = 7;
  test1(1, 0) = 6;
  test1(1, 1) = 3;
  test1(1, 2) = 4;
  test1(2, 0) = 5;
  test1(2, 1) = -2;
  test1(2, 2) = -3;

  test1 = test1.InverseMatrix();
  EXPECT_EQ(test1(0, 0), 1);
  EXPECT_EQ(test1(0, 1), -1);
  EXPECT_EQ(test1(0, 2), 1);
  EXPECT_EQ(test1(1, 0), -38);
  EXPECT_EQ(test1(1, 1), 41);
  EXPECT_EQ(test1(1, 2), -34);
  EXPECT_EQ(test1(2, 0), 27);
  EXPECT_EQ(test1(2, 1), -29);
  EXPECT_EQ(test1(2, 2), 24);
}

}  // namespace test