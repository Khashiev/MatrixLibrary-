#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <utility>

class S21Matrix {
 private:
  int _rows, _cols;
  double** _matrix;

  // helpers
  void Complements(S21Matrix& result);
  void ResMatrix(S21Matrix& result);
  double GetMinor(int i, int j);
  // double DetTwo() const;

 public:
  // constructors
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  // accessors and mutators
  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows);
  void SetCols(int cols);

  // methods
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  // operators
  double& operator()(int rows, int cols);
  double& operator()(int rows, int cols) const;
  S21Matrix& operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(const double num, const S21Matrix& matrix);
};

#endif
