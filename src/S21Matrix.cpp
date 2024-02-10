#include "S21Matrix.hpp"

S21Matrix::S21Matrix() : _rows(2), _cols(2) {
  _matrix = new double *[_rows];

  for (int i = 0; i < _rows; i++) {
    _matrix[i] = new double[_cols]();
  }
}

S21Matrix::S21Matrix(int rows, int cols) : _rows(rows), _cols(cols) {
  if (rows < 0 || cols < 0) {
    throw std::length_error("invalid matrix size");
  }

  _matrix = new double *[_rows];

  for (int i = 0; i < _rows; i++) {
    _matrix[i] = new double[_cols]();

    if (_matrix[i] == nullptr) {
      throw std::runtime_error("memory error");
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : S21Matrix(other._rows, other._cols) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = other._matrix[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : _rows(other._rows), _cols(other._cols), _matrix(other._matrix) {
  other._rows = 0;
  other._cols = 0;
  other._matrix = nullptr;
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < _rows; i++) {
    delete[] _matrix[i];
  }

  delete[] _matrix;
}

int S21Matrix::GetRows() const noexcept { return _rows; }

int S21Matrix::GetCols() const noexcept { return _cols; }

void S21Matrix::SetRows(int rows) {
  if (rows == _rows) return;
  if (rows <= 0) {
    throw std::length_error("invalid rows size");
  }

  double **tmp = new double *[rows];
  int minRows = std::min(rows, _rows);

  for (int i = 0; i < minRows; i++) {
    tmp[i] = new double[_cols];

    for (int j = 0; j < _cols; j++) {
      tmp[i][j] = _matrix[i][j];
    }
  }

  for (int i = minRows; i < rows; i++) {
    tmp[i] = new double[_cols]();
  }

  this->~S21Matrix();
  _rows = rows;
  _matrix = tmp;
}

void S21Matrix::SetCols(int cols) {
  if (cols == _cols) return;
  if (cols < 0) {
    throw std::length_error("invalid cols size");
  }

  for (int i = 0; i < _rows; i++) {
    double *tmp = new double[cols]();

    for (int j = 0; i < cols; j++) {
      if (j < _cols) tmp[j] = _matrix[i][j];
    }

    delete[] _matrix[i];
    _matrix[i] = tmp;
  }

  _cols = cols;
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool res = true;

  if (_rows != other._rows || _cols != other._cols) {
    res = false;
    return res;
  }

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      if (fabs(_matrix[i][j] - other._matrix[i][j]) > 1e-7) {
        res = false;
        return res;
      }
    }
  }

  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range("Incorrect input");
  }

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] += other._matrix[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range("Incorrect input");
  }

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] -= other._matrix[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (_cols != other._rows) {
    throw std::out_of_range("Incorrect input");
  }

  S21Matrix res(_rows, other._cols);

  for (int i = 0; i < res._rows; i++) {
    for (int j = 0; j < res._cols; j++) {
      for (int k = 0; k < _cols; k++) {
        res._matrix[i][j] += _matrix[i][k] * other._matrix[k][j];
      }
    }
  }

  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(_cols, _rows);

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      res._matrix[j][i] = _matrix[i][j];
    }
  }

  return res;
}

double S21Matrix::Determinant() {
  if (_rows != _cols) {
    throw std::logic_error("Invalid matrix size");
  }

  if (_rows == 1) {
    return _matrix[0][0];
  }

  if (_rows == 2) {
    return _matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1];
  }

  double result = 1;
  double number = 0;
  S21Matrix temp(*this);

  for (int k = 0; k < _rows; k++) {
    for (int i = k; i < _rows; i++) {
      number = temp(i, k);
      if (fabs(number) >= 1e-6) {
        for (int j = 0; j < _cols; j++) {
          temp(i, j) /= number;
        }
        result *= number;
      } else if (i == k) {
        return 0.0;
      }
    }

    for (int i = k + 1; i < _rows; i++) {
      if (fabs(temp(i, k) >= 1e-6)) {
        for (int j = k; j < _cols; j++) {
          temp(i, j) -= temp(k, j);
        }
      }
    }
  }

  result *= llround(temp(_rows - 1, _cols - 1) * 1e7) / 1e7;
  result = llround(result * 1e7) / 1e7;

  return result;
}

// double S21Matrix::DetTwo() const {
//   return _matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1];
// }

double S21Matrix::GetMinor(int i, int j) {
  double minor = 0.0;

  if (_rows == 2) {
    for (int m = 0; m < _rows; m++) {
      for (int n = 0; n < _cols; n++) {
        if (m != i && n != j) {
          minor = _matrix[m][n];
        }
      }
    }
  } else {
    int x = 0;
    int y = 0;
    S21Matrix res(_rows - 1, _cols - 1);

    for (int m = 0; m < _rows; m++) {
      if (m != i) {
        for (int n = 0; n < _cols; n++) {
          if (n != j) {
            res(x, y) = _matrix[m][n];
            y++;
          }
        }
        x++;
        y = 0;
      }
    }

    minor = res.Determinant();

    // if (res._cols > 2) {
    //   minor = res.Determinant();
    // } else {
    //   minor = res.DetTwo();
    // }
  }

  return minor;
}

void S21Matrix::ResMatrix(S21Matrix &result) {
  for (int i = 0; i < result._rows; i++) {
    for (int j = 0; j < result._cols; j++) {
      result(i, j) = GetMinor(i, j);
    }
  }
}

void S21Matrix::Complements(S21Matrix &result) {
  ResMatrix(result);

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      result(i, j) = result(i, j) * pow(-1, i + j);
    }
  }
}

S21Matrix S21Matrix::CalcComplements() {
  if (_rows != _cols) {
    throw std::logic_error("Invalid matrix size");
  }

  S21Matrix result(*this);

  if (_rows == 1) {
    result(0, 0) = 1;
  } else {
    Complements(result);
  }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result(_rows, _cols);
  double det = Determinant();

  if (fabs(det) < 1e-6) {
    throw std::logic_error("determinant of matrix = 0");
  } else {
    if (_rows == 1) {
      result(0, 0) = 1 / _matrix[0][0];
    } else {
      S21Matrix complement(_rows, _cols);
      complement = CalcComplements();
      result = complement.Transpose();
      result.MulNumber(1 / det);
    }
  }

  return result;
}

double &S21Matrix::operator()(int rows, int cols) {
  if (rows < 0 || cols < 0 || rows >= _rows || cols >= _cols) {
    throw std::out_of_range("Invalid element");
  }

  return _matrix[rows][cols];
}

double &S21Matrix::operator()(int rows, int cols) const {
  if (rows < 0 || cols < 0 || rows >= _rows || cols >= _cols) {
    throw std::out_of_range("Invalid element");
  }

  return _matrix[rows][cols];
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix(other._rows, other._cols);

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = other._matrix[i][j];
    }
  }

  return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
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

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

//////
S21Matrix operator*(const double num, const S21Matrix &matrix) {
  S21Matrix res(matrix);
  res.MulNumber(num);
  return res;
}