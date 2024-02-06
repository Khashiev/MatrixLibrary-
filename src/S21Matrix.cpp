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
