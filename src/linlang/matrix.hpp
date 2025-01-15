#pragma once

#include <cstddef>
#include <initializer_list>
#include <string>
#include <tuple>
#include <vector>

namespace math {
namespace linlang {
template <typename T> class Matrix {

public:
  Matrix(size_t rows, size_t cols);
  Matrix(const Matrix &other);
  Matrix(Matrix &&other) noexcept;
  Matrix(std::initializer_list<std::initializer_list<T>> init);

  // Factory Methods
  static Matrix Indentify(size_t size);
  static Matrix Zeros(size_t rows, size_t cols);
  static Matrix Ones(size_t rows, size_t cols);
  static Matrix Random(size_t rows, size_t cols, T min = 0, T max = 1);
  static Matrix Rotation2D(T angle);
  static Matrix Rotation3D(T angle, const std::vector<T> &axis);
  static Matrix Translation3D(const std::vector<T> &translation);

  // Basic Operations

  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &other) noexcept;
  T &operator()(size_t row, size_t col);
  const T &operator()(size_t row, size_t col) const;

  // Arithmetics
  Matrix operator+(const Matrix &other) const;
  Matrix operator-(const Matrix &other) const;
  Matrix operator*(const Matrix &other) const;
  Matrix operator*(T scalar) const;

  Matrix transpose() const;
  T determinant() const;
  Matrix inverse() const;
  std::tuple<Matrix, Matrix, Matrix> svd() const;
  std::tuple<Matrix, Matrix> eigenDecomposition() const;
  T trace() const;
  size_t rank() const;
  Matrix psuedoInverse() const;

  // Decompositions
  std::tuple<Matrix, Matrix> lu() const;
  std::tuple<Matrix, Matrix> qr() const;
  Matrix cholesky() const;

  std::vector<T> solve(const std::vector<T> &b) const;

  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  bool isSquare() const { return rows_ == cols_; }
  bool isSymmetric() const;
  bool isOrthogonal() const;
  bool isPositiveDefinite() const;
  T norm(const std::string &type = "frobenius") const;
  T condition() const;

private:
  size_t rows_;
  size_t cols_;
  std::vector<T> data_;

  void householderReduction();
  void checkDimension(const Matrix &other) const;
  size_t index(size_t row, size_t col) const;
};

using Matrixf = Matrix<float>;
using Matrixd = Matrix<double>;

} // namespace linlang

} // namespace math
