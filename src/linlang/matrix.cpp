#include "matrix.hpp"
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace math {
namespace linlang {

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, T(0)) {
  if (rows == 0 || cols == 0)
    throw std::invalid_argument("Matrix dimensions must be positive");
}

template <typename T>
Matrix<T>::Matrix(const Matrix &other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init) {
  rows_ = init.size();
  if (rows_ == 0)
    throw std::invalid_argument("Empty initializer list");

  cols_ = init.begin()->size();
  data_.reserve(rows_ * cols_);
  for (const auto &row : init) {
    if (row.size() != cols_)
      throw std::invalid_argument("Irregular matrix initialization");
    data_.insert(data_.end(), row.begin(), row.end());
  }
}

template <typename T> Matrix<T> Matrix<T>::Indentify(size_t size) {
  Matrix<T> result(size, size);
  for (size_t i = 0; i < size; i++)
    result(i, i) = T(1);
  return result;
}

template <typename T> Matrix<T> Matrix<T>::Rotation2D(T angle) {
  Matrix<T> R(2, 2);
  T cos_theta = std::cos(angle);
  T sin_theta = std::sin(angle);
  R(0, 0) = cos_theta;
  R(0, 1) = -sin_theta;
  R(1, 0) = sin_theta;
  R(1, 1) = -cos_theta;

  return R;
}

template <typename T>
Matrix<T> Matrix<T>::Rotation3D(T angle, const std::vector<T> &axis) {
  Matrix<T> R(3, 3);
  T cos_theta = std::cos(angle);
  T sin_theta = std::sin(angle);
  T one_minus_cos = T(1) - cos_theta;

  std::vector<T> a = axis.normalized();
  R(0, 0) = cos_theta + a.x() * a.x() * one_minus_cos;
  R(0, 1) = a.x() * a.y() * one_minus_cos - a.z() * sin_theta;
  R(0, 2) = a.x() * a.z() * one_minus_cos + a.y() * sin_theta;

  R(1, 0) = a.y() * a.x() * one_minus_cos + a.z() * sin_theta;
  R(1, 1) = cos_theta + a.y() * a.y() * one_minus_cos;
  R(1, 2) = a.y() * a.z() * one_minus_cos - a.x() * sin_theta;

  R(2, 0) = a.z() * a.x() * one_minus_cos - a.y() * sin_theta;
  R(2, 1) = a.z() * a.y() * one_minus_cos + a.x() * sin_theta;
  R(2, 2) = cos_theta + a.z() * a.z() * one_minus_cos;

  return R;
}

} // namespace linlang

} // namespace math
