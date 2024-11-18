#pragma once

#include <new>
#include <utility>

template <typename T>
class FixedArray {
 public:
  FixedArray() : _size(0), _data(nullptr) {}
  /// size = 0 if allocation fails
  FixedArray(size_t fixedSize) { Allocate(fixedSize); }

  /// @return -1 if allocation fails.
  int Allocate(size_t size) {
    if (_data != nullptr) {
      delete[] _data;
      _data = nullptr;
    }
    _size = 0;
    try {
      _data = new T[size];
    } catch (std::bad_alloc& ba) {
      return -1;
    }
    _size = size;
    return 0;
  }

  size_t size() const { return _size; }
  ~FixedArray() {
    if (_data != nullptr) {
      delete[] _data;
    }
    _data = nullptr;
    _size = 0;
  }
  T& operator[](size_t i) { return _data[i]; }
  const T& operator[](size_t i) const { return _data[i]; }

  FixedArray(const FixedArray& other) : _size(other._size) {
    if (_size > 0) {
      _data = (_size ? new T[_size] : nullptr);
      std::copy(other._data, other._data + _size, _data);
    }
  }

  FixedArray(FixedArray&& other) noexcept : FixedArray() { swap(*this, other); }

  friend void swap(FixedArray& first, FixedArray& second) {
    using std::swap;
    swap(first._size, second._size);
    swap(first._data, second._data);
  }

  FixedArray& operator=(const FixedArray& other) {
    FixedArray temp(other);
    swap(*this, temp);

    return *this;
  }

  const T& back() const { return _data[_size - 1]; }
  T& back() { return _data[_size - 1]; }
  bool empty() const { return _data == nullptr; }

  using const_iterator = T const*;
  using iterator = T*;

  const_iterator cbegin() const { return _data; };
  const_iterator cend() const { return _data + _size; };
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }
  iterator begin() { return _data; };
  iterator end() { return _data + _size; };

 private:
  T* _data = nullptr;
  size_t _size;
};