#ifndef EVAL_H
#define EVAL_H

#include <stdint.h>
#include <utility>

// typedef pair<double, long long> my_pair;

class Eval
{
private:
  uint64_t n;
  double sum;
  // double sum_sqr;

public:
  Eval();  // constructor
  ~Eval(); // destructor

  /// @brief Returns the number of samples 
  uint64_t get_n() { return n; };
  /// @brief Returns the sum of the samples
  double get_sum() { return sum; };
  /// @brief Sets the number of samples
  void set_n(uint64_t _n) { n = _n; };
  /// @brief Sets the sum of the samples
  void set_sum(double s) { sum = s; };

  /// @brief Resets the Eval to n = 0 and sum = 0
  void reset();
  /// @brief Returns the mean of the samples
  double get_mean();

  // double get_variance();

  /// @brief Adds a sample to the Eval, raising n by 1
  Eval &operator<<(const double d);
  
  // Eval& operator<<(const my_pair p);
};

#endif
