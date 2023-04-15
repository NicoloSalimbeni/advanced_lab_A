// Copyright 2023 nicolò salimbeni andrea de vita
#ifndef AnUtil_h
#define AnUtil_h
/*
  This class contains general functions and object than can be usefull
  in any data analysis, they are not related with this specific
  experiment.
*/

#include <RtypesCore.h>

#include <cmath>
#include <string>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TVectorD.h"

class AnUtil {
 public:
  AnUtil();
  ~AnUtil();

  // template functions must be declared and defined in the same file
  // put all of them below here
  template <class T>
  static void Reorder(T* v, int n, const std::string order = "ascending") {
    // order a vector, array, whatever with ascending or descending order
    T tmp;
    if (order == "ascending") {
      for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
          if (v[j] < v[i]) {
            tmp  = v[i];
            v[i] = v[j];
            v[j] = tmp;
          }
        }
      }
    }

    if (order == "descending") {
      for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
          if (v[j] > v[i]) {
            tmp  = v[i];
            v[i] = v[j];
            v[j] = tmp;
          }
        }
      }
    }
  }

  template <class T>
  static Double_t Mean(const T& v, int n) {
    // compute the mean value of a list v
    Double_t sum = 0;
    for (int i = 0; i < n; i++) {
      sum += v[i];
    }
    return sum / n;
  }

  template <class T>
  static Double_t Rms(const T& v, int n) {
    // compute the RMS value of a list v
    Double_t mean = Mean(v, n);
    Double_t sum  = 0;
    for (int i = 0; i < n; i++) {
      sum += std::pow(v[i] - mean, 2);
    }

    // pay attenction to the n-1 at the denominator
    return std::sqrt(sum / (n - 1));
  }
};

#endif