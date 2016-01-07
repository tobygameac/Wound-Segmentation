#pragma once

#include <cmath>

namespace WoundSegmentation {

  namespace ImageProcessing {

    template<class T>
    class RGBA {

    public:

      RGBA() : r_(0), g_(0), b_(0), a_(1) {
      }

      RGBA(T r, T g, T b) : r_(r), g_(g), b_(b), a_(1) {
      }

      RGBA(T r, T g, T b, T a) : RGBA(r, g, b), a_(a) {
      }

      RGBA<T> Square() const {
        return RGBA<T>(r_ * r_, g_ * g_, b_ * b_, a_);
      }

      RGBA<T> Sqrt() const {
        return RGBA<T>(sqrt(r_), sqrt(g_), sqrt(b_), a_);
      }

      T Magnitude() const {
        return r_ * r_ + g_ * g_ + b_ * b_;
      }

      T NormalizedMagnitude() const {
        return sqrt((r_ * r_ + g_ * g_ + b_ * b_) / 3.0);
      }

      const RGBA<T> &operator =(const RGBA<T> &other) {
        r_ = other.r_;
        g_ = other.g_;
        b_ = other.b_;
        a_ = other.a_;

        return *this;
      }

      bool operator ==(const RGBA<T> &other) const {
        if (r_ != other.r_) {
          return false;
        }

        if (g_ != other.g_) {
          return false;
        }

        if (b_ != other.b_) {
          return false;
        }

        if (a_ != other.a_) {
          return false;
        }

        return true;
      }

      T r_, g_, b_, a_;
    };
  }
}