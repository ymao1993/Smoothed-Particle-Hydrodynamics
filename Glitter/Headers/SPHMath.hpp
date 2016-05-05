/**
 * Math utilities used by SPHSimulator
 *
 * modified from CMU 15462's utility class:
 * https://github.com/462cmu/CMU462
 *
 */

#ifndef SPHMath_h
#define SPHMath_h

#include <ostream>
#include <cmath>

namespace SPHSim
{
	const float PI = 3.1415926;
	const float EPSILON = 0.001;

	struct vec3 {
	  // components
	  double x, y, z;

	  /**
	   * Constructor.
	   * Initializes tp vector (0,0,0).
	   */
	  vec3() : x(0.0), y(0.0), z(0.0) {}

	  /**
	   * Constructor.
	   * Initializes to vector (x,y,z).
	   */
	  vec3(double x, double y, double z) : x(x), y(y), z(z) {}

	  /**
	   * Constructor.
	   * Initializes from existing vector
	   */
	  vec3(const vec3 &v) : x(v.x), y(v.y), z(v.z) {}

	  // returns reference to the specified component (0-based indexing: x, y, z)
	  inline double &operator[](const int &index) { return (&x)[index]; }

	  // returns const reference to the specified component (0-based indexing: x, y,
	  // z)
	  inline const double &operator[](const int &index) const {
	    return (&x)[index];
	  }

	  inline bool operator==(const vec3 &v) const {
	    return v.x == x && v.y == y && v.z == z;
	  }

	  // negation
	  inline vec3 operator-(void) const { return vec3(-x, -y, -z); }

	  // addition
	  inline vec3 operator+(const vec3 &v) const {
	    return vec3(x + v.x, y + v.y, z + v.z);
	  }

	  // subtraction
	  inline vec3 operator-(const vec3 &v) const {
	    return vec3(x - v.x, y - v.y, z - v.z);
	  }

	  // right scalar multiplication
	  inline vec3 operator*(const double &c) const {
	    return vec3(x * c, y * c, z * c);
	  }

	  // scalar division
	  inline vec3 operator/(const double &c) const {
	    const double rc = 1.0 / c;
	    return vec3(rc * x, rc * y, rc * z);
	  }

	  // addition / assignment
	  inline void operator+=(const vec3 &v) {
	    x += v.x;
	    y += v.y;
	    z += v.z;
	  }

	  // subtraction / assignment
	  inline void operator-=(const vec3 &v) {
	    x -= v.x;
	    y -= v.y;
	    z -= v.z;
	  }

	  // scalar multiplication / assignment
	  inline void operator*=(const double &c) {
	    x *= c;
	    y *= c;
	    z *= c;
	  }

	  // scalar division / assignment
	  inline void operator/=(const double &c) { (*this) *= (1. / c); }

	  /**
	   * Returns Euclidean length.
	   */
	  inline double norm(void) const { return sqrt(x * x + y * y + z * z); }

	  /**
	   * Returns Euclidean length squared.
	   */
	  inline double norm2(void) const { return x * x + y * y + z * z; }

	  /**
	   * Returns unit vector.
	   */
	  inline vec3 unit(void) const {
	    double rNorm = 1. / sqrt(x * x + y * y + z * z);
	    return vec3(rNorm * x, rNorm * y, rNorm * z);
	  }

	  /**
	   * Divides by Euclidean length.
	   */
	  inline void normalize(void) { (*this) /= norm(); }

	}; // class vec3

	// left scalar multiplication
	inline vec3 operator*(const double &c, const vec3 &v) {
	  return vec3(c * v.x, c * v.y, c * v.z);
	}

	// dot product (a.k.a. inner or scalar product)
	inline double dot(const vec3 &u, const vec3 &v) {
	  return u.x * v.x + u.y * v.y + u.z * v.z;
	}

	// cross product
	inline vec3 cross(const vec3 &u, const vec3 &v) {
	  return vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
	                  u.x * v.y - u.y * v.x);
	}

	// prints components
	inline std::ostream &operator<<(std::ostream &os, const vec3 &v)
	{
	  os << "(" << v.x << "," << v.y << "," << v.z << ")";
	  return os;
	}

	/**
	 * clamp p if out of bounds
	 */
	inline void box_clamp(vec3 &p, const vec3 &min, const vec3& max)
	{
		p.x = p.x < min.x ? min.x : p.x > max.x ? max.x : p.x;
		p.y = p.y < min.y ? min.y : p.y > max.y ? max.y : p.y;
		p.z = p.z < min.z ? min.z : p.z > max.z ? max.z : p.z;
		
		return;
	}

	/**
	 * clamp p and reflect q if out of bounds.
	 */
	inline void box_clamp_and_reflect(vec3 &p, vec3 &q, const vec3 &min, const vec3& max, float damping = 0)
	{
		if(p.x < min.x) { p.x = min.x; q.x *= -(1.f-damping); }
		if(p.x > max.x) { p.x = max.x; q.x *= -(1.f-damping); }
		if(p.y < min.y) { p.y = min.y; q.y *= -(1.f-damping); }
		if(p.y > max.y) { p.y = max.y; q.y *= -(1.f-damping); }
		if(p.z < min.z) { p.z = min.z; q.z *= -(1.f-damping); }
		if(p.z > max.z) { p.z = max.z; q.z *= -(1.f-damping); }

		return;
	}

}

#endif