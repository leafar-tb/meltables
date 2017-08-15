#pragma once

#include <iostream>
#include <functional>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <libalglib/integration.h>

#include "igl/viewer/Viewer.h"

#define RAD_TO_DEG (180.0 / M_PI)
#define DEG_TO_RAD (M_PI / 180.0)

typedef Eigen::Vector3d Point;
typedef Eigen::Vector2d Point2D;
typedef std::vector<Point> Polyline;

double length(const Polyline* pl);
double length(const Polyline& pl);

typedef struct Mesh_s {
	Mesh_s(size_t vCount = 0, size_t fCount = 0) :
			vertices(vCount, 3), faces(fCount, 3) {
	}
	Mesh_s(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces) :
			vertices(vertices), faces(faces) {
	}

	Eigen::MatrixXd vertices;
	Eigen::MatrixXi faces;
} Mesh;

Mesh& translate(Mesh& mesh, const Point& dP);
Mesh& rotate(Mesh& mesh, const Eigen::Matrix3d& rot);
void saveMesh(std::string name, const Mesh& mesh);
//a new mesh with all vertices and faces from m1 and m2 combined
Mesh merge(const Mesh& m1, const Mesh& m2);

//allows integration of functions using the gauss-legendre method
//this has the advantage, that the points and weights can be stored and reused
//
//Value needs to support multiplication with double and summation
template<typename Value>
class GausLegendreIntegrator {
public:
	typedef std::function<Value(double)> FunctionType;

	GausLegendreIntegrator(size_t degree) {
		alglib::ae_int_t info;
		alglib::gqgenerategausslegendre(degree, info, points, weights);
		if (info < 0) {
			std::cerr << "alglib::gqgenerategausslegendre() failed! Error code: " << info;
		}
	}

	//integrate f over [lower, upper]
	//you can pass an explicit zero value, in case it is not provided through the default constructor
	Value integrate(FunctionType f, double lower, double upper, Value zero) {
		Value integral = zero;
		for (int j = 0; j < points.length(); ++j) {
			//[-1,1] -> [lower, upper]
			double param = (upper - lower) * (points[j] / 2) + (lower + upper) / 2;
			//factor also from interval change
			integral += (upper - lower) * 0.5 * weights[j] * f(param);
		}
		return integral;
	}

	//shortcut, when the default constructor of Value gives the zero
	Value integrate(FunctionType f, double lower, double upper) {
		return integrate(f, lower, upper, defaultVal);
	}

	Value operator()(FunctionType f, double lower, double upper) {
		return integrate(f, lower, upper, defaultVal);
	}

	Value operator()(FunctionType f, double lower, double upper, Value zero) {
		return integrate(f, lower, upper, zero);
	}

private:
	alglib::real_1d_array points;
	alglib::real_1d_array weights;
	Value defaultVal;
};

template<typename val>
val clamp(val value, val min, val max) {
	if (value < min)
		return min;
	if (value > max)
		return max;
	return value;
}

Point worldToScreen(igl::viewer::Viewer *viewer, Point worldPoint);
Point2D worldToScreenFlat(igl::viewer::Viewer *viewer, Point worldPoint);

//angle in DEG
double angle(Point vec1, Point vec2);

// returns a normalized vector, orthogonal to the input, if possible
// returns the Zero vector otherwise
Point getOrthogonal(Point& vec);

// a rotation object that corresponds to the rotation of vector from onto vector to
Eigen::Matrix3d getRotation(Point from, Point to);

//b1, b2 should be orthogonal and normalized
//angle is in RAD
double getAngleInBasis(const Point& v, const Point& b1, const Point& b2);
