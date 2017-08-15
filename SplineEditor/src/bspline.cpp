#include "bspline.h"
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

using std::vector;
using alglib::real_1d_array;
using alglib::real_2d_array;

typedef Eigen::Matrix<float, 6, 3> Matrix6x3;
typedef Eigen::VectorXf VectorX;

ts::BSpline constructBSpline(Polyline& controlPoints, size_t degree) {
	ts::BSpline spline = ts::BSpline(degree, 3, controlPoints.size(), TS_CLAMPED);
	auto splinePoints = spline.ctrlp();
	for (size_t i = 0; i < controlPoints.size(); ++i) {
		splinePoints[3 * i] = controlPoints[i](0);
		splinePoints[3 * i + 1] = controlPoints[i](1);
		splinePoints[3 * i + 2] = controlPoints[i](2);
	}
	spline.setCtrlp(splinePoints);
	return spline;
}

Point evaluate(const ts::BSpline& spline, float val) {
	auto result = spline.evaluate(val).result();
	return Point(result[0], result[1], result[2]);
}

//works only for arc-length parameterization
float kappa(ts::BSpline& d1, ts::BSpline& d2, float x) {
	return evaluate(d1, x).cross(evaluate(d2, x)).norm() / pow(evaluate(d1, x).norm(), 3);
}

//curvature from local tangents
float kappa2(ts::BSpline& spline, ts::BSpline& d1, float x) {
	//FIXME use better delta
	float delta = x < (spline.knots().front() + spline.knots().back()) / 2 ? 0.01 : -0.01;

	Point v1 = evaluate(d1, x).normalized();
	Point v2 = evaluate(d1, x + delta).normalized();
	return v1.dot(v2) / (evaluate(spline, x) - evaluate(spline, x + delta)).norm();
}

//calculate from curvature radius
float kappa3(ts::BSpline& spline, ts::BSpline& d1, float x) {
	//FIXME use better delta
	float delta = 0.01;

	Point v1, v2, v3;
	v2 = evaluate(spline, x);

	if (x - delta >= spline.knots().front()) {
		v1 = evaluate(spline, x - delta);
	} else {
		v1 = v2 - delta * evaluate(d1, x).normalized();
	}
	if (x + delta <= spline.knots().back()) {
		v3 = evaluate(spline, x + delta);
	} else {
		v3 = v2 + delta * evaluate(d1, x).normalized();
	}

	//construct the bisectors to find the circle's center
	Point planeNormal = (v1 - v2).cross(v2 - v3);
	Point a1 = (v1 + v2) / 2;
	Point b1 = planeNormal.cross(v1 - v2);
	Point a2 = (v2 + v3) / 2;
	Point b2 = planeNormal.cross(v2 - v3);

	if (abs(b1.dot(b2) / (b1.norm() * b2.norm())) == 1) {
		return 0;
	}

	Matrix6x3 M = Matrix6x3::Zero();
	VectorX values(6);
	for (size_t i = 0; i < 3; ++i) {
		size_t i1 = (i + 1) % 3, i2 = (i + 2) % 3;

		M(i, i1) = -b1[i2];
		M(i, i2) = b1[i1];
		values[i] = b1[i1] * a1[i2] - b1[i2] * a1[i1];

		M(i + 3, i1) = -b2[i2];
		M(i + 3, i2) = b2[i1];
		values[i + 3] = b2[i1] * a2[i2] - b2[i2] * a2[i1];
	}
	//VectorX result = M.fullPivLu().solve(values);
	Eigen::FullPivLU<Matrix6x3> lu(M);
	VectorX result = lu.solve(values);

	return 1 / (Point(result[0], result[1], result[2]) - v2).norm();
}

float curvature(ts::BSpline& spline) {
	ts::BSpline d1 = spline.derive();
	ts::BSpline d2 = d1.derive();

	GausLegendreIntegrator<float> integral(spline.deg());
	//auto fcurv = std::bind(kappa, d1, d2, std::placeholders::_1);
	//auto fcurv = std::bind(kappa2, spline, d1, std::placeholders::_1);
	auto fcurv = std::bind(kappa3, spline, d1, std::placeholders::_1);
	size_t num_samples = spline.order() * spline.nCtrlp();

	float kMin = spline.knots().front();
	float kMax = spline.knots().back();

	float result = 0;

	for (size_t i = 0; i < num_samples; ++i) {
		result += integral(fcurv, kMin + (i / (double) num_samples) * (kMax - kMin),
				kMin + ((i + 1) / (double) num_samples) * (kMax - kMin));
	}

	return result;
}

//---------------------------------------------------------
//------------------ ArcLengthBSpline ---------------------
//---------------------------------------------------------

Point ArcLengthBSpline::getPointAtNormalizedLength(double param) {
	return getPointAtLength(param * length);
}

Point ArcLengthBSpline::getDerivativeAtNormalizedLength(double param) {
	return getDerivativeAtLength(param * length);
}

Point ArcLengthBSpline::getPointAtLength(double evalLength) {
	return evaluate(generalSpline, lengthToGeneralParam(evalLength));
}

double firstDerivative(const alglib::spline1dinterpolant& s, const double x) {
	double ds0, ds1, ds2;
	alglib::spline1ddiff(s, x, ds0, ds1, ds2);
	return ds1;
}

Point ArcLengthBSpline::getDerivativeAtLength(double evalLength) {
	return evaluate(generalSplineD1, lengthToGeneralParam(evalLength))
			* firstDerivative(inverseLength, clamp(evalLength, 0.0, length));
}

double ArcLengthBSpline::lengthToGeneralParam(double evalLength) {
	return alglib::spline1dcalc(inverseLength, clamp(evalLength, 0.0, length));
}

//functor to give the length of the derivative of the spline
class SplineLengthFunction {
public:
	SplineLengthFunction(BSpline& spline) :
			spline(spline.derive()) {
	}

	double operator()(double param) {
		return evaluate(spline, param).norm();
	}

private:
	BSpline spline;
};

//based on John W. Peterson: Arc Length Parameterization of Spline Curves
void ArcLengthBSpline::generateMapping() {
	float kMin = generalSpline.knots().front();
	float kMax = generalSpline.knots().back();

	int samples = generalSpline.nCtrlp() * generalSpline.order();

	real_1d_array lengths, knots;
	lengths.setlength(samples + 1);
	knots.setlength(samples + 1);
	lengths[0] = 0;
	knots[0] = 0;

	SplineLengthFunction arcLength(generalSpline);
	GausLegendreIntegrator<double> glIntegrate(degree);
	for (double i = 0; i < samples; ++i) {
		double a = kMin + (i / samples) * (kMax - kMin);
		double b = kMin + ((i + 1) / samples) * (kMax - kMin);
		lengths[i + 1] = lengths[i] + glIntegrate(arcLength, a, b);
		knots[i + 1] = b;
	}
	length = lengths[samples];

//or better? alglib::spline1dbuildakima(lengths, knots, inverseLength);
	alglib::spline1dbuildcubic(lengths, knots, samples + 1, 2, 0, 2, 0, inverseLength); //natural bound conditions
}
