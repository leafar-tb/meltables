#pragma once

#include <vector>

#include <libalglib/interpolation.h>
#include "tinyspline/tinysplinecpp.h"
#include "util.h"

using ts::BSpline;

BSpline constructBSpline(Polyline& controlPoints, size_t degree);

Point evaluate(const BSpline& spline, float val);

float curvature(BSpline& spline);

//a BSpline with arc length parameterization
class ArcLengthBSpline {
public:
	ArcLengthBSpline(Polyline controlPoints, size_t degree) {
		this->controlPoints = controlPoints;
		this->degree = degree;
		generalSpline = constructBSpline(this->controlPoints, degree);
		generalSplineD1 = generalSpline.derive();
		length = -1;

		generateMapping();
	}

	//clamped to [0,length]
	Point getPointAtLength(double length);
	Point getDerivativeAtLength(double length);

	//normalized to [0,1]
	Point getPointAtNormalizedLength(double param);
	Point getDerivativeAtNormalizedLength(double param);

	double getLength() {
		return length;
	}

	size_t getDegree() {
		return degree;
	}

private:
	Polyline controlPoints;
	size_t degree;
	double length;

	//the underlying spline curve without special parameterization
	BSpline generalSpline;
	BSpline generalSplineD1;

	//spline approx of the inverse length function
	alglib::spline1dinterpolant inverseLength;

	//get the parameter in [0,1] for the underlying spline, that corresponds to length
	double lengthToGeneralParam(double evalLength);

	//calculate the length and a mapping from length to general parameter
	void generateMapping();

};
