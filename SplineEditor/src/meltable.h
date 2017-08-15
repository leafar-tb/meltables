#pragma once

#include "bspline.h"
#include "util.h"
#include <vector>
#include <utility>

struct MeltableConstraints {
	size_t numSegments = 10; //desired number of segments
	double minTheta = 20; //in degrees
	double maxTheta = 90; //in degrees

	//parameters for the exported model
	double targetLength = 10;
	double radius = 0.25;
	double jointThickness = 0.1;

	size_t samplingLOD = 100; //number of samples for meltable calculation
	size_t cylinderLOD = 32; //number of points to approx the circles
};

//see paper 'Meltables: Fabrication of Complex 3D Curves by Melting'
class Meltable {
public:
	Meltable(ArcLengthBSpline& spline, MeltableConstraints& constraints);

	Polyline getPoints() {
		return q;
	}

	double length() {
		return _length;
	}

	size_t numSegments() {
		return q.size() - 1;
	}

	Point segment(size_t seg) {
		if (seg >= numSegments())
			return Point::Zero();
		return q[seg + 1] - q[seg];
	}

	//see calculateBases()
	std::pair<Point, Point> segmentBase(size_t seg) {
		if (seg >= numSegments())
			return std::make_pair(Point::Zero(), Point::Zero());
		return segmentBases[seg];
	}

	Mesh generateMesh();

	Mesh generatePreview();

	//-----
	// some static functions we want to reuse
	//-----

	//contains for each segment a pair of normalized vectors that are orthogonal to it and each other
	//these represent the direction of the x- and y- axis if the meltable were straightened along the z-axis
	static std::vector<std::pair<Point, Point> > calculateBases(Polyline& q);

	//here you can define which vectors to use as base for the first segment
	static std::vector<std::pair<Point, Point> > calculateBases(Polyline& q, const std::pair<Point, Point>& firstBase);

	static Mesh generateMesh(const Polyline& q, const std::vector<std::pair<Point, Point> >& segmentBases,
			const MeltableConstraints& constraints, bool preview = false);

private:
	MeltableConstraints constraints;

	//the points q_i of the meltable
	Polyline q;

	double _length;

	std::vector<std::pair<Point, Point> > segmentBases;

	//private implementation class
	class Calculation;

};

