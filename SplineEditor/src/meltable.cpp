#include "meltable.h"

#include <cmath>
#include <limits>
#include <map>
#include <iostream>

using namespace std;
using namespace alglib;

//--------------------------------------------------
//             actual implementation
//--------------------------------------------------

//see section 3 'Method' in paper
class Meltable::Calculation {
public:
	Calculation(Meltable* m, ArcLengthBSpline& spline) :
			meltable(m), spline(spline), samples(m->constraints.samplingLOD), ds(spline.getLength() / (samples - 1)) {
		for (size_t i = 0; i < samples; ++i) {
			points.push_back(spline.getPointAtLength(length(i)));
		}
	}

	Polyline run() {
		size_t bestI;
		double bestErr = numeric_limits<double>::infinity();

		for (size_t joints = 0; joints < meltable->constraints.numSegments; ++joints) {
			auto result = getOrCalcF(joints, samples - 1);
			cout << joints << ": " << result.first << endl;
			if (result.first < bestErr) {
				bestI = joints;
				bestErr = result.first;
			}
		}

		cout << endl << bestI << ": " << bestErr << endl << endl;
		return getQVec(bestI, samples - 1);
	}

private:
	size_t samples;
	double ds;
	Meltable* meltable;
	ArcLengthBSpline spline;
	vector<Point> points;

	typedef pair<size_t, size_t> key_type;
	typedef pair<double, size_t> entry_type;

	struct CacheCompare {
		bool operator()(key_type v1, key_type v2) {
			if (v1.first != v2.first) {
				return v1.first < v2.first;
			}
			return v1.second < v2.second;
		}
	};

	//(i,s) -> (f_i(s), s_i)
	map<key_type, entry_type, CacheCompare> cache;

	//ensure that the cache contains the entry and return it
	entry_type getOrCalcF(size_t i, size_t s) {
		if (!hasKey(i, s)) {
			cache.insert(make_pair(make_pair(i, s), f(i, s)));
		}
		return cache.find(make_pair(i, s))->second;
	}

	//approximate the spline with i joints up to length s
	//returns the error and the position of the next joint
	entry_type f(size_t i, size_t s) {
		if (i == 0) {
			return make_pair(h((size_t) 0, (size_t) 0, s), 0);
		}

		size_t si = s - 1;
		entry_type best = make_pair(numeric_limits<double>::infinity(), si);

		while (si > i) {
			auto tmp = getOrCalcF(i - 1, si);
			if (checkAngleConstraint(tmp.second, si, s) && checkJointWidth(tmp.second, si, s)) {
				double err = tmp.first + h(i, si, s);
				if (err < best.first) {
					best = make_pair(err, si);
				}
			}
			si--;
		}
		return best;
	}

	//the approximation-error over the intervall [s,t], if [0,s] is approximated with i joints
	double h(size_t i, size_t s, size_t t) {
		Point start = i == 0 ? points[0] : getQ(i - 1, s);
		Point tang = tangent(s, t);

		double integral = 0;
		for (size_t u = s; u < t; ++u) {
			integral += (points[u] - (start + (length(u) - length(s)) * tang)).squaredNorm() * ds;
		}
		return integral;
	}

	//are the segments wide enough for the joint's angle?
	bool checkJointWidth(size_t s1, size_t s2, size_t s3) {
		double a = angle(tangent(s1, s2), tangent(s2, s3));
		double width = (meltable->constraints.radius * 2 - meltable->constraints.jointThickness)
				* tan(DEG_TO_RAD * a / 2);
		return (targetLength(s2) - targetLength(s1)) / 2 > width && (targetLength(s3) - targetLength(s2)) / 2 > width;
	}

	//is the angle between segments [s1,s2] and [s2,s3] within the given constraints?
	bool checkAngleConstraint(size_t s1, size_t s2, size_t s3) {
		double a = angle(tangent(s1, s2), tangent(s2, s3));
		return meltable->constraints.minTheta <= a && meltable->constraints.maxTheta >= a;
	}

	//get the points q_i of the meltable that approximates the spline with i joints over [0,s]
	Polyline getQVec(size_t i, size_t s) {
		if (i == 0) {
			Polyline q;
			q.push_back(points[0]);
			q.push_back(points[0] + length(s) * tangent(0, s));
			return q;
		}

		auto tmp = getOrCalcF(i, s);
		auto q = getQVec(i - 1, tmp.second);
		q.push_back(q.back() + (length(s) - length(tmp.second)) * tangent(tmp.second, s));
		return q;
	}

	//the last point of the meltable that approximates the spline with i joints over [0,s]
	Point getQ(size_t i, size_t s) {
		return getQVec(i, s).back();
	}

	//is the key present in the cache?
	bool hasKey(size_t i, size_t s) {
		return cache.find(make_pair(i, s)) != cache.end();
	}

	//tangent of the spline at the midpoint of the segment [s1,s2]
	Point tangent(size_t s1, size_t s2) {
		return spline.getDerivativeAtLength((length(s1) + length(s2)) / 2).normalized();
	}

	//length corresponding to the sampling index
	double length(size_t index) {
		return (spline.getLength() * index) / (samples - 1);
	}

	//length corresponding to the sampling index - scaled to targetLength
	double targetLength(size_t index) {
		return (meltable->constraints.targetLength * index) / (samples - 1);
	}

};

vector<pair<Point, Point> > Meltable::calculateBases(Polyline& q) {
	vector < pair<Point, Point> > bases;

	Point segment0 = q[1] - q[0];
	Point base1 = getOrthogonal(segment0);
	Point base2 = segment0.cross(base1).normalized();
	return calculateBases(q, make_pair(base1, base2));
}

vector<pair<Point, Point> > Meltable::calculateBases(Polyline& q, const pair<Point, Point>& firstBase) {
	vector < pair<Point, Point> > bases;

	Point base1 = firstBase.first;
	Point base2 = firstBase.second;
	bases.push_back(firstBase);

	for (size_t i = 0; i < q.size() - 2; ++i) {
		Point seg1 = q[i + 1] - q[i]; //the two segments that form the joint
		Point seg2 = q[i + 2] - q[i + 1];

		//move the base vectors along
		auto rot = getRotation(seg1, seg2);
		base1 = rot * base1;
		base2 = rot * base2;
		bases.push_back(make_pair(base1, base2));
	}

	return bases;
}

//call calculation
Meltable::Meltable(ArcLengthBSpline& spline, MeltableConstraints& constraints) :
		constraints(constraints), _length(spline.getLength()) {
	Calculation cal(this, spline);
	q = cal.run();
	segmentBases = Meltable::calculateBases(q);
}
