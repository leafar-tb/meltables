#include "meltable.h"

#include <cmath>
#include <iostream>
#include <Eigen/Geometry>

#include "igl/copyleft/cgal/CSGTree.h"

using namespace std;

Mesh cylinder(int polyCount, double length, double radius = .5) {
	Mesh mesh(2 * polyCount + 2, 4 * polyCount);
	if (polyCount < 3) {
		polyCount = 3;
	}

	for (int i = 0; i < polyCount; ++i) {
		//top/bottom circle
		double angle = i * 2 * M_PI / polyCount;
		mesh.vertices.row(i) << cos(angle) * radius, sin(angle) * radius, 0.;
		mesh.vertices.row(i + polyCount) << cos(angle) * radius, sin(angle) * radius, length;
	}
	//center vertices
	mesh.vertices.row(2 * polyCount) << 0., 0., 0.;
	mesh.vertices.row(2 * polyCount + 1) << 0., 0., length;

	for (int i = 0; i < polyCount; ++i) {
		//top/bottom faces
		mesh.faces.row(i) << 2 * polyCount, (i + 1) % polyCount, i;
		mesh.faces.row(i + polyCount) << 2 * polyCount + 1, i + polyCount, ((i + 1) % polyCount) + polyCount;
		//side faces
		mesh.faces.row(i + 2 * polyCount) << i, (i + 1) % polyCount, i + polyCount;
		mesh.faces.row(i + 3 * polyCount) << (i + 1) % polyCount, ((i + 1) % polyCount) + polyCount, i + polyCount;
	}

	return mesh;
}

Mesh wedge(double angle, double size = 1) {
	Mesh mesh(6, 8);
	mesh.vertices.row(0) << -size / 2, size / 2, 0;
	mesh.vertices.row(1) << -size / 2, -size / 2, 0;

	double halfOpening = tan(DEG_TO_RAD * angle / 2) * size;
	mesh.vertices.row(2) << size / 2, size / 2, halfOpening;
	mesh.vertices.row(3) << size / 2, size / 2, -halfOpening;
	mesh.vertices.row(4) << size / 2, -size / 2, halfOpening;
	mesh.vertices.row(5) << size / 2, -size / 2, -halfOpening;

	//sides
	mesh.faces.row(0) << 0, 2, 3;
	mesh.faces.row(1) << 1, 5, 4;
	//top/bottom
	mesh.faces.row(2) << 0, 4, 2;
	mesh.faces.row(3) << 1, 4, 0;
	mesh.faces.row(4) << 0, 3, 1;
	mesh.faces.row(5) << 1, 3, 5;
	//front
	mesh.faces.row(6) << 3, 2, 4;
	mesh.faces.row(7) << 5, 3, 4;

	return mesh;
}

double plLength(const Polyline& pl) {
	return length(pl);
}

//half the width of the joint in the mid of the rod
double jointMidHalfWidth(double angle, const MeltableConstraints& constraints) {
	return (constraints.radius - constraints.jointThickness) * tan(DEG_TO_RAD * angle / 2);
}

Mesh Meltable::generateMesh(const Polyline& q, const vector<pair<Point, Point> >& segmentBases,
		const MeltableConstraints& constraints, bool preview) {
	using igl::copyleft::cgal::CSGTree;
	using igl::MeshBooleanType;

	cout << "generating mesh " << flush;
	vector<double> turnAngles; //in which direction a joint faces; in RAD
	vector<double> jointAngles; //the closing angle of a joint; in DEG
	vector<double> partialLengths; //summed target-lengths of the segments up to an index
	double targetScaling = constraints.targetLength / plLength(q);

	Point segment0 = q[1] - q[0];
	partialLengths.push_back(segment0.norm() * targetScaling);

	for (size_t i = 0; i < q.size() - 2; ++i) {
		Point seg1 = q[i + 1] - q[i]; //the two segments that form the joint
		Point seg2 = q[i + 2] - q[i + 1];
		jointAngles.push_back(angle(seg1, seg2));

		//get the angle of seg2 in this system
		auto basis = segmentBases[i];
		turnAngles.push_back(getAngleInBasis(seg2, basis.first, basis.second));

		//we want to achieve targetLength along the mid of the rod
		//, so we need a correction for the length 'lost' when cutting out the joints
		//we only use this for the preview mesh, to match the curve, otherwise targetLength means length of the rod
		double correction = preview ? jointMidHalfWidth(jointAngles.back(), constraints) : 0;
		partialLengths.back() += correction;
		partialLengths.push_back(partialLengths.back() + correction + targetScaling * seg2.norm());
	}
	cout << "." << flush;

	CSGTree wedges;
	//double targetScaling = constraints.targetLength / partialLengths.back();
	for (size_t seg = 0; seg < turnAngles.size(); ++seg) {
		Mesh wedgeM = wedge(jointAngles[seg], constraints.radius * 2);

		//translate to proper height and move a bit to leave the joint connected
		translate(wedgeM, Point(constraints.jointThickness, 0, partialLengths[seg]));
		//rotate into the correct direction
		rotate(wedgeM, Eigen::AngleAxisd(turnAngles[seg], Point::UnitZ()).toRotationMatrix());

		wedges = CSGTree(CSGTree(wedgeM.vertices, wedgeM.faces), wedges, MeshBooleanType::MESH_BOOLEAN_TYPE_UNION);
	}
	cout << "." << flush;

	if (preview) { //preview of the folded object: cut the single segments and combine them
		Mesh result;
		for (size_t seg = 0; seg < partialLengths.size(); ++seg) {
			double from = (seg == 0 ? 0 : partialLengths[seg - 1]);
			double to = partialLengths[seg];
			Mesh tmpMesh = cylinder(constraints.cylinderLOD, to - from, constraints.radius);
			translate(tmpMesh, Point(0, 0, from)); //move segment upwards
			//cut out the wedges
			CSGTree segment(CSGTree(tmpMesh.vertices, tmpMesh.faces), wedges, MeshBooleanType::MESH_BOOLEAN_TYPE_MINUS);
			//re-extract into a mesh
			tmpMesh = Mesh(segment.cast_V<Eigen::MatrixXd>(), segment.F());
			translate(tmpMesh,
					Point(0, 0, -from - jointMidHalfWidth(seg == 0 ? 0 : jointAngles[seg - 1], constraints))); //move it back down
					//(with a correction for the cut out joint)
			//to rotate it into the right direction
			//first align z
			Eigen::Matrix3d rotation = getRotation(Point::UnitZ(), q[1] - q[0]);
			//rotate around new z, to align x (and y)
			double ra = getAngleInBasis(segmentBases[0].first, rotation * Point::UnitX(), rotation * Point::UnitY());
			rotation = Eigen::AngleAxisd(ra, rotation * Point::UnitZ()).toRotationMatrix() * rotation;
			//now apply the rotations along the previous segments
			for (size_t i = 0; i < seg; ++i) {
				rotation = getRotation(q[i + 1] - q[i], q[i + 2] - q[i + 1]) * rotation;
			}
			rotate(tmpMesh, rotation);
			tmpMesh.vertices /= targetScaling; //rescale to editor-size
			translate(tmpMesh, q[seg]); //move it to its actual position
			result = merge(result, tmpMesh); //merge into result
		}
		cout << ". done" << endl;
		return result;
	} else { //export for printing: one rod with the wedges cut out
		Mesh cyl = cylinder(constraints.cylinderLOD, constraints.targetLength, constraints.radius);
		CSGTree rod(CSGTree(cyl.vertices, cyl.faces), wedges, MeshBooleanType::MESH_BOOLEAN_TYPE_MINUS);
		cout << ". done" << endl;

		return Mesh(rod.cast_V<Eigen::MatrixXd>(), rod.F());
	}
}

Mesh Meltable::generateMesh() {
	return Meltable::generateMesh(q, segmentBases, constraints);
}

Mesh Meltable::generatePreview() {
	return Meltable::generateMesh(q, segmentBases, constraints, true);
}
