#include "util.h"
#include "igl/project.h"
#include "igl/writeOBJ.h"

double length(const Polyline& pl) {
	return length(&pl);
}

double length(const Polyline* pl) {
	double result = 0;
	for (int i = 0; i < pl->size() - 1; ++i) {
		result += (pl->at(i) - pl->at(i + 1)).norm();
	}
	return result;
}

Mesh& translate(Mesh& mesh, const Point& dP) {
	auto row = dP.transpose();
	mesh.vertices.rowwise() += row;
	return mesh;
}

Mesh& rotate(Mesh& mesh, const Eigen::Matrix3d& rot) {
	auto rotT = rot.transpose();
	for (int row = 0; row < mesh.vertices.rows(); ++row) {
		mesh.vertices.row(row) *= rotT;
	}
	return mesh;
}

void saveMesh(std::string name, const Mesh& mesh) {
	igl::writeOBJ("meltable.obj", mesh.vertices, mesh.faces);
}

Mesh merge(const Mesh& m1, const Mesh& m2) {
	Mesh result(m1.vertices.rows() + m2.vertices.rows(), m1.faces.rows() + m2.faces.rows());
	//concatenate all vertex rows
	result.vertices << m1.vertices, m2.vertices;
	//concatenate all face rows (shift the indices of m2)
	result.faces << m1.faces, m2.faces.array() + m1.vertices.rows();
	return result;
}

Point worldToScreen(igl::viewer::Viewer *viewer, Point worldPoint) {
	Point tmp = igl::project(Eigen::Vector3f(worldPoint(0), worldPoint(1), worldPoint(2)),
			(Eigen::Matrix<float, 4, 4>)(viewer->core.view * (viewer->core.model)), viewer->core.proj,
			viewer->core.viewport).cast<double>();
	return Point(tmp(0), viewer->core.viewport(3) - tmp(1), tmp(2));
}

Point2D worldToScreenFlat(igl::viewer::Viewer *viewer, Point worldPoint) {
	Point tmp = worldToScreen(viewer, worldPoint);
	return Point2D(tmp(0), tmp(1));
}

double angle(Point vec1, Point vec2) {
	return RAD_TO_DEG * acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm()));
}

Point getOrthogonal(Point& vec) {
	if (vec(0) != 0 || vec(1) != 0) {
		return Point(vec(1), -vec(0), 0).normalized();
	} else if (vec(2) != 0) {
		return Point(0, vec(2), 0).normalized();
	} else {
		return Point::Zero();
	}
}

Eigen::Matrix3d getRotation(Point from, Point to) {
	Point axis = from.cross(to);
	//if they are collinear, choose any orthogonal axis
	if (axis.norm() == 0) {
		axis = getOrthogonal(from);
	}

	Point b1 = from.normalized();
	Point b2 = axis.cross(b1).normalized();
	double angle = getAngleInBasis(to, b1, b2);
	return Eigen::AngleAxisd(angle, axis.normalized()).toRotationMatrix();
}

double getAngleInBasis(const Point& v, const Point& b1, const Point& b2) {
	return atan2(b2.dot(v), b1.dot(v));
}
