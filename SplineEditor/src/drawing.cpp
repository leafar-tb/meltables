#include "drawing.h"
#include "bspline.h"

using Eigen::Vector3d;

void drawSpline(igl::viewer::Viewer& viewer, const ts::BSpline& spline, size_t segments,
		const Eigen::RowVector3d& color) {
	float kMin = spline.knots().front();
	float kMax = spline.knots().back();

	Point lastPoint = evaluate(spline, kMin);
	for (size_t i = 1; i <= segments; ++i) {
		Point nextPoint = evaluate(spline, kMin + (i / (double) segments) * (kMax - kMin));
		viewer.data.add_edges(lastPoint.transpose().cast<double>(), nextPoint.transpose().cast<double>(), color);
		lastPoint = nextPoint;
	}
}

void drawPL(igl::viewer::Viewer & viewer, const Polyline &pl) {
	drawPL(viewer, pl, Eigen::RowVector3d(0.1, 0.1, 0.1));
}

void drawPL(igl::viewer::Viewer& viewer, const Polyline& pl, const Eigen::RowVector3d& color) {
	for (size_t i = 0; i < pl.size() - 1; ++i) {
		viewer.data.add_edges(pl[i].transpose().cast<double>(), pl[i + 1].transpose().cast<double>(), color);
	}
}
