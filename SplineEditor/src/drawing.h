#include <Eigen/Core>

#include "tinyspline/tinysplinecpp.h"
#include "util.h"
#include "igl/viewer/Viewer.h"

void drawSpline(igl::viewer::Viewer& viewer, const ts::BSpline& spline, size_t segments = 20,
		const Eigen::RowVector3d& color = Eigen::RowVector3d(1.0, 0.0, 0.0));
void drawPL(igl::viewer::Viewer & viewer, const Polyline &pl);
void drawPL(igl::viewer::Viewer & viewer, const Polyline &pl, const Eigen::RowVector3d& color);
