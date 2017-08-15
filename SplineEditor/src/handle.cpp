#include "handle.hpp"

#include <GL/glut.h>
#include <Eigen/Core>
#include <iostream>

#include "igl/viewer/Viewer.h"
#include "igl/project_to_line.h"

bool Handle::pre_draw() {
	if (isReadyForUpdate()) {
		viewer->data.add_edges((*reference - axisVector(AXIS_X)).transpose().cast<double>(),
				(*reference + axisVector(AXIS_X)).transpose().cast<double>(), Eigen::RowVector3d(1, 0, 0));
		viewer->data.add_edges((*reference - axisVector(AXIS_Y)).transpose().cast<double>(),
				(*reference + axisVector(AXIS_Y)).transpose().cast<double>(), Eigen::RowVector3d(0, 1, 0));
		viewer->data.add_edges((*reference - axisVector(AXIS_Z)).transpose().cast<double>(),
				(*reference + axisVector(AXIS_Z)).transpose().cast<double>(), Eigen::RowVector3d(0, 0, 1));
	}
	return false;
}

bool Handle::mouse_move(int mouse_x, int mouse_y) {
	if (isReadyForUpdate()) {
		if (activeAxis == NONE) {
			return false;
		}

		//project the axis to a vector in screen space
		Point2D axisOnScreen = worldToScreenFlat(viewer, *reference + axisVector(activeAxis))
				- worldToScreenFlat(viewer, *reference);
		Point2D mouseDelta = Point2D(mouse_x, mouse_y) - Point2D(prev_mouse_x, prev_mouse_y);
		float worldDeltaScale = mouseDelta.dot(axisOnScreen) / axisOnScreen.squaredNorm();
		*reference = *reference + worldDeltaScale * axisVector(activeAxis);

		prev_mouse_x = mouse_x;
		prev_mouse_y = mouse_y;
		return true;
	}
	return false;
}

bool Handle::mouse_down(int button, int modifier) {
	if (isReadyForUpdate() && button == GLUT_LEFT_BUTTON) {
		Axis nearestAxis = NONE;
		double minScreenDist = 1000; //clicks too far away will be ignored
		Point2D mouse(viewer->current_mouse_x, viewer->current_mouse_y);
		Point2D axisStart, axisEnd;
		double t, sqrDist;

		axisStart = worldToScreenFlat(viewer, *reference + axisVector(AXIS_X));
		axisEnd = worldToScreenFlat(viewer, *reference - axisVector(AXIS_X));
		igl::project_to_line(mouse(0), mouse(1), 0. //
				, axisStart(0), axisStart(1), 0. //
				, axisEnd(0), axisEnd(1), 0. //
				, t, sqrDist);
		if (sqrDist < minScreenDist) {
			nearestAxis = AXIS_X;
			minScreenDist = sqrDist;
		}

		axisStart = worldToScreenFlat(viewer, *reference + axisVector(AXIS_Y));
		axisEnd = worldToScreenFlat(viewer, *reference - axisVector(AXIS_Y));
		igl::project_to_line(mouse(0), mouse(1), 0. //
				, axisStart(0), axisStart(1), 0. //
				, axisEnd(0), axisEnd(1), 0. //
				, t, sqrDist);
		if (sqrDist < minScreenDist) {
			nearestAxis = AXIS_Y;
			minScreenDist = sqrDist;
		}

		axisStart = worldToScreenFlat(viewer, *reference + axisVector(AXIS_Z));
		axisEnd = worldToScreenFlat(viewer, *reference - axisVector(AXIS_Z));
		igl::project_to_line(mouse(0), mouse(1), 0. //
				, axisStart(0), axisStart(1), 0. //
				, axisEnd(0), axisEnd(1), 0. //
				, t, sqrDist);
		if (sqrDist < minScreenDist) {
			nearestAxis = AXIS_Z;
			//no longer needed minScreenDist = sqrDist;
		}

		activeAxis = nearestAxis;
		prev_mouse_x = viewer->current_mouse_x;
		prev_mouse_y = viewer->current_mouse_y;
		return activeAxis != NONE;
	}
	return false;
}

bool Handle::mouse_up(int button, int modifier) {
	if (button == GLUT_LEFT_BUTTON) {
		activeAxis = NONE;
		return true;
	}
	return false;
}

Point Handle::axisVector(Axis axis) {
	switch (axis) {
	case AXIS_X:
		return Point(handleSize, 0, 0);
	case AXIS_Y:
		return Point(0, handleSize, 0);
	case AXIS_Z:
		return Point(0, 0, handleSize);
	case NONE:
		return Point::Zero();
	}
	return Point::Zero();
}
