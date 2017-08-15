#pragma once

#include <functional>

#include "util.h"
#include "igl/viewer/ViewerPlugin.h"

class Handle: public igl::viewer::ViewerPlugin {
public:
	Handle() {
		viewer = nullptr;
		reference = nullptr;
		enabled = false;
		activeAxis = NONE;
		handleSize = 0.5f;
	}

	bool isEnabled() const {
		return enabled;
	}

	void setEnabled(bool enabled) {
		this->enabled = enabled;
		activeAxis = NONE;
	}

	void setReference(Point* reference) {
		if (activeAxis == NONE) { //do not allow change while editing
			this->reference = reference;
		}
	}

	//is the handle enabled and has it a camera and point reference to work with?
	bool isReadyForUpdate() {
		return isEnabled() && reference != nullptr && viewer != nullptr;
	}

	bool mouse_move(int mouse_x, int mouse_y);
	bool mouse_down(int button, int modifier);
	bool mouse_up(int button, int modifier);
	bool pre_draw();

	void setHandleSize(float handleSize) {
		this->handleSize = handleSize;
	}

private:
	enum Axis {
		NONE, AXIS_X, AXIS_Y, AXIS_Z
	};

	Point* reference;
	bool enabled;
	Axis activeAxis;
	float handleSize;
	int prev_mouse_x = 0, prev_mouse_y = 0;

	Point axisVector(Axis axis);
};
