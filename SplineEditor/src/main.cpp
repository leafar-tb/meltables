#include <iostream>

#include "igl/viewer/Viewer.h"
#include "nanogui/formhelper.h"
#include "igl/project.h"

#include "handle.hpp"
#include "meltable.h"
#include "branched_meltable.h"
#include "util.h"

int main(int argc, char * argv[]) {
	using namespace std;
	using namespace igl;

	igl::viewer::Viewer viewer;
	viewer.core.orthographic = true;
	viewer.core.show_lines = false;

	Polyline initPL;
	initPL.push_back(Point(0, 0, 0));
	initPL.push_back(Point(1, 1, 0));
	initPL.push_back(Point(2, -1, 0));
	initPL.push_back(Point(3, 0, 0));
	initPL.push_back(Point(4, -2, 0));
	initPL.push_back(Point(5, 1, 0));

	BranchedMeltable branched(initPL);
	MeltableConstraints& constraints = branched.constraints; //FIXME shouldn't be modifiable after sprouting

	// save the meltable mesh
	const auto exportMeltable = [&] () {
		branched.generateMesh();
		saveMesh("meltable", branched.exportMesh);
	};

	// widget for editing the polylines
	Handle handle;
	viewer.plugins.push_back(&handle);
	// show/hide the handles
	const auto toggleHandles = [&] () {
		handle.setEnabled(!handle.isEnabled());
		if(!handle.isEnabled()) { /*after disabling the handles*/
			branched.applyConstraints();
		}
	};

	bool showExport = false;
	bool showPreview = false;
	bool showEditor = true;
	// refresh displayed data
	const auto & redraw = [&]() {
		viewer.data.clear();
		viewer.data.face_based = true;
		if(showEditor || handle.isEnabled()) {
			branched.draw(viewer);
		}
		if(!handle.isEnabled()) {
			if(showPreview) {
				viewer.data.set_mesh(branched.previewMesh.vertices, branched.previewMesh.faces);
			} else if(showExport) {
				viewer.data.set_mesh(branched.exportMesh.vertices, branched.exportMesh.faces);
			}
		}
	};
	redraw();

	// init viewer menu
	viewer.callback_init = [&](igl::viewer::Viewer& viewer)->bool {
		/* create a new window, we do not need the standard gui */
		viewer.ngui->window()->dispose();
		viewer.ngui->addWindow(Eigen::Vector2i(10,10),"Menu");

		viewer.ngui->addGroup("Shape Constraints");
		{
			auto segmentsWidget = viewer.ngui->addVariable("Number of Segments", constraints.numSegments);
			segmentsWidget->setMinValue(1);
			segmentsWidget->setSpinnable(true);
			auto minTWidget = viewer.ngui->addVariable("Minimum Joint Angle", constraints.minTheta);
			minTWidget->setMinMaxValues(0, 160);
			minTWidget->setSpinnable(true);
			minTWidget->setUnits("°");
			minTWidget->setValueIncrement(5);
			auto maxTWidget = viewer.ngui->addVariable("Maximum Joint Angle", constraints.maxTheta);
			maxTWidget->setMinMaxValues(0, 160);
			maxTWidget->setSpinnable(true);
			maxTWidget->setValueIncrement(5);
			maxTWidget->setUnits("°");
		}

		viewer.ngui->addGroup("Mesh Constraints");
		{
			auto lengthWidget = viewer.ngui->addVariable("Length", constraints.targetLength);
			lengthWidget->setMinValue(0.0001);
			auto radiusWidget = viewer.ngui->addVariable("Radius", constraints.radius);
			radiusWidget->setMinValue(0.0001);
			auto thicknessWidget = viewer.ngui->addVariable("Joint Thickness", constraints.jointThickness);
			thicknessWidget->setMinValue(0.0001);
		}

		viewer.ngui->addGroup("Quality");
		{
			auto samplingWidget = viewer.ngui->addVariable("Sampling LOD", constraints.samplingLOD);
			samplingWidget->setMinValue(50);
			samplingWidget->setSpinnable(true);
			samplingWidget->setValueIncrement(25);
			auto cylinderWidget = viewer.ngui->addVariable("Cylinder LOD", constraints.cylinderLOD);
			cylinderWidget->setMinValue(8);
			cylinderWidget->setSpinnable(true);
		}

		viewer.ngui->addGroup("Actions");
		{
			viewer.ngui->addButton("Toggle Edit-Mode", toggleHandles);
			viewer.ngui->addButton("Add Control-Point", [&]() {branched.addPoint(); redraw();});
			viewer.ngui->addButton("Remove Control-Point", [&]() {branched.removePoint(); redraw();});
			viewer.ngui->addButton("Refresh Meltable", [&]() {branched.refreshMeltable(); redraw();});
			viewer.ngui->addButton("Refresh Meshes", [&]() {branched.generateMesh(); redraw();});
			viewer.ngui->addVariable<bool>("Show Export Mesh",
					[&](bool val) {
						showExport = val;
						redraw();
					},[&]() {
						return showExport;
					});
			viewer.ngui->addVariable<bool>("Show Preview Mesh",
					[&](bool val) {
						showPreview = val;
						redraw();
					},[&]() {
						return showPreview;
					});
			viewer.ngui->addVariable<bool>("Show Editor",
					[&](bool val) {
						showEditor = val;
						redraw();
					},[&]() {
						return showEditor;
					});
			viewer.ngui->addButton("Export Mesh", exportMeltable);
		}

		static bool right = true;
		static int segment = 0;
		/*create a second window for the branching menu*/
		viewer.ngui->addWindow(Eigen::Vector2i(viewer.core.viewport(2)-200, 10), "Branching");
		{
			auto segWidget = viewer.ngui->addVariable("Segment", segment);
			segWidget->setMinValue(0);
			segWidget->setSpinnable(true);
			viewer.ngui->addVariable("Right", right);
			viewer.ngui->addButton("Sprout", [&]() {branched.sprout(segment, right); redraw();});
			viewer.ngui->addButton("Prune", [&]() {branched.prune(segment, right); redraw();});
			viewer.ngui->addButton("Prune All", [&]() {branched.pruneAll(); redraw();});
			viewer.ngui->addButton("Select Next", [&]() {branched.selectNext(); redraw();});
			viewer.ngui->addButton("Apply Constraints", [&]() {branched.applyConstraints(); redraw();});

			auto angleWidget = viewer.ngui->addVariable<double>("Branch Angle",
					[&](double val) {
						branched.setBranchingAngle(val);
						redraw();
					},[&]() {
						return branched.getBranchingAngle();
					});
			angleWidget->setMinMaxValues(0, 360);
			angleWidget->setSpinnable(true);
			angleWidget->setValueIncrement(5);
			angleWidget->setUnits("°");
		}

		viewer.screen->performLayout();
		return false;
	};

	// setup key-bindings
	viewer.callback_key_pressed = [&](igl::viewer::Viewer &viewer, unsigned char key, int mods)->bool {
		switch(key) {
			case 'h':
			toggleHandles();
			break;
			case 'e':
			exportMeltable();
			break;
			case 'm':
			branched.refreshMeltable();
			break;
			case 'r':
			//just trigger redraw
			break;
			default:
			return false;
		}
		redraw();
		return true;
	};

	// always select the point closest to the mouse for editing
	viewer.callback_mouse_move = [&](igl::viewer::Viewer &viewer, int mouseX, int mouseY)->bool {
		Point* closest;
		float minDist = 5000;
		Point2D mouse(mouseX, mouseY);
		for(auto it = branched.editablesBegin(); it != branched.editablesEnd(); ++it) {
			Point2D screenpoint = worldToScreenFlat(&viewer, *it);
			float dist = (mouse-Point2D(screenpoint(0),screenpoint(1))).norm();
			if(dist < minDist) {
				minDist = dist;
				closest = &(*it);
			}
		}
		handle.setReference(closest);
		if(handle.isEnabled()) {
			redraw();
		}
		return false;
	};

	viewer.launch();
}
