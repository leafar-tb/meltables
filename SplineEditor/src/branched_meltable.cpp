#include "branched_meltable.h"
#include <iostream>
#include <cmath>

#include "igl/copyleft/cgal/CSGTree.h"

using namespace std;

#define SPLINE_DEGREE 3
#define SELECTED_COLOR Eigen::RowVector3d(1.0, 0.0, 0.0)
#define NORMAL_COLOR Eigen::RowVector3d(0.0, 1.0, 0.0)

//---------------------------------
//        Helper/Utility
//---------------------------------
Point directionFromAngle(pair<Point, Point> basis, double angleDEG) {
	return basis.first * sin(DEG_TO_RAD * angleDEG) + basis.second * cos(DEG_TO_RAD * angleDEG);
}

//length up to the mid of a segment of the meltable
double lengthToSegmentMid(Meltable* m, int segment) {
	double l = m->segment(segment).norm() / 2;
	for (int i = 0; i < segment; ++i) {
		l += m->segment(i).norm();
	}
	return l;
}

//rotate the given range around the first point by the given matrix
void rotate(Polyline::iterator begin, Polyline::iterator end, const Eigen::Matrix3d& rot) {
	for (auto it = begin; it != end; ++it) {
		*it = (rot * (*it - *begin)) + *begin;
	}
}

void translate(Polyline::iterator begin, Polyline::iterator end, const Point& dP) {
	for (auto it = begin; it != end; ++it) {
		*it = *it + dP;
	}
}

//modify the polyline, such that the angle constraints and constraints on joint width are met
//does change the polyline
void enforceConstraints(Polyline* plP, const MeltableConstraints& constraints) {
	Polyline& pl = *plP;
	double plLength = length(&pl);
	double plScale = constraints.targetLength / plLength;

	for (int i = 0; i < pl.size() - 2; ++i) {
		Point seg1 = pl[i + 1] - pl[i];
		Point seg2 = pl[i + 2] - pl[i + 1];
		double hingeAngle = angle(seg1, seg2);

		//---------------------
		//Step 1: adjust angles
		//---------------------
		Point rotAxis; //the axis we will rotate around
		if (seg1.normalized().dot(seg2.normalized()) > 0.999) { // if the segments are (almost) parallel
			rotAxis = getOrthogonal(seg1); //choose any axis
			hingeAngle = 0;
		} else if (seg1.normalized().dot(seg2.normalized()) < -0.999) { //(opposite directions)
			rotAxis = getOrthogonal(seg1);
			hingeAngle = 180;
		} else {
			rotAxis = seg1.cross(seg2).normalized();
		}
		//assuming reasonable min/max
		double dA = clamp(hingeAngle, constraints.minTheta, constraints.maxTheta) - hingeAngle;
		auto rotation = Eigen::AngleAxisd(DEG_TO_RAD * dA, rotAxis).toRotationMatrix();
		rotate(pl.begin() + (i + 1), pl.end(), rotation);

		//update angle
		hingeAngle = clamp(hingeAngle, constraints.minTheta, constraints.maxTheta);
		//----------------------
		//Step 2: adjust lengths
		//----------------------
		double hingeWidth = (constraints.radius * 2 - constraints.jointThickness) * tan(DEG_TO_RAD * hingeAngle / 2);
		//joint too wide -> make segment longer
		if ((seg1.norm() * plScale) / 2 < hingeWidth) {
			double dL = hingeWidth - ((seg1.norm() * plScale) / 2);
			translate(pl.begin() + (i + 1), pl.end(), seg1.normalized() * dL);
		}
		if ((seg2.norm() * plScale) / 2 < hingeWidth) {
			double dL = hingeWidth - ((seg2.norm() * plScale) / 2);
			translate(pl.begin() + (i + 2), pl.end(), seg2.normalized() * dL);
		}
	}
}

//---------------------------------
//          WorkingSet
//---------------------------------
BranchedMeltable::WorkingSet::WorkingSet(Polyline& pl) :
		pl(pl) {
}

void BranchedMeltable::WorkingSet::draw(igl::viewer::Viewer& viewer, const Eigen::RowVector3d& color) {
	drawSpline(viewer, constructBSpline(pl, SPLINE_DEGREE), 40, Eigen::RowVector3d(0.0, 0.0, 1.0));
	drawPL(viewer, pl, Eigen::RowVector3d(0.0, 0.0, 0.7));
	if (meltable != nullptr) {
		auto mpl = meltable->getPoints();
		drawPL(viewer, mpl, color);
	}
}

void BranchedMeltable::WorkingSet::refreshMeltable(MeltableConstraints constraints) {
	delete (meltable);
	ArcLengthBSpline albs(pl, SPLINE_DEGREE);
	meltable = new Meltable(albs, constraints);
}

//---------------------------------
//        BranchedMeltable
//---------------------------------
BranchedMeltable::BranchedMeltable(Polyline& pl) :
		trunk(pl) {
	branchingAngle = 0;
	selectedPosition = std::make_pair(-1, true);
}

void BranchedMeltable::draw(igl::viewer::Viewer& viewer) {
	if (isTrunkSelected()) {
		trunk.draw(viewer, SELECTED_COLOR);
	} else {
		trunk.draw(viewer, NORMAL_COLOR);
	}

	for (auto branchEntry = branches.begin(); branchEntry != branches.end(); ++branchEntry) {
		if (selected() == &(branchEntry->second)) {
			drawPL(viewer, branchEntry->second, SELECTED_COLOR);
		} else {
			drawPL(viewer, branchEntry->second, NORMAL_COLOR);
		}
	}
}

void BranchedMeltable::sprout(int segment, bool rightBranch) {
	if (branches.empty()) {
		// no branches -> trunk is selected -> refresh before branching off it
		refreshMeltable();
	}
	// out of bounds check
	if (segment < 0 || segment >= trunk.meltable->numSegments()) {
		return;
	}

	BranchPosition pos = make_pair(segment, rightBranch);
	if (get(pos) != nullptr) //branch already exists
		return;

	// check if segment is long enough for branching FIXME does not consider joints
	if (trunk.meltable->segment(segment).norm() * targetScale() <= 2 * constraints.radius) {
		cerr << "segment is too short to support a branch";
		return;
	}

	Polyline pl; // starts at the midpoint of the segment and goes out straight
	pl.push_back(trunk.meltable->getPoints()[segment] + 0.5 * trunk.meltable->segment(segment));
	Point extension = directionFromAngle(trunk.meltable->segmentBase(segment),
			branchingAngle + (rightBranch ? 0 : 180));
	pl.push_back(pl.back() + extension * branchBaseLength());
	extension *= trunk.meltable->segment(segment).norm();
	for (int i = 0; i < 4; ++i) {
		pl.push_back(pl.back() + extension);
	}

	branches.insert(make_pair(pos, pl));
	if (selectedPosition.first == -1) { // if trunk was selected, change to the new branch
		selectedPosition = pos;
	}
	applyConstraints(pos);
}

void BranchedMeltable::setBranchingAngle(double branchingAngle) {
	double dAngle = branchingAngle - this->branchingAngle;
	//rotate all branches according to the change in angle
	for (auto entryI = branches.begin(); entryI != branches.end(); ++entryI) {
		Branch& branch = entryI->second;
		auto rotation = Eigen::AngleAxisd(DEG_TO_RAD * dAngle,
				trunk.meltable->segment(entryI->first.first).normalized()).toRotationMatrix();
		rotate(branch.begin(), branch.end(), rotation);
	}

	this->branchingAngle = branchingAngle;
}

void BranchedMeltable::generateMesh() {
	using igl::copyleft::cgal::CSGTree;
	using igl::MeshBooleanType;

	if (trunk.meltable == nullptr) {
		refreshMeltable(); //trunk is selected
	}

	//union of all branches
	CSGTree branchMeshes;
	for (auto entryI = branches.begin(); entryI != branches.end(); ++entryI) {
		BranchPosition pos = entryI->first;
		applyConstraints(pos); //update may change length a bit, so we get constraints afterwards
		Polyline& pl = entryI->second;
		//choose local base (for hinge direction) s.t. we can preserve orientation more easily
		Point by = -trunk.meltable->segment(pos.first).normalized(); //y along the segment we branch off
		Point bx = by.cross(pl[1] - pl[0]).normalized(); //z is along first segment of the branch and x=y.cross(z)
		Mesh branchMesh = Meltable::generateMesh(pl, Meltable::calculateBases(pl, make_pair(bx, by)),
				branchConstraints(pos));

		//need to rotate and move branches to correct position
		double angleRAD = DEG_TO_RAD * (branchingAngle + (pos.second ? 0 : 180));
		//first rotate around z, to align hinge direction
		rotate(branchMesh, Eigen::AngleAxisd(angleRAD, Point::UnitZ()).toRotationMatrix());
		//then tilt the branch from upright(along z) to sideways, pointing in the branching angle
		rotate(branchMesh, getRotation(Point::UnitZ(), Point(sin(angleRAD), cos(angleRAD), 0)));
		//we want the length according to target size
		double l = lengthToSegmentMid(trunk.meltable, pos.first) * targetScale();
		translate(branchMesh, Point(0, 0, l));
		branchMeshes = CSGTree(CSGTree(branchMesh.vertices, branchMesh.faces), branchMeshes,
				MeshBooleanType::MESH_BOOLEAN_TYPE_UNION);
	}
	//union branches with trunk
	Mesh trunkMesh = trunk.meltable->generateMesh();
	CSGTree result(CSGTree(trunkMesh.vertices, trunkMesh.faces), branchMeshes,
			MeshBooleanType::MESH_BOOLEAN_TYPE_UNION);
	exportMesh = Mesh(result.cast_V<Eigen::MatrixXd>(), result.F());

	//merge trunk and branch previews
	previewMesh = trunk.meltable->generatePreview();
	for (auto entryI = branches.begin(); entryI != branches.end(); ++entryI) {
		Polyline& pl = entryI->second;
		previewMesh = merge(previewMesh,
				Meltable::generateMesh(pl, Meltable::calculateBases(pl), branchConstraints(entryI->first), true));
	}
}

//------------
//------------

void BranchedMeltable::selectNext() {
	applyConstraints(); //leave branch in a valid state

	if (branches.empty()) {
		selectedPosition = make_pair(-1, true);
	} else {
		if (selectedPosition.second) { // if right is selected
			auto segBranches = get(selectedPosition.first);
			if (segBranches.second != nullptr) { // and there is a left
				selectedPosition.second = false; //select left
				return;
			}
		}

		size_t segmentCount = trunk.meltable->numSegments();
		for (int i = 1; i < segmentCount; ++i) {
			auto segBranches = get((selectedPosition.first + i) % segmentCount);
			if (segBranches.first != nullptr) {
				selectedPosition = make_pair((selectedPosition.first + i) % segmentCount, true);
				return;
			}
			if (segBranches.second != nullptr) {
				selectedPosition = make_pair((selectedPosition.first + i) % segmentCount, false);
				return;
			}
		}
	}
}

void BranchedMeltable::addPoint() {
	if (isTrunkSelected()) {
		trunk.pl.push_back(2 * trunk.pl.back() - trunk.pl[trunk.pl.size() - 2]);
	} else {
		selected()->push_back(2 * selected()->back() - (*selected())[selected()->size() - 2]);
	}

}

void BranchedMeltable::removePoint() {
	if (isTrunkSelected()) {
		if (trunk.pl.size() > SPLINE_DEGREE + 1) {
			trunk.pl.pop_back();
		}
	} else {
		if (selected()->size() > 2)
			selected()->pop_back();
	}
}

void BranchedMeltable::refreshMeltable() {
	if (isTrunkSelected()) {
		trunk.refreshMeltable(constraints);
	}
}

Polyline::iterator BranchedMeltable::editablesBegin() {
	if (isTrunkSelected())
		return trunk.pl.begin();
	else
		return selected()->begin() + 2;
}

Polyline::iterator BranchedMeltable::editablesEnd() {
	if (isTrunkSelected())
		return trunk.pl.end();
	else
		return selected()->end();
}

auto BranchedMeltable::get(int segment)-> pair<Branch*, Branch*> {
	return make_pair(get(make_pair(segment, true)) //
			, get(make_pair(segment, false)));
}

auto BranchedMeltable::get(const BranchPosition& pos)->Branch* {
	auto branch = branches.find(pos);
	return branch == branches.end() ? nullptr : &(branch->second);
}

auto BranchedMeltable::selected()-> Branch* {
	if (selectedPosition.first == -1)
		return nullptr;
	else
		return &(branches.find(selectedPosition)->second);
}

void BranchedMeltable::prune(int segment, bool rightBranch) {
	branches.erase(make_pair(segment, rightBranch));
	if (selectedPosition.first == segment && selectedPosition.second == rightBranch) {
		selectNext();
	}
}

void BranchedMeltable::pruneAll() {
	branches.clear();
	selectNext();
}

bool BranchedMeltable::isTrunkSelected() {
	return selectedPosition.first == -1;
}

double BranchedMeltable::targetScale() {
	return constraints.targetLength / trunk.meltable->length();
}

MeltableConstraints BranchedMeltable::branchConstraints(const BranchPosition& pos) {
	Branch* branch = get(pos);
	if (branch == nullptr)
		return constraints;

	MeltableConstraints bc = constraints;
	bc.targetLength = targetScale() * length(branch);
	return bc;
}

void BranchedMeltable::applyConstraints() {
	if (isTrunkSelected())
		return;
	else
		applyConstraints (selectedPosition);
}

void BranchedMeltable::applyConstraints(const BranchPosition& pos) {
	if (get(pos) != nullptr) { //may be a branch that was just pruned
		enforceConstraints(get(pos), branchConstraints(pos));
	}
}

double BranchedMeltable::branchBaseLength() {
	return 3 * constraints.radius / targetScale();
}
