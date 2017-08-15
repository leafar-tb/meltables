#include <vector>
#include <utility>
#include <map>

#include "igl/viewer/Viewer.h"

#include "util.h"
#include "bspline.h"
#include "meltable.h"
#include "drawing.h"

// A branched meltable consists of several meltables, which are organized in so called working sets.
// A working set stores the information necessary to work with a single meltable, such as the underlying spline.
// One working set forms the basis or trunk of the branched meltable; other meltables can branch off it.
// Branching is (for now) only possible in the middle of a trunk-segment and can be to the left or right.
// The direction of left and right are controlled by the branchingAngle.
//
// One of the working sets is considered selected and can be modified via the respective methods.
// The trunk can only be modified, if there are no branches on it.
class BranchedMeltable {
public:
	BranchedMeltable(Polyline& pl);

	void draw(igl::viewer::Viewer& viewer);
	void generateMesh();

	void sprout(int segment, bool rightBranch);
	void prune(int segment, bool rightBranch);
	void pruneAll();

	// set the next branch to active
	void selectNext();

	void refreshMeltable();
	void addPoint();
	void removePoint();

	//to a branch
	void applyConstraints();

	// iterator over the editable points of the currently selected working set
	Polyline::iterator editablesBegin();
	// end-iterator of the editable points of the currently selected working set
	Polyline::iterator editablesEnd();

	MeltableConstraints constraints;

	double getBranchingAngle() const {
		return branchingAngle;
	}

	void setBranchingAngle(double branchingAngle);

	Mesh exportMesh;
	Mesh previewMesh;

private:
	struct WorkingSet {
		WorkingSet(Polyline& pl);

		void draw(igl::viewer::Viewer& viewer, const Eigen::RowVector3d& color);
		void refreshMeltable(MeltableConstraints constraints);

		Polyline pl;
		Meltable* meltable = nullptr;
	};

	//at which segment and if it is to the right
	typedef std::pair<int, bool> BranchPosition;
	typedef Polyline Branch;

	WorkingSet trunk;

	//in DEG
	double branchingAngle;

	std::map<BranchPosition, Branch> branches;

	BranchPosition selectedPosition;

	//constraints.targetLength / trunk.length
	double targetScale();

	//minimal length of the first branch segment (s.t. it doesn't intersect with the trunk)
	double branchBaseLength();

	MeltableConstraints branchConstraints(const BranchPosition& pos);
	void applyConstraints(const BranchPosition& pos);

	std::pair<Branch*, Branch*> get(int segment);
	Branch* get(const BranchPosition& pos);
	Branch* selected();
	bool isTrunkSelected();
};
