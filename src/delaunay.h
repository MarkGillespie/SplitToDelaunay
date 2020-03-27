#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"
#include "my_assert.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <deque>
#include <set>

using std::cout;
using std::endl;
using namespace geometrycentral;
using namespace geometrycentral::surface;

inline std::string toString(size_t n, size_t digits) {
    size_t existingDigits = ceil(log10(n + 0.01));
    my_assert(existingDigits < digits, std::to_string(n) +
                                           " is too big to fit in " +
                                           std::to_string(digits) + " digits");
    size_t nZeros = digits - existingDigits;
    return std::string(nZeros, '0') + std::to_string(n);
}

// Intervals must have lower end as first component.
// If interval.first > interval.second, the interval is treated as empty
using Interval = std::pair<double, double>;
using Disk     = std::pair<Vector2, double>;


/*
 * Class to make a mesh delaunay via edge splits. At the moment, it only works
 * on meshes without boundary
 *
 */
class Splitter {
  public:
    Splitter(VertexPositionGeometry& geo_);

    // Make the input geometry Delaunay using edge splits
    // Return which edge of the original mesh each edge lies on
    // New edges which split faces (and thus do not lie on any original mesh
    // edge) have value -1
    EdgeData<int> splitGeometry(bool verbose = false);

    void flipFlatEdgesToDelaunay(std::set<Edge> edges);
    void flipAllFlatEdgesToDelaunay();

    VertexPositionGeometry& geo;
    double rho_v, rho_e;
    EdgeData<std::vector<double>> edgeSamplePoints;
    EdgeData<char> flatEdge;
    EdgeData<int> parentIndex;

    std::pair<Halfedge, Halfedge> splitEdge(Edge e);

    bool isDelaunay(Edge e);

    Interval computeSplitInterval(Edge e);
    Interval computeBoundarySplitInterval(Edge e);

    static std::vector<std::vector<size_t>> intervalCombinations;
    static std::vector<std::vector<size_t>> genIntervalCombinations();
    static bool empty(Interval i);
    static double distTo(Interval i, double d);

    // Staps little to big's closes endpoint and then returns the interval big -
    // little
    static Interval complement(Interval little, Interval big);

    static Interval intersect(Interval a, Interval b);
    static Interval intersectIndices(const std::vector<size_t>& inds,
                                     const std::array<Interval, 5>& I);
    static size_t closestToInterval(Interval I, const std::vector<double>& pts);

    // Computes the intersection of the disk with the x-axis. Represents the
    // resulting interval with its lower end first. Returns (1, -1) if the
    // intersection is empty
    static Interval diskInterval(Disk disk);
    static Disk circumcenter(Vector2 v1, Vector2 v2, Vector2 v3);
    std::vector<double> edgePoints(Edge e);


    /*
     *  v5_____________v1___________ v4
     *    \           / \           /
     *     \         /   \         /
     *      \   D   /     \   C   /
     *       \     /   A   \     /
     *        \   /         \   /
     *         \ /           \ /
     *        v2 ----- e --->-  v0
     *         / \           / \
     *        /   \         /   \
     *       /     \   B   /     \
     *      /   E   \     /   F   \
     *     /         \   /         \
     *    /___________\ /___________\
     *   v6           v3             v7
     */
    // Compute the intervals corresponding to the circumdisks of
    // triangles in the neighborhood, but not incident on edge e (Triangles C,
    // D, E, F) in the butterfly
    // Some of these triangles might not exist of A or B are boundary faces
    // In that case, I give empty intervals
    std::array<Interval, 5> computeInteriorIntervals(Edge e);
    std::array<Interval, 3> computeBoundaryIntervals(Edge e);

    std::array<Vector2, 8> layOutButterfly(Edge e);
    std::array<Vector2, 5> layOutCroissant(Edge e);
};
