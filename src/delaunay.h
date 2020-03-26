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

class Splitter {
  public:
    Splitter(VertexPositionGeometry& geo_);
    void splitGeometry(bool verbose = false);

    std::vector<Edge> flipFlatEdgesToDelaunay(std::set<Edge> edges);

    VertexPositionGeometry& geo;
    double rho_v, rho_e;
    EdgeData<std::vector<double>> edgeSamplePoints;
    EdgeData<char> flatEdge;

    std::pair<Halfedge, Halfedge> splitEdge(Edge e);

    bool isDelaunay(Edge e);
    static std::vector<std::vector<size_t>> intervalCombinations;
    static std::vector<std::vector<size_t>> genIntervalCombinations();
    static bool empty(Interval i);
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
    std::array<Vector2, 8> layOutButterfly(Edge e);
};
