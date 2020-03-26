#include "delaunay.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector3.h"
#include "test_utils.h"

class SplitterTest : public ::testing::Test {
  public:
    static std::unique_ptr<double> globalDouble;
    double localDouble;

  protected:
    static void SetUpTestSuite() {
        double* ten  = new double();
        *ten         = 10;
        globalDouble = std::unique_ptr<double>(ten);
    }

    void SetUp() override { localDouble = 10; }
};

std::unique_ptr<double> SplitterTest::globalDouble = nullptr;

/*
 *
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
TEST_F(SplitterTest, ButterflyLayoutWorks) {

    std::vector<std::vector<size_t>> faces(6);
    faces[0] = std::vector<size_t>{2, 0, 1};
    faces[1] = std::vector<size_t>{0, 2, 3};
    faces[2] = std::vector<size_t>{0, 4, 1};
    faces[3] = std::vector<size_t>{2, 1, 5};
    faces[4] = std::vector<size_t>{2, 6, 3};
    faces[5] = std::vector<size_t>{0, 3, 7};

    HalfedgeMesh mesh(faces);

    VertexData<Vector3> posData(mesh);
    posData[0] = Vector3{2, 0, 0};
    posData[1] = Vector3{1.5, 1, 0};
    posData[2] = Vector3{0, 0, 0};
    posData[3] = Vector3{0.75, -1.2, 0};
    posData[4] = Vector3{5, 5, 0};
    posData[5] = Vector3{-1, 2, 0};
    posData[6] = Vector3{-0.2, -0.1, 0};
    posData[7] = Vector3{3, -0.25, 0};

    Edge e = mesh.edge(0);
    EXPECT_TRUE(e.halfedge().vertex() == mesh.vertex(2));
    EXPECT_TRUE(e.halfedge().twin().vertex() == mesh.vertex(0));

    VertexPositionGeometry geo(mesh, posData);
    Splitter splitter(geo);

    std::array<Vector2, 8> butterfly = splitter.layOutButterfly(e);

    auto toPlane = [](Vector3 v) { return Vector2{v.x, v.y}; };
    for (size_t iV = 0; iV < 8; ++iV) {
        EXPECT_VEC2_NEAR(butterfly[iV], toPlane(posData[iV]), 1e-8);
    }
}

TEST_F(SplitterTest, circumcenterSameDistanceFromPoints) {
    Vector2 v1{1, 0};
    Vector2 v2{12, 6};
    Vector2 v3{8, 9};
    Vector2 circumcenter;
    double circumradius;

    std::tie(circumcenter, circumradius) = Splitter::circumcenter(v1, v2, v3);

    EXPECT_NEAR((v1 - circumcenter).norm(), circumradius, 1e-8);
    EXPECT_NEAR((v2 - circumcenter).norm(), circumradius, 1e-8);
    EXPECT_NEAR((v3 - circumcenter).norm(), circumradius, 1e-8);
}

TEST_F(SplitterTest, diskLineIntersection) {
    Disk d     = std::make_pair(Vector2{1, 1}, 2);
    Interval i = Splitter::diskInterval(d);

    EXPECT_NEAR(i.first, 1 - sqrt(3), 1e-8);
    EXPECT_NEAR(i.second, 1 + sqrt(3), 1e-8);

    d = std::make_pair(Vector2{0, 2}, 1);
    i = Splitter::diskInterval(d);

    EXPECT_TRUE(Splitter::empty(i));
}
