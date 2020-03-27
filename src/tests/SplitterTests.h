#include "Interval.h"
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


TEST_F(SplitterTest, circumcenterSameDistanceFromPoints) {
    Vector2 v1{1, 0};
    Vector2 v2{12, 6};
    Vector2 v3{8, 9};
    Vector2 circumcenter;
    double circumradius;

    std::tie(circumcenter, circumradius) = circumcircle(v1, v2, v3);

    EXPECT_NEAR((v1 - circumcenter).norm(), circumradius, 1e-8);
    EXPECT_NEAR((v2 - circumcenter).norm(), circumradius, 1e-8);
    EXPECT_NEAR((v3 - circumcenter).norm(), circumradius, 1e-8);
}

TEST_F(SplitterTest, diskLineInterval) {
    Disk d     = std::make_pair(Vector2{1, 1}, 2);
    Interval i = diskInterval(d);

    EXPECT_NEAR(i.first, 1 - sqrt(3), 1e-8);
    EXPECT_NEAR(i.second, 1 + sqrt(3), 1e-8);

    d = std::make_pair(Vector2{0, 2}, 1);
    i = diskInterval(d);

    EXPECT_TRUE(empty(i));
}

TEST_F(SplitterTest, IntervalComplement) {
    Interval big    = std::make_pair(-1, 2);
    Interval little = std::make_pair(0, 2);
    Interval comp   = complement(little, big);
    EXPECT_NEAR(comp.first, -1, 1e-8);
    EXPECT_NEAR(comp.second, 0, 1e-8);

    little = std::make_pair(0, 3);
    comp   = complement(little, big);
    EXPECT_NEAR(comp.first, -1, 1e-8);
    EXPECT_NEAR(comp.second, 0, 1e-8);

    little = std::make_pair(-1, 0);
    comp   = complement(little, big);
    EXPECT_NEAR(comp.first, 0, 1e-8);
    EXPECT_NEAR(comp.second, 2, 1e-8);

    little = std::make_pair(-2, 0);
    comp   = complement(little, big);
    EXPECT_NEAR(comp.first, 0, 1e-8);
    EXPECT_NEAR(comp.second, 2, 1e-8);

    little = std::make_pair(-2, 3);
    comp   = complement(little, big);
    EXPECT_TRUE(empty(comp));

    little = std::make_pair(1, -1);
    comp   = complement(little, big);
    EXPECT_NEAR(comp.first, -1, 1e-8);
    EXPECT_NEAR(comp.second, 2, 1e-8);

    little = std::make_pair(-1.1, -1);
    comp   = complement(little, big);
    EXPECT_NEAR(comp.first, -1, 1e-8);
    EXPECT_NEAR(comp.second, 2, 1e-8);
}
