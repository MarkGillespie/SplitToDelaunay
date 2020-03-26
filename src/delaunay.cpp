#include "delaunay.h"

std::vector<std::vector<size_t>> Splitter::intervalCombinations =
    Splitter::genIntervalCombinations();

std::vector<std::vector<size_t>> Splitter::genIntervalCombinations() {
    std::vector<std::vector<size_t>> combos;
    combos.reserve(14);
    for (size_t iSkip = 1; iSkip < 5; ++iSkip) {
        std::vector<size_t> skippedOne;
        for (size_t i = 1; i < 5; ++i) {
            if (i != iSkip) {
                skippedOne.push_back(i);
            }
        }
        combos.push_back(skippedOne);
    }
    for (size_t iInclude = 1; iInclude < 5; ++iInclude) {
        for (size_t jInclude = 1; jInclude < iInclude; ++jInclude) {
            combos.push_back(std::vector<size_t>{iInclude, jInclude});
        }
    }

    for (size_t i = 1; i < 5; ++i) {
        combos.push_back(std::vector<size_t>{i});
    }

    return combos;
}

Splitter::Splitter(VertexPositionGeometry& geo_) : geo(geo_) {
    geo.requireEdgeLengths();
    geo.requireCornerAngles();
    double lMin     = geo.edgeLengths[geo.mesh.edge(0)];
    double thetaMin = M_PI;

    for (Edge e : geo.mesh.edges()) {
        lMin = fmin(lMin, geo.edgeLengths[e]);
    }

    for (Corner c : geo.mesh.corners()) {
        thetaMin = fmin(thetaMin, geo.cornerAngles[c]);
    }

    double s = sin(thetaMin);
    rho_v    = fmin((lMin * s) / (0.5 + s), lMin / 2);
    rho_e    = 2 * rho_v * s;

    flatEdge = EdgeData<char>(geo.mesh, false);
}

void Splitter::flipFlatEdgesToDelaunay(std::array<Edge, 6> edges) {
    std::deque<Edge> edgesToCheck;
    EdgeData<char> inQueue(geo.mesh, false);
    for (Edge e : edges) {
        if (flatEdge[e]) {
            edgesToCheck.push_back(e);
            inQueue[e] = true;
        }
    }

    while (!edgesToCheck.empty()) {

        // Get the top element from the queue of possibily non-Delaunay edges
        Edge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e] = false;

        if (!isDelaunay(e)) {
            bool wasFlipped = geo.mesh.flip(e);

            if (!wasFlipped) continue;

            // Handle the aftermath of a flip

            // Add neighbors to queue, as they may need flipping now
            Halfedge he   = e.halfedge();
            Halfedge heN  = he.next();
            Halfedge heT  = he.twin();
            Halfedge heTN = heT.next();
            std::array<Edge, 4> neighEdges{heN.edge(), heN.next().edge(),
                                           heTN.edge(), heTN.next().edge()};
            for (Edge nE : neighEdges) {
                if (flatEdge[nE] && !inQueue[nE]) {
                    edgesToCheck.push_back(nE);
                    inQueue[nE] = true;
                }
            }
        }
    }
    geo.refreshQuantities();
}

void Splitter::splitGeometry(bool verbose) {
    Halfedge newHe1, newHe2;

    std::vector<Vector3> samplePositions;
    edgeSamplePoints = EdgeData<std::vector<double>>(geo.mesh);
    for (Edge e : geo.mesh.edges()) {
        edgeSamplePoints[e] = edgePoints(e);
        geo.requireVertexPositions();
        Vector3 p = geo.vertexPositions[e.halfedge().vertex()];
        Vector3 q = geo.vertexPositions[e.halfedge().twin().vertex()];
        double l  = (p - q).norm();
        for (double d : edgeSamplePoints[e]) {
            samplePositions.push_back(p * (d / l) + q * (1 - d / l));
        }
    }

    EdgeData<size_t> eIdx = geo.mesh.getEdgeIndices();

    std::deque<DynamicEdge> edgesToCheck;
    EdgeData<char> inQueue(geo.mesh, false);
    for (Edge e : geo.mesh.edges()) {
        if (flatEdge[e]) {
            edgesToCheck.push_back(DynamicEdge(e));
            inQueue[e] = true;
        }
    }

    size_t nSplits = 0;
    bool done      = false;
    while (!edgesToCheck.empty()) {
        // Get the top element from the queue of possibily non-Delaunay edges
        DynamicEdge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e.decay()] = false;

        if (!edgeSamplePoints[e].empty() && !isDelaunay(e)) {
            done = false;

            std::tie(newHe1, newHe2) = splitEdge(e);

            // Handle the aftermath of a flip
            nSplits++;

            std::array<Edge, 6> neighEdges{newHe1.next().edge(),
                                           newHe1.next().next().edge(),
                                           newHe2.next().edge(),
                                           newHe2.twin().next().edge(),
                                           newHe2.twin().next().next().edge(),
                                           newHe1.twin().next().edge()};

            flipFlatEdgesToDelaunay(neighEdges);
            for (Edge nE : neighEdges) {
                if (!flatEdge[nE] && !inQueue[nE]) {
                    edgesToCheck.emplace_back(nE);
                    inQueue[nE] = true;
                }
            }

            eIdx = geo.mesh.getEdgeIndices();
        }
    }


    geo.refreshQuantities();

    double isOverallDelaunay = true;
    EdgeData<double> edgeIsDelaunay(geo.mesh);
    bool printed = false;
    for (Edge e : geo.mesh.edges()) {
        isOverallDelaunay = isOverallDelaunay && isDelaunay(e);
        if (!printed && !isDelaunay(e)) {
            cout << "Edge " << eIdx[e] << " is not Delaunay" << endl;
            cout << "It is " << (flatEdge[e] ? "" : "not ") << "flat" << endl;
            cout << "It has " << (edgeSamplePoints[e].empty() ? "no " : "some ")
                 << "sample points" << endl;
            printed = true;
        }
        edgeIsDelaunay[e] = isDelaunay(e) ? 1 : -1;
    }

    if (verbose) {
        cout << "Finished splitting to delaunay. Used " << nSplits
             << " splits. Mesh is" << (isOverallDelaunay ? "" : " not")
             << " Delaunay" << endl;
    }
}

// Positions of "Delaunay sampling points" along this edge
std::vector<double> Splitter::edgePoints(Edge e) {
    geo.requireEdgeLengths();
    double edgeLen = geo.edgeLengths[e];

    if (edgeLen < 2 * rho_v + 1e-12) {
        return std::vector<double>{edgeLen / 2.0};
    } else {
        double middleLen        = edgeLen - 2 * rho_v;
        size_t nMiddleIntervals = ceil(middleLen / rho_e);
        double middleSpacing    = middleLen / nMiddleIntervals;
        std::vector<double> samplePoints{rho_v};
        for (size_t i = 1; i <= nMiddleIntervals; ++i) {
            samplePoints.push_back(rho_v + i * middleSpacing);
        }

        double farPtErr =
            abs((edgeLen - samplePoints[samplePoints.size() - 1]) - rho_v);

        my_assert(abs((edgeLen - samplePoints[samplePoints.size() - 1]) -
                      rho_v) < 1e-8,
                  "Error, far point in the wrong place");

        return samplePoints;
    }
}

// TODO: binary search instead
size_t Splitter::closestToInterval(Interval I, const std::vector<double>& pts) {
    double minDist      = I.second - I.first;
    size_t minDistIndex = 0;
    for (size_t iP = 0; iP < pts.size(); ++iP) {
        double pt = pts[iP];
        double dist;
        if (pt >= I.first && pt <= I.second) {
            return iP;
        } else if (pt < I.first) {
            dist = I.first - pt;
        } else { // pt > I.second
            dist = pt - I.second;
        }
        if (dist < minDist) {
            minDist      = dist;
            minDistIndex = iP;
        }
    }

    my_assert(pts.size() > 0, "How can this happen?");
    my_assert(minDistIndex >= 0 && minDistIndex < pts.size(),
              "list index out of range");
    // cout << "\t not in interval, just close" << endl;
    return minDistIndex;
}

std::pair<Halfedge, Halfedge> Splitter::splitEdge(Edge e) {
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
    std::array<Vector2, 8> v         = layOutButterfly(e);
    std::vector<double> samplePoints = edgeSamplePoints[e];

    Interval cInt     = diskInterval(circumcenter(v[0], v[4], v[1]));
    Interval dInt     = diskInterval(circumcenter(v[5], v[2], v[1]));
    Interval eInt     = diskInterval(circumcenter(v[3], v[2], v[6]));
    Interval fInt     = diskInterval(circumcenter(v[3], v[7], v[0]));
    Interval leftInt  = diskInterval(circumcenter(v[1], v[2], v[3]));
    Interval rightInt = diskInterval(circumcenter(v[0], v[1], v[3]));

    auto complement = [&](Interval I) {
        if (empty(I)) {
            return std::make_pair(0.0, v[0].x);
        } else if (I.first < 1e-8) {
            return std::make_pair(I.second, v[0].x);
        } else if (I.second > v[0].x - 1e-8) {
            return std::make_pair(0.0, I.first);
        } else {
            my_assert(false,
                      "Interval should have one endpoint at edge endpoint");
        }
    };

    // I1 ... I4 are the complements of the circumcircle intervals
    std::array<Interval, 5> I{intersect(leftInt, rightInt), complement(cInt),
                              complement(dInt), complement(eInt),
                              complement(fInt)};
    auto intersectAll = [&](std::vector<size_t> inds) {
        if (inds.empty()) {
            return std::make_pair(1.0, -1.0);
        }
        Interval out = I[0];
        for (size_t iI = 0; iI < inds.size(); ++iI) {
            out = intersect(out, I[inds[iI]]);
        }
        return out;
    };

    Interval bestIntersector = intersectAll(std::vector<size_t>{0, 1, 2, 3, 4});
    bool allIntersect        = !empty(bestIntersector);
    size_t comboIndex        = 0;
    while (empty(bestIntersector)) {
        bestIntersector = intersectAll(intervalCombinations[comboIndex++]);
    }

    geo.requireEdgeLengths();
    size_t splitPointIdx = closestToInterval(bestIntersector, samplePoints);
    double splitDist     = samplePoints[splitPointIdx];
    double splitBary     = 1 - splitDist / geo.edgeLengths[e];

    my_assert(splitBary > 0 && splitBary < 1, "splitBary is invalid");

    std::vector<double> firstHalfSamplePoints, secondHalfSamplePoints;
    for (size_t iS = 0; iS < splitPointIdx; ++iS) {
        firstHalfSamplePoints.push_back(samplePoints[iS]);
    }
    for (size_t iS = splitPointIdx + 1; iS < samplePoints.size(); ++iS) {
        secondHalfSamplePoints.push_back(samplePoints[iS] - splitDist);
    }

    // Compute new vertex position before splitting e
    geo.requireVertexPositions();
    Vector3 p2     = geo.vertexPositions[e.halfedge().vertex()];
    Vector3 p0     = geo.vertexPositions[e.halfedge().twin().vertex()];
    Vector3 newPos = splitBary * p2 + (1 - splitBary) * p0;

    // ------------ e -------------->
    // -- newHe1 --> . --- newHe2 -->
    Halfedge newHe2                 = geo.mesh.splitEdgeTriangular(e);
    Halfedge newHe1                 = newHe2.next().next().twin().next().next();
    edgeSamplePoints[newHe1.edge()] = firstHalfSamplePoints;
    edgeSamplePoints[newHe2.edge()] = secondHalfSamplePoints;

    flatEdge[newHe1.edge()] = false;
    flatEdge[newHe2.edge()] = false;

    flatEdge[newHe1.next().edge()]            = true;
    flatEdge[newHe2.twin().next().edge()]     = true;
    geo.inputVertexPositions[newHe2.vertex()] = newPos;

    // Only need a few angles here
    geo.refreshQuantities();

    return std::make_pair(newHe1, newHe2);
}

bool Splitter::isDelaunay(Edge e) { return geo.edgeCotanWeight(e) > -1e-12; }

bool Splitter::empty(Interval i) { return i.first >= i.second; }
Interval Splitter::intersect(Interval a, Interval b) {
    return std::make_pair(fmax(a.first, b.first), fmin(a.second, b.second));
}

// Computes the intersection of the disk with the x-axis. Represents the
// resulting interval with its lower end first. Returns (1, -1) if the
// intersection is empty
Interval Splitter::diskInterval(Disk disk) {
    double height = disk.first.y;
    double r      = disk.second;
    if (abs(height) > r) {
        return std::make_pair(1, -2);
    }
    double theta = asin(-height / r);
    return std::minmax(disk.first.x - r * cos(theta),
                       disk.first.x + r * cos(theta));
}

Disk Splitter::circumcenter(Vector2 v1, Vector2 v2, Vector2 v3) {
    double a                = (v2 - v3).norm();
    double b                = (v3 - v1).norm();
    double c                = (v1 - v2).norm();
    double a2               = a * a;
    double b2               = b * b;
    double c2               = c * c;
    Vector3 circumcenterLoc = {a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2),
                               c2 * (a2 + b2 - c2)};
    double s                = sum(circumcenterLoc);
    circumcenterLoc /= s; // normalize barycentric coordinates

    Vector2 circumcenter = circumcenterLoc.x * v1 + circumcenterLoc.y * v2 +
                           circumcenterLoc.z * v3;

    double circumradius = (v1 - circumcenter).norm();
    return std::make_pair(circumcenter, circumradius);
}

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
 *  I refer to edges by the adjacent faces. E.g. e is edge AB
 *     For edges on the boundary, I never use the top and
 *     bottom edges, so I just refer to them by face.
 *     E.g. edge Ddiag is the left edge of face D.
 *  I refer to corners by the face and vertex. E.g. corner A2
 *     is the corner at the base of edge e in face A
 */
std::array<Vector2, 8> Splitter::layOutButterfly(Edge e) {
    std::array<Vector2, 8> butterfly;
    geo.requireEdgeLengths();
    geo.requireCornerAngles();
    const EdgeData<double>& len   = geo.edgeLengths;
    const CornerData<double>& ang = geo.cornerAngles;

    double lAB     = len[e];
    double lAD     = len[e.halfedge().next().next().edge()];
    double lBE     = len[e.halfedge().twin().next().edge()];
    double thetaA0 = ang[e.halfedge().next().corner()];
    double thetaA2 = ang[e.halfedge().corner()];
    double thetaB0 = ang[e.halfedge().twin().corner()];
    double thetaB2 = ang[e.halfedge().twin().next().corner()];

    butterfly[0] = Vector2{lAB, 0};
    butterfly[1] = Vector2{cos(thetaA2), sin(thetaA2)} * lAD;
    butterfly[2] = Vector2{0, 0};
    butterfly[3] = Vector2{cos(thetaB2), -sin(thetaB2)} * lBE;

    double lCdiag = len[e.halfedge().next().twin().next().edge()];
    double lDdiag = len[e.halfedge().next().next().twin().next().next().edge()];
    double lEdiag = len[e.halfedge().twin().next().twin().next().edge()];
    double lFdiag =
        len[e.halfedge().twin().next().next().twin().next().next().edge()];

    double thetaC0 = ang[e.halfedge().next().twin().next().corner()];
    double thetaD2 = ang[e.halfedge().next().next().twin().corner()];
    double thetaE2 = ang[e.halfedge().twin().next().twin().next().corner()];
    double thetaF0 = ang[e.halfedge().twin().next().next().twin().corner()];

    butterfly[4] = butterfly[0] + Vector2{cos(M_PI - thetaC0 - thetaA0),
                                          sin(M_PI - thetaC0 - thetaA0)} *
                                      lCdiag;
    butterfly[5] =
        Vector2{cos(thetaA2 + thetaD2), sin(thetaA2 + thetaD2)} * lDdiag;

    butterfly[6] =
        Vector2{cos(thetaB2 + thetaE2), -sin(thetaB2 + thetaE2)} * lEdiag;

    butterfly[7] = butterfly[0] + Vector2{cos(M_PI - thetaF0 - thetaB0),
                                          -sin(M_PI - thetaF0 - thetaB0)} *
                                      lFdiag;

    return butterfly;
}
