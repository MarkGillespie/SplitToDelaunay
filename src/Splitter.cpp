#include "Splitter.h"

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

    // Set default parent index to -1 so that newly created edges have parent
    // index -1
    parentIndex           = EdgeData<int>(geo.mesh, -1);
    EdgeData<size_t> eIdx = geo.mesh.getEdgeIndices();
    for (Edge e : geo.mesh.edges()) {
        parentIndex[e] = eIdx[e];
    }
}

void Splitter::flipFlatEdgesToDelaunay(std::set<Edge> edges) {
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
}

EdgeData<int> Splitter::splitGeometry(bool verbose) {
    Halfedge newHe1, newHe2;

    HalfedgeData<int> originalFaceIndices(geo.mesh);
    std::vector<std::set<Edge>> edgesPerFace(geo.mesh.nFaces(),
                                             std::set<Edge>{});

    FaceData<size_t> fIdx = geo.mesh.getFaceIndices();
    for (Halfedge he : geo.mesh.halfedges()) {
        if (he.isInterior()) {
            originalFaceIndices[he] = fIdx[he.face()];
            edgesPerFace[originalFaceIndices[he]].insert(he.edge());
        } else {
            originalFaceIndices[he] = -1;
        }
    }

    std::deque<Edge> edgesToCheck;
    EdgeData<char> inQueue(geo.mesh, false);
    for (Edge e : geo.mesh.edges()) {
        edgesToCheck.push_back(e);
        inQueue[e] = true;
    }

    size_t nSplits = 0;
    while (!edgesToCheck.empty()) {

        // Get the top element from the queue of possibily non-Delaunay edges
        Edge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e] = false;

        // TODO: This doesn't quite do the length constraints properly
        if (isDelaunay(e) || geo.edgeLength(e) < 2 * rho_e) {
            continue;
        }

        int faceIndex     = originalFaceIndices[e.halfedge()];
        int faceIndexTwin = originalFaceIndices[e.halfedge().twin()];

        // Split the edge
        std::tie(newHe1, newHe2) = splitEdge(e);

        // Handle the aftermath of a split
        nSplits++;

        // Update parent face pointers
        originalFaceIndices[newHe1]        = faceIndex;
        originalFaceIndices[newHe2]        = faceIndex;
        originalFaceIndices[newHe1.twin()] = faceIndexTwin;
        originalFaceIndices[newHe2.twin()] = faceIndexTwin;

        // Add new edges to the appropriate face lists, and
        // identify the edges that lie in the faces of the original mesh
        // which are incident on edge e
        // These are the edges that we might need to flip
        std::set<Edge> faceEdges;
        if (faceIndex >= 0) {
            edgesPerFace[faceIndex].insert(newHe1.edge());
            edgesPerFace[faceIndex].insert(newHe2.edge());
            edgesPerFace[faceIndex].insert(newHe1.next().edge());

            faceEdges.insert(edgesPerFace[faceIndex].begin(),
                             edgesPerFace[faceIndex].end());
        }
        if (faceIndexTwin >= 0) {
            edgesPerFace[faceIndexTwin].insert(newHe1.edge());
            edgesPerFace[faceIndexTwin].insert(newHe2.edge());
            edgesPerFace[faceIndexTwin].insert(newHe2.twin().next().edge());

            faceEdges.insert(edgesPerFace[faceIndexTwin].begin(),
                             edgesPerFace[faceIndexTwin].end());
        }

        // Flip any flippable edges to Delaunay and push non-Delaunay
        // non-flippable edges onto the queue
        flipFlatEdgesToDelaunay(faceEdges);
        for (Edge nE : faceEdges) {
            if (!flatEdge[nE] && !inQueue[nE] && !isDelaunay(nE)) {
                edgesToCheck.push_back(nE);
                inQueue[nE] = true;
            }
        }
    }

    geo.refreshQuantities();
    if (verbose) {
        cout << "Finished splitting to delaunay. Used " << nSplits << " splits."
             << endl;
        checkMesh();
    }

    return parentIndex;
}

void Splitter::checkMesh() {
    EdgeData<size_t> eIdx = geo.mesh.getEdgeIndices();

    double isOverallDelaunay = true;
    bool printedFlat         = false;
    bool printedEssential    = false;

    size_t nFlatBad          = 0;
    size_t nFlatBadFlippable = 0;
    for (Edge e : geo.mesh.edges()) {
        isOverallDelaunay = isOverallDelaunay && isDelaunay(e);
        if (!printedFlat && flatEdge[e] && !isDelaunay(e)) {
            cout << "Edge " << eIdx[e] << " is not Delaunay" << endl;
            cout << "It is flat" << endl;
            bool flippable = geo.mesh.flip(e);
            geo.mesh.flip(e);
            cout << "It was " << (flippable ? "" : "not ") << "flippable"
                 << endl;
            printedFlat = true;
        }
        if (!printedEssential && !flatEdge[e] && !isDelaunay(e)) {
            cout << "Edge " << eIdx[e] << " is not Delaunay" << endl;
            cout << "It is not flat" << endl;
            cout << "It is " << (geo.edgeLength(e) >= 2 * rho_e ? "" : "not ")
                 << "long enough to split" << endl;
            printedEssential = true;
        }
        if (flatEdge[e] && !isDelaunay(e)) {
            nFlatBad += 1;
            if (geo.mesh.flip(e)) {
                nFlatBadFlippable += 1;
                geo.mesh.flip(e);
            }
        }
    }
    cout << "Mesh is" << (isOverallDelaunay ? "" : " not") << " Delaunay"
         << endl;

    if (!isOverallDelaunay) {
        cout << "# non-delaunay flat edges: " << nFlatBad << "\t of those, "
             << nFlatBadFlippable << " are flippable" << endl;
    }
}

std::pair<Halfedge, Halfedge> Splitter::splitEdge(Edge e) {

    Interval splitI;
    if (e.isBoundary()) {
        splitI = computeBoundarySplitInterval(e);
    } else {
        splitI = computeInteriorSplitInterval(e);
    }
    double len       = geo.edgeLength(e);
    double mid       = len / 2;
    double splitDist = 0;

    if (splitI.first < mid && splitI.second > mid) {
        splitDist = mid;
    } else if (splitI.first >= mid) {
        // This method of enforcing spacing isn't quite right, but it seems to
        // work
        splitDist = fmin(splitI.first, len - rho_e);
    } else if (splitI.second <= mid) {
        splitDist = fmax(rho_e, splitI.second);
    }
    double splitBary = 1 - splitDist / geo.edgeLength(e);

    // Compute new vertex position before splitting edge e
    Vector3 p2         = geo.inputVertexPositions[e.halfedge().vertex()];
    Vector3 p0         = geo.inputVertexPositions[e.halfedge().twin().vertex()];
    Vector3 newPos     = splitBary * p2 + (1 - splitBary) * p0;
    int oldParentIndex = parentIndex[e];

    // ------------ e -------------->
    // -- newHe1 --> . --- newHe2 -->
    Halfedge newHe2 = geo.mesh.splitEdgeTriangular(e);
    Halfedge newHe1 = newHe2.next().next().twin().next().next();

    flatEdge[newHe1.edge()]    = false;
    flatEdge[newHe2.edge()]    = false;
    parentIndex[newHe1.edge()] = oldParentIndex;
    parentIndex[newHe2.edge()] = oldParentIndex;

    flatEdge[newHe1.next().edge()]        = true;
    flatEdge[newHe2.twin().next().edge()] = true;

    geo.inputVertexPositions[newHe2.vertex()] = newPos;

    return std::make_pair(newHe1, newHe2);
}

// Finds the interval that you should make your edge split in
Interval Splitter::computeInteriorSplitInterval(Edge e) {
    std::array<Interval, 5> I = computeInteriorIntervals(e);

    Interval bestIntersector =
        intersectIndices(std::vector<size_t>{0, 1, 2, 3, 4}, I);

    // IntervalCombinations is a list of subsets of {0, ..., 4} ordered from
    // largest cardinality to smallest. All of the subsets in
    // intervalCombinations contain 0.
    // This lets us find the largest possible subset of {1, ..., 4} whose
    // intersection with interval 0 is nonempty
    // We keep trying subsets until we find one with a nonempty intersection,
    // and then return that one

    // TODO: What if there are multiple possible answers of the same size? Right
    // now, we pick an arbitrary one (well, the first lexicographically) but
    // maybe we could pick the best one in some way. For the boundary case with
    // fewer possibilities, I try to take an interval containing the edge's
    // midpoint if possible

    size_t comboIndex = 0;
    while (empty(bestIntersector)) {
        bestIntersector =
            intersectIndices(intervalCombinations[comboIndex++], I);
    }

    my_assert(!empty(bestIntersector), "At least I0 isn't empty");
    return bestIntersector;
}

Interval Splitter::computeBoundarySplitInterval(Edge e) {
    std::array<Interval, 3> I = computeBoundaryIntervals(e);

    Interval all = intersect(I[0], intersect(I[1], I[2]));
    if (!empty(all)) {
        return all;
    } else {
        Interval intersectC = intersect(I[0], I[1]);
        Interval intersectD = intersect(I[0], I[2]);

        if (empty(intersectC)) {
            if (empty(intersectD)) {
                return I[0];
            } else {
                return intersectD;
            }
        } else if (empty(intersectD)) {
            return intersectC;
        } else {

            // If I0 ∩ I1 and I0 ∩ I2 are nonempty, but I0 ∩ I1 ∩ I2 is empty,
            // then we could return either I0 ∩ I1 or I0 ∩ I2. We choose to
            // return the one which is closest to containing the midpoint of our
            // edge. If both contain the midpoint, we arbitrarily choose the
            // first one
            // My hope is that this leads to "more balanced" splits
            double len   = geo.edgeLength(e);
            double midPt = len / 2;
            double cDist = distTo(intersectC, midPt);
            double dDist = distTo(intersectD, midPt);
            if (cDist <= dDist) {
                return intersectC;
            } else {
                return intersectD;
            }
        }
    }
}

bool Splitter::isDelaunay(Edge e) {
    // All boundary edges are Delaunay. If the edge is not a boundary edge, then
    // we check if its cotan weight is positive, up to an arbitrary epsilon
    return e.isBoundary() || geo.edgeCotanWeight(e) > -1e-12;
}


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
 *  I refer to edges by the adjacent faces. E.g. e is edge AB
 *     For edges on the boundary, I never use the top and
 *     bottom edges, so I just refer to them by face.
 *     E.g. edge Ddiag is the left edge of face D.
 *  I refer to corners by the face and vertex. E.g. corner A2
 *     is the corner at the base of edge e in face A
 *
 * Compute the intervals corresponding to the circumdisks of
 * triangles in the neighborhood, but not incident on edge e (Triangles C,
 * D, E, F) in the butterfly
 * Some of these triangles might not exist of A or
 *  B are boundary faces
 */
std::array<Interval, 5> Splitter::computeInteriorIntervals(Edge e) {
    std::array<Interval, 5> intervals;
    std::array<Vector2, 8> v; // Vertex positions

    double lAB     = geo.edgeLength(e);
    double lAD     = geo.edgeLength(e.halfedge().next().next().edge());
    double lBE     = geo.edgeLength(e.halfedge().twin().next().edge());
    double thetaA0 = geo.cornerAngle(e.halfedge().next().corner());
    double thetaA2 = geo.cornerAngle(e.halfedge().corner());
    double thetaB0 = geo.cornerAngle(e.halfedge().twin().corner());
    double thetaB2 = geo.cornerAngle(e.halfedge().twin().next().corner());

    v[0] = Vector2{lAB, 0};
    v[1] = Vector2{cos(thetaA2), sin(thetaA2)} * lAD;
    v[2] = Vector2{0, 0};
    v[3] = Vector2{cos(thetaB2), -sin(thetaB2)} * lBE;

    Interval leftInt  = diskInterval(circumcircle(v[1], v[2], v[3]));
    Interval rightInt = diskInterval(circumcircle(v[0], v[1], v[3]));
    intervals[0]      = intersect(leftInt, rightInt);

    // intervals[0] is guaranteed not to be empty, but sometimes it's just a
    // single point. Here, we widen it slightly so that it doesn't fail tests
    // later
    if (empty(intervals[0])) {
        intervals[0] =
            std::make_pair(leftInt.second - 1e-6, rightInt.first + 1e-6);
    }

    Interval edgeInt = std::make_pair(v[2].x, v[0].x);
    if (e.halfedge().next().twin().isInterior()) {
        // Face C interior face
        double lCdiag =
            geo.edgeLength(e.halfedge().next().twin().next().edge());
        double thetaC0 =
            geo.cornerAngle(e.halfedge().next().twin().next().corner());
        v[4] = v[0] + Vector2{cos(M_PI - thetaC0 - thetaA0),
                              sin(M_PI - thetaC0 - thetaA0)} *
                          lCdiag;
        intervals[1] = diskInterval(circumcircle(v[0], v[4], v[1]));
    } else {
        intervals[1] = std::make_pair(1, -1);
    }
    if (e.halfedge().next().next().twin().isInterior()) {
        // Face D interior face
        double lDdiag = geo.edgeLength(
            e.halfedge().next().next().twin().next().next().edge());
        double thetaD2 =
            geo.cornerAngle(e.halfedge().next().next().twin().corner());
        v[5] = Vector2{cos(thetaA2 + thetaD2), sin(thetaA2 + thetaD2)} * lDdiag;
        intervals[2] = diskInterval(circumcircle(v[2], v[1], v[5]));
    } else {
        intervals[2] = std::make_pair(1, -1);
    }
    if (e.halfedge().twin().next().twin().isInterior()) {
        // Face E interior face
        double lEdiag =
            geo.edgeLength(e.halfedge().twin().next().twin().next().edge());
        double thetaE2 =
            geo.cornerAngle(e.halfedge().twin().next().twin().next().corner());
        v[6] =
            Vector2{cos(thetaB2 + thetaE2), -sin(thetaB2 + thetaE2)} * lEdiag;
        intervals[3] = diskInterval(circumcircle(v[6], v[3], v[2]));
    } else {
        intervals[3] = std::make_pair(1, -1);
    }
    if (e.halfedge().twin().next().next().twin().isInterior()) {
        // Face F interior face
        double lFdiag = geo.edgeLength(
            e.halfedge().twin().next().next().twin().next().next().edge());
        double thetaF0 =
            geo.cornerAngle(e.halfedge().twin().next().next().twin().corner());
        v[7] = v[0] + Vector2{cos(M_PI - thetaF0 - thetaB0),
                              -sin(M_PI - thetaF0 - thetaB0)} *
                          lFdiag;
        intervals[4] = diskInterval(circumcircle(v[7], v[0], v[3]));
    } else {
        intervals[4] = std::make_pair(1, -1);
    }


    // The intervals we computed so far were the intersections of circumdisks
    // with the edge
    // For all except the first interval, we now have to take their complement
    // in the edge

    Interval Iedge = std::make_pair(0, v[0].x);

    // Floating point might mess up interval endpoints. We fix them here
    // Note that it's possible that e.g. cInt has v0 as its left endpoint rather
    // than the right endpoint as I assume here That's okay, since I only care
    // about the intersection of cInt with the edge, which is empty in that case
    // anyway
    intervals[1].second = v[0].x;
    intervals[2].first  = v[2].x;
    intervals[3].first  = v[2].x;
    intervals[4].second = v[0].x;

    intervals[1] = complement(intervals[1], Iedge);
    intervals[2] = complement(intervals[2], Iedge);
    intervals[3] = complement(intervals[3], Iedge);
    intervals[4] = complement(intervals[4], Iedge);

    return intervals;
}

/*
 *
 *  v4_____________v1___________ v3
 *    \           / \           /
 *     \         /   \         /
 *      \   D   /     \   C   /
 *       \     /   A   \     /
 *        \   /         \   /
 *         \ /           \ /
 *        v2 ----- e --->-  v0
 *  I refer to edges by the adjacent faces. E.g. e is edge AB
 *     For edges on the boundary, I never use the top and
 *     bottom edges, so I just refer to them by face.
 *     E.g. edge Ddiag is the left edge of face D.
 *  I refer to corners by the face and vertex. E.g. corner A2
 *     is the corner at the base of edge e in face A
 */
std::array<Interval, 3> Splitter::computeBoundaryIntervals(Edge e) {
    std::array<Interval, 3> intervals;
    std::array<Vector2, 5> v; // Vertex positions

    double lAB     = geo.edgeLength(e);
    double lAD     = geo.edgeLength(e.halfedge().next().next().edge());
    double lBE     = geo.edgeLength(e.halfedge().twin().next().edge());
    double thetaA0 = geo.cornerAngle(e.halfedge().next().corner());
    double thetaA2 = geo.cornerAngle(e.halfedge().corner());

    v[0] = Vector2{lAB, 0};
    v[1] = Vector2{cos(thetaA2), sin(thetaA2)} * lAD;
    v[2] = Vector2{0, 0};

    // The left and right intervals work differently now
    // LeftInt goes from v2 to the point on e where a line perpendicular to
    // v2-v1 coming from v1 would intersect it
    Vector2 e21 = v[1] - v[2];
    Vector2 e21Perp{e21.y, -e21.x};
    e21Perp *= e21.y / -e21Perp.y;
    Vector2 intervalVec = e21 + e21Perp;
    my_assert(abs(intervalVec.y) < 1e-8,
              "You did some algebra wrong. This vector should lie along edge e "
              "(the x axis)");
    Interval leftInt = std::make_pair(0, intervalVec.x);

    Vector2 e01 = v[1] - v[2];
    Vector2 e01Perp{-e01.y, e01.x};
    e01Perp *= e01.y / -e01Perp.y;
    intervalVec = e01 + e01Perp;
    my_assert(abs(intervalVec.y) < 1e-8,
              "You did some algebra wrong. This vector should lie along edge e "
              "(the x axis)");
    Interval rightInt = std::make_pair(intervalVec.x, v[0].x);

    intervals[0] = intersect(leftInt, rightInt);

    // intervals[0] is guaranteed not to be empty, but sometimes it's just a
    // single point. Here, we widen it slightly so that it doesn't fail tests
    // later
    if (empty(intervals[0])) {
        intervals[0] =
            std::make_pair(leftInt.second - 1e-6, rightInt.first + 1e-6);
    }

    if (e.halfedge().next().twin().isInterior()) {
        // Face C is a real face
        double lCdiag =
            geo.edgeLength(e.halfedge().next().twin().next().edge());
        double thetaC0 =
            geo.cornerAngle(e.halfedge().next().twin().next().corner());

        v[3] = v[0] + Vector2{cos(M_PI - thetaC0 - thetaA0),
                              sin(M_PI - thetaC0 - thetaA0)} *
                          lCdiag;
        intervals[1] = diskInterval(circumcircle(v[0], v[3], v[1]));
    } else {
        intervals[1] = std::make_pair(1, -1);
    }


    if (e.halfedge().next().twin().isInterior()) {
        // Face D is a real face
        double lDdiag = geo.edgeLength(
            e.halfedge().next().next().twin().next().next().edge());
        double thetaD2 =
            geo.cornerAngle(e.halfedge().next().next().twin().corner());
        v[4] = Vector2{cos(thetaA2 + thetaD2), sin(thetaA2 + thetaD2)} * lDdiag;
        intervals[2] = diskInterval(circumcircle(v[2], v[1], v[4]));
    } else {
        intervals[2] = std::make_pair(1, -1);
    }

    // The intervals we computed so far were the intersections of circumdisks
    // with the edge
    // For all except the first interval, we now have to take their complement
    // in the edge

    Interval Iedge = std::make_pair(0, v[0].x);

    // Floating point might mess up interval endpoints. We fix them here
    // Note that it's possible that e.g. cInt has v0 as its left endpoint rather
    // than the right endpoint as I assume here That's okay, since I only care
    // about the intersection of cInt with the edge, which is empty in that case
    // anyway
    intervals[1].second = v[0].x;
    intervals[2].first  = v[2].x;

    intervals[1] = complement(intervals[1], Iedge);
    intervals[2] = complement(intervals[2], Iedge);

    return intervals;
}

std::vector<std::vector<size_t>> Splitter::intervalCombinations =
    Splitter::genIntervalCombinations();

std::vector<std::vector<size_t>> Splitter::genIntervalCombinations() {
    std::vector<std::vector<size_t>> combos;
    combos.reserve(14);
    for (size_t iSkip = 1; iSkip < 5; ++iSkip) {
        std::vector<size_t> skippedOne;
        for (size_t i = 0; i < 5; ++i) {
            if (i != iSkip) {
                skippedOne.push_back(i);
            }
        }
        combos.push_back(skippedOne);
    }
    for (size_t iInclude = 1; iInclude < 5; ++iInclude) {
        for (size_t jInclude = 1; jInclude < iInclude; ++jInclude) {
            combos.push_back(std::vector<size_t>{0, iInclude, jInclude});
        }
    }

    for (size_t i = 1; i < 5; ++i) {
        combos.push_back(std::vector<size_t>{0, i});
    }
    combos.push_back(std::vector<size_t>{0});

    return combos;
}
