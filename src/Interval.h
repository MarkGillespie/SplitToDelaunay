#pragma once

#include "geometrycentral/utilities/vector2.h"

#include "my_assert.h"

#include <array>
#include <vector>

using std::cout;
using std::endl;
using namespace geometrycentral;

// Intervals must have lower end as first component.
// If interval.first > interval.second, the interval is treated as empty
using Interval = std::pair<double, double>;
using Disk     = std::pair<Vector2, double>;

bool empty(Interval i);
double distTo(Interval i, double d);

// Staps little to big's closes endpoint and then returns the interval big -
// little
Interval complement(Interval little, Interval big);

Interval intersect(Interval a, Interval b);
Interval intersectIndices(const std::vector<size_t>& inds,
                          const std::array<Interval, 5>& I);

// Computes the intersection of the disk with the x-axis. Represents the
// resulting interval with its lower end first. Returns (1, -1) if the
// intersection is empty
Interval diskInterval(Disk disk);

Disk circumcircle(Vector2 v1, Vector2 v2, Vector2 v3);
