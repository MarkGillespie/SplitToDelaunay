#include "Interval.h"

bool empty(Interval i) {
    return (i.first >= i.second) || std::isnan(i.first) || std::isnan(i.second);
}

double distTo(Interval i, double d) {
    my_assert(!empty(i), "Cannot find distance to empty interval");
    if (d >= i.first && d <= i.second) {
        return 0;
    } else if (d < i.first) {
        return i.first - d;
    } else {
        return d - i.second;
    }
}

Interval complement(Interval little, Interval big) {
    if (empty(little)) {
        return big;
    } else {
        double leftDist  = little.first - big.first;
        double rightDist = big.second - little.second;

        if (leftDist > 1e-4 && rightDist > 1e-4) {
            cout << "little: " << little.first << ", " << little.second << endl;
            cout << "big   : " << big.first << ", " << big.second << endl;
            my_assert(false,
                      "little should have one endpoint at big's endpoint");
        }

        if (leftDist > rightDist) {
            return std::make_pair(big.first, fmin(little.first, big.second));
        } else {
            return std::make_pair(fmax(little.second, big.first), big.second);
        }
    }
}

Interval intersect(Interval a, Interval b) {
    return std::make_pair(fmax(a.first, b.first), fmin(a.second, b.second));
}

// Takes a set of indices and a set of intervals. Computes the intersection of
// the intervals with given indices. E.g
//         intersectIndices({0, 2, 4}, {I0, I1, I2, I3, I4}) = I0 ∩ I2 ∩ I4
Interval intersectIndices(const std::vector<size_t>& inds,
                          const std::array<Interval, 5>& I) {

    if (inds.empty()) {
        return std::make_pair(1.0, -1.0);
    }
    Interval out = I[0];
    for (size_t iI = 0; iI < inds.size(); ++iI) {
        out = intersect(out, I[inds[iI]]);
    }
    return out;
}

// Computes the intersection of the disk with the x-axis. Represents the
// resulting interval with its lower end first. Returns (1, -1) if the
// intersection is empty
Interval diskInterval(Disk disk) {
    double height = disk.first.y;
    double r      = disk.second;
    if (abs(height) > r) {
        return std::make_pair(1, -2);
    }
    double theta = asin(-height / r);
    return std::minmax(disk.first.x - r * cos(theta),
                       disk.first.x + r * cos(theta));
}

Disk circumcircle(Vector2 v1, Vector2 v2, Vector2 v3) {
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
