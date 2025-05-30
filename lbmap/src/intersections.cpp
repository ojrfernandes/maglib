#include "intersections.h"

// Curve class constructor
curve::curve(const std::string &filename, const size_t numPoints) {
    // Load curve from file
    if (!loadFromFile(filename, numPoints)) {
        std::cerr << "Error: Could not load curve from file " << filename << std::endl;
    }
}

// Bounding box class constructor
boundingBox::boundingBox(segment s) {
    // Find the minimum and maximum values of R and Z
    minR = std::min(s.p1.R, s.p2.R);
    maxR = std::max(s.p1.R, s.p2.R);
    minZ = std::min(s.p1.Z, s.p2.Z);
    maxZ = std::max(s.p1.Z, s.p2.Z);
}

// Check if a point is contained within the bounding box
bool boundingBox::contains(const point &p) const {
    return (p.R >= minR && p.R <= maxR && p.Z >= minZ && p.Z <= maxZ);
}

// Calculate the Euclidean distance between two points
double point::distanceTo(const point &other) const {
    return std::sqrt(std::pow(R - other.R, 2) + std::pow(Z - other.Z, 2));
}

// Calculate the angle between two points with respect to a given origin point
double point::angleTo(const point &other, const point &origin) const {
    // Calculate the vectors from the origin to the two points
    double x1 = this->R - origin.R;
    double y1 = this->Z - origin.Z;
    double x2 = other.R - origin.R;
    double y2 = other.Z - origin.Z;

    // Calculate the dot product
    double dotProduct = x1 * x2 + y1 * y2;

    // Calculate the magnitudes of the vectors
    double mag1 = std::sqrt(x1 * x1 + y1 * y1);
    double mag2 = std::sqrt(x2 * x2 + y2 * y2);

    // Calculate the cosine of the angle
    double cosTheta = dotProduct / (mag1 * mag2);

    // Clamp the cosine value to the range [-1, 1] to avoid NaN from acos
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

    // Calculate the angle in radians
    double angle = std::acos(cosTheta);

    return angle;
}

// Helper function to find the orientation of the ordered triplet (p, q, r)
// Returns 0 if p, q and r are collinear, 1 if clockwise, 2 if counterclockwise.
int segment::orientation(const point &p, const point &q, const point &r) {
    double val = (q.Z - p.Z) * (r.R - q.R) - (q.R - p.R) * (r.Z - q.Z);
    if (val == 0)
        return 0;             // Collinear
    return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
}

// Checks if point q lies on the segment pr
bool segment::onSegment(const point &p, const point &q, const point &r) {
    if (q.R <= std::max(p.R, r.R) && q.R >= std::min(p.R, r.R) &&
        q.Z <= std::max(p.Z, r.Z) && q.Z >= std::min(p.Z, r.Z)) {
        return true;
    }
    return false;
}

// Checks if segments s1 and s2 intersect
bool segment::doIntersect(const segment &otherSegment) {
    point p1 = this->p1;
    point q1 = this->p2;
    point p2 = otherSegment.p1;
    point q2 = otherSegment.p2;

    // Find the four orientations needed for the general and special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1))
        return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1))
        return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2))
        return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2))
        return true;

    // No intersection
    return false;
}

// Find the intersection point of two segments
point segment::findIntersectionWith(const segment &otherSegment) const {
    segment s1 = *this;
    segment s2 = otherSegment;

    // Coefficients of the line equations
    double A1 = s1.p2.Z - s1.p1.Z;
    double B1 = s1.p1.R - s1.p2.R;
    double C1 = A1 * s1.p1.R + B1 * s1.p1.Z;

    // Coefficients of the line equations
    double A2 = s2.p2.Z - s2.p1.Z;
    double B2 = s2.p1.R - s2.p2.R;
    double C2 = A2 * s2.p1.R + B2 * s2.p1.Z;

    // Calculate the determinant
    double det = A1 * B2 - A2 * B1;
    if (det == 0) {
        throw std::runtime_error("No intersection: segments are parallel or collinear.");
    }

    // Find the intersection point
    double R = (B2 * C1 - B1 * C2) / det;
    double Z = (A1 * C2 - A2 * C1) / det;
    return point(R, Z);
}

// Check if the bounding boxes of two segments overlap
bool boundingBox::overlaps(const boundingBox &other) const {
    return (maxR >= other.minR && minR <= other.maxR && maxZ >= other.minZ && minZ <= other.maxZ);
}

// Function to read points from a .dat file into a vector of point structs
bool curve::loadFromFile(const std::string &filename, const size_t numPoints) {
    // Open file
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // Get the number of lines in the file
    if (numPoints == 0) {
        std::string line;
        int numLines = 0;
        while (std::getline(file, line)) {
            numLines++;
        }
        // Reserve space in the vector for all data points
        this->curvePoints.reserve(numLines);

        // Reset file stream position to the beginning
        file.clear(); // Clear EOF flag
        file.seekg(0, std::ios::beg);
    } else {
        this->curvePoints.reserve(numPoints);
    }

    // Write data to the vector
    double R, Z;
    size_t count = 0;
    while (file >> R >> Z) {
        this->curvePoints.emplace_back(R, Z);
        count++;
        // If numPoints is set and we've read enough points, stop reading
        if (numPoints != 0 && count >= numPoints) {
            break;
        }
    }

    return true;
}

// Method to find the intersection points with another curve
std::vector<point> curve::intersectionsWith(curve &otherCurve) {
    std::vector<point> intersections;

    // Loop over consecutive pairs of points from the two sets
    for (size_t i = 0; i < this->curvePoints.size() - 1; ++i) {
        // for (size_t i = this->curvePoints.size() - 1; i > 0; --i) {
        segment thisSegment(this->curvePoints[i], this->curvePoints[i + 1]);

        for (size_t j = 0; j < otherCurve.curvePoints.size() - 1; ++j) {
            // for (size_t j = otherCurve.curvePoints.size() - 1; j > 0; --j) {
            segment otherSegment(otherCurve.curvePoints[j], otherCurve.curvePoints[j + 1]);

            // Check if the bounding boxes of the two segments overlap
            if (!boundingBox(thisSegment).overlaps(boundingBox(otherSegment))) {
                continue; // Skip the intersection check
            } else if (thisSegment.doIntersect(otherSegment)) {
                // Find the actual intersection point (not done yet)
                point intersection = thisSegment.findIntersectionWith(otherSegment);
                intersections.push_back(intersection);
            }
        }
    }
    return intersections;
}