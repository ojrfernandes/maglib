#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H
#define INTERSECTIONS_V 241114 // version (yy.mm.dd)

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Point and Segment structs to represent points and segments in 2D space
class point {
  public:
    // Coordinates of the point
    double R, Z;

    // operator ==
    bool operator==(const point &other) const {
        return (std::abs(R - other.R) < tolerance &&
                std::abs(Z - other.Z) < tolerance);
    }
    // Default constructor no arguments
    point() = default;
    // Class constructor
    point(double r, double z) : R(r), Z(z) {}
    // Method to calculate the Euclidean distance between two points
    double distanceTo(const point &other) const;

  private:
    // Customizable tolerance for floating-point comparison
    static constexpr double tolerance = 1e-14;
};

class segment {
  public:
    // Endpoints of the segment
    point p1, p2;

    // Default constructor
    segment() = default;
    // Class constructor
    segment(const point &point1, const point &point2) : p1(point1), p2(point2) {}
    // Checks if segments s1 and s2 intersect
    bool doIntersect(const segment &otherSegment);
    // Find the intersection point of two segments
    point findIntersectionWith(const segment &otherSegment) const;
    // Helper function to find the orientation of the ordered triplet (p, q, r)
    // Returns 0 if p, q and r are collinear, 1 if clockwise, 2 if counterclockwise.
    int orientation(const point &p, const point &q, const point &r);
    // Checks if point q lies on the segment pr
    bool onSegment(const point &p, const point &q, const point &r);
};

class boundingBox {
  public:
    // Minimum and maximum R and Z coordinates of the bounding box
    double minR, maxR, minZ, maxZ;

    // Default constructor
    boundingBox() = default;
    // Class constructor
    boundingBox(segment s);
    // Check if the bounding boxes of two segments overlap
    bool overlaps(const boundingBox &other) const;
    // Check if a point is contained within the bounding box
    bool contains(const point &p) const;
};

class curve {
  public:
    // Vector of points representing the curve
    std::vector<point> curvePoints;

    // Default constructor
    curve() = default;
    // Class constructor load from file
    curve(const std::string &filename);
    // Method to find the intersection points with another curve
    std::vector<point> intersectionsWith(curve &other);

  private:
    // Function to read points from a .dat file into a vector of point structs
    bool loadFromFile(const std::string &filename, const size_t numPoints = 0);
};

#endif