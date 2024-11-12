#ifndef LOBES_H
#define LOBES_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Point and Segment structs to represent points and segments in 2D space
class point {
  public:
    double R, Z;

    // Class constructor
    point(double r, double z) : R(r), Z(z) {}

    // Method to calculate the Euclidean distance between two points
    double distanceTo(const point &other) const;
};

class segment {
  public:
    point p1, p2;

    // Class constructor
    segment(const point &point1, const point &point2) : p1(point1), p2(point2) {}

    bool doIntersect(const segment &otherSegment);
    point findIntersectionWith(const segment &otherSegment) const;

  private:
    int orientation(const point &p, const point &q, const point &r);
    bool onSegment(const point &p, const point &q, const point &r);
};

class boundingBox {
  public:
    double minR, maxR, minZ, maxZ;

    // Class constructor
    boundingBox(segment s);

    bool overlaps(const boundingBox &other) const;
};

class curve {
  public:
    std::vector<point> curvePoints;

    // Class constructor
    curve(const std::string &filename);

    // Method to find the intersection points with another curve
    std::vector<point> intersectionsWith(const curve &other) const;

  private:
    // Method to read points from a file
    bool loadFromFile(const std::string &filename, const size_t numPoints = 0);
};

#endif