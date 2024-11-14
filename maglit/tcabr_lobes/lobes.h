#ifndef LOBES_H
#define LOBES_H
#define LOBES_V 241114 // version (yy.mm.dd)

#include "intersections.h"

class lobe {
  public:
    // Default constructor
    lobe() = default;
    // class constructor with boundary points
    lobe(const point &p1, const point &p2);

    // intersection points baoundin the lobe
    point pBoundary_1, pBoundary_2;
    // center coordinate of the lobe
    point center;
    // equilibrium and perturbed curves bounding the lobe
    curve cBoundary_eq, cBoundary_ptb, curveBoundary;
    // area and perimeter of the lobe
    double area, perimeter;

    // get boundary curve of the lobe between two points for equilibrium and perturbed sets
    void getBoundaries(curve &equilibrium, curve &perturbed);
    // get perimeter of the lobe
    void getPerimeter();
    // get area of the lobe
    void getArea();
};

#endif