#ifndef LOBES_H
#define LOBES_H
#define LOBES_V 250120 // version (yy.mm.dd)

#include "intersections.h"

class lobe {
  public:
    // Default constructor
    lobe() = default;
    // class constructor with boundary points
    lobe(const point &p1, const point &p2, curve &equilibrium, curve &perturbed, const point &magAxis, const point &xPoint);

    point pBoundary_1, pBoundary_2, midpoint, referencePoint; // intersection points bounding the lobe
    curve cBoundary_eq, cBoundary_ptb, curveBoundary;         // equilibrium and perturbed curves bounding the lobe
    double area, perimeter, hParameter, referenceAngle;       // area, perimeter, H-parameter and poloidal angle coordinate of the lobe

  private:
    // get the boundary curve of the lobe between two points for equilibrium and perturbed sets
    void getBoundaries(curve &equilibrium, curve &perturbed);
    // get the perimeter of the lobe
    void getPerimeter();
    // get the area of the lobe
    void getArea();
    // get the h parameter of the lobe
    void getHParameter();
    // get the midpoint of the lobe
    void getMidpoint();
    // get the angle between the lobe's reference pointy and the X-point with respect to the magnetic axis
    void getReferenceAngle();

    curve equilibrium, perturbed;
    point magAxis, xPoint;
};

#endif