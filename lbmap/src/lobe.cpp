#include "lobe.h"

// Class constructor
lobe::lobe(const point &p1, const point &p2, curve &equilibrium, curve &perturbed, const point &mAxis, const point &xPoint) : pBoundary_1(p1), pBoundary_2(p2), equilibrium(equilibrium), perturbed(perturbed), magAxis(mAxis), xPoint(xPoint) {
    // has to be called first
    this->getBoundaries(this->equilibrium, this->perturbed);
    // has to be called after getBoundaries
    this->getMidpoint();
    this->getPerimeter();
    this->getArea();
    this->getHParameter();
    this->getReferenceAngle(); // has to be called after getMidpoint
}

// get midpoint of the lobe over equilibrium curve
void lobe::getMidpoint() {
    // check if curveBoundary has at least 3 points
    if (this->cBoundary_eq.curvePoints.size() < 3) {
        std::cerr << "Error: The lobe's boudary curve has less than 3 points." << std::endl;
        return;
    }
    // calculate the midpoint
    double R = 0.0, Z = 0.0;
    for (size_t i = 0; i < cBoundary_eq.curvePoints.size() - 1; ++i) {
        R += cBoundary_eq.curvePoints[i].R;
        Z += cBoundary_eq.curvePoints[i].Z;
    }
    R /= cBoundary_eq.curvePoints.size() - 1;
    Z /= cBoundary_eq.curvePoints.size() - 1;

    this->midpoint = point(R, Z);
}

// get boundary curve of the lobe between two points for equilibrium and perturbed sets
void lobe::getBoundaries(curve &equilibrium, curve &perturbed) {
    point p1 = this->pBoundary_1;
    point p2 = this->pBoundary_2;
    // Start with the boundary points
    this->cBoundary_eq.curvePoints.push_back(p1);
    this->cBoundary_ptb.curvePoints.push_back(p1);
    // Set flags to check if point i is on the boundary
    bool on_boundary_eq = false;
    bool on_boundary_ptb = false;

    // check all segments of the equilibrium curve to find in which bounding box p1 is
    for (size_t i = 0; i < equilibrium.curvePoints.size() - 1; ++i) {
        segment s(equilibrium.curvePoints[i], equilibrium.curvePoints[i + 1]);
        boundingBox bb(s);
        if (!on_boundary_eq && bb.contains(p1)) {
            // add the segment to the boundary curve
            this->cBoundary_eq.curvePoints.push_back(s.p2);
            on_boundary_eq = true;
        } else if (on_boundary_eq) {
            // add the segment to the boundary curve
            if (!bb.contains(p2)) {
                this->cBoundary_eq.curvePoints.push_back(s.p2);
            } else {
                break;
            }
        }
    }
    cBoundary_eq.curvePoints.push_back(p2);

    // check all segments of the perturbed curve to find in which bounding box p1 is
    for (size_t i = 0; i < perturbed.curvePoints.size() - 1; ++i) {
        segment s(perturbed.curvePoints[i], perturbed.curvePoints[i + 1]);
        boundingBox bb(s);
        if (!on_boundary_ptb && bb.contains(p1)) {
            // add the segment to the boundary curve
            this->cBoundary_ptb.curvePoints.push_back(s.p2);
            on_boundary_ptb = true;
        } else if (on_boundary_ptb) {
            // add the segment to the boundary curve
            if (!bb.contains(p2)) {
                this->cBoundary_ptb.curvePoints.push_back(s.p2);
            } else {
                break;
            }
        }
    }
    cBoundary_ptb.curvePoints.push_back(p2);

    // fill curveBoundary with the boundary points of cBoundary_eq and cBoundary_ptb forming a closed curve
    this->curveBoundary.curvePoints.reserve(cBoundary_eq.curvePoints.size() + cBoundary_ptb.curvePoints.size() - 2);
    // Insert all elements from cBoundary_eq.curvePoints
    curveBoundary.curvePoints.insert(curveBoundary.curvePoints.end(),
                                     cBoundary_eq.curvePoints.begin(),
                                     cBoundary_eq.curvePoints.end());

    // Insert all elements from cBoundary_ptb.curvePoints, excluding the first and last
    if (cBoundary_ptb.curvePoints.size() > 2) {
        curveBoundary.curvePoints.insert(curveBoundary.curvePoints.end(),
                                         std::next(cBoundary_ptb.curvePoints.rbegin()), // skip the last element
                                         cBoundary_ptb.curvePoints.rend());
    }
}

// get perimeter of the lobe
void lobe::getPerimeter() {
    // check if curveBoundary has at least 3 points
    if (curveBoundary.curvePoints.size() < 3) {
        this->perimeter = 0;
        std::cerr << "Error: The lobe's boudary curve has less than 3 points." << std::endl;
        return;
    }
    // check if curveBoundary is closed
    if (curveBoundary.curvePoints.front().distanceTo(curveBoundary.curvePoints.back()) > 1e-14) {
        this->perimeter = 0;
        std::cerr << "Error: The lobe's boudary curve is not closed." << std::endl;
    }
    // calculate the perimeter
    this->perimeter = 0;
    for (size_t i = 0; i < curveBoundary.curvePoints.size() - 1; ++i) {
        this->perimeter += curveBoundary.curvePoints[i].distanceTo(curveBoundary.curvePoints[i + 1]);
    }
}

// get area of the lobe by shoelace formula
void lobe::getArea() {
    // check if curveBoundary has at least 3 points
    if (curveBoundary.curvePoints.size() < 3) {
        this->area = 0;
        std::cerr << "Error: The lobe's boudary curve has less than 3 points." << std::endl;
        return;
    }
    // check if curveBoundary is closed
    if (curveBoundary.curvePoints.front().distanceTo(curveBoundary.curvePoints.back()) > 1e-14) {
        this->area = 0;
        std::cerr << "Error: The lobe's boudary curve is not closed." << std::endl;
    }
    // calculate the area
    this->area = 0.0;
    for (size_t i = 0; i < curveBoundary.curvePoints.size() - 1; ++i) {
        this->area += curveBoundary.curvePoints[i].R * curveBoundary.curvePoints[i + 1].Z -
                      curveBoundary.curvePoints[i + 1].R * curveBoundary.curvePoints[i].Z;
    }
    this->area = std::abs(this->area) / 2.0;
}

// get h parameter of the lobe
void lobe::getHParameter() {
    double baseLength = this->pBoundary_1.distanceTo(this->pBoundary_2);
    this->hParameter = 2 * this->area / baseLength;
}

// get the angle between the lobe's midpoint and the X-point with respect to the magnetic axis
void lobe::getReferenceAngle() {
    // this->referencePoint.R = (this->pBoundary_2.R - this->pBoundary_1.R) / 2.0;
    // this->referencePoint.Z = (this->pBoundary_2.Z - this->pBoundary_1.Z) / 2.0;
    this->referenceAngle = midpoint.angleTo(this->xPoint, this->magAxis);
}