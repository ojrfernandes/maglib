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
    void getBoundaries(const point &p1, const point &p2, curve &equilibrium, curve &perturbed);
};

// Class constructor
lobe::lobe(const point &p1, const point &p2) : pBoundary_1(p1), pBoundary_2(p2) {
    this->center = point((p1.R + p2.R) / 2, (p1.Z + p2.Z) / 2);
}

// get boundary curve of the lobe between two points for equilibrium and perturbed sets
void lobe::getBoundaries(const point &p1, const point &p2, curve &equilibrium, curve &perturbed) {
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
                                         cBoundary_ptb.curvePoints.rend());             // skip last element
    }
}

int main() {
    std::string equilibriumFile = "/home/jfernandes/Software/maglib/maglit/tcabr_manifolds/manifolds/equilibrium/eq022_S.dat";
    std::string perturbedFile = "/home/jfernandes/Software/maglib/maglit/tcabr_manifolds/manifolds/response/eq022_S.dat";

    curve equilibrium(equilibriumFile);
    curve perturbed(perturbedFile);

    std::cout << "Equilibrium points: " << equilibrium.curvePoints.size() << std::endl;
    std::cout << "Perturbed points: " << perturbed.curvePoints.size() << std::endl;

    // perturbed.xPoint = point(0.46950445, -0.21327174);
    std::vector<point> intersection = perturbed.intersectionsWith(equilibrium);

    // Write the intersection points to a file
    std::ofstream outFile("intersection.dat");
    for (const auto &point : intersection) {
        outFile << point.R << " " << point.Z << std::endl;
    }

    lobe lobetest;
    lobetest.getBoundaries(intersection[407], intersection[408], equilibrium, perturbed);

    // Write the boundary points to a file
    std::ofstream outFile2("boundary.dat");
    for (const auto &point : lobetest.curveBoundary.curvePoints) {
        outFile2 << point.R << " " << point.Z << std::endl;
    }

    return 0;
}