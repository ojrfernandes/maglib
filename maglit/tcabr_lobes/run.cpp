#include "lobes.h"
#include <iomanip>

int main() {
    std::string equilibriumFile = "/home/jfernandes/Software/maglib/maglit/tcabr_manifolds/manifolds/equilibrium/eq022_S.dat";
    std::string perturbedFile = "/home/jfernandes/Software/maglib/maglit/tcabr_manifolds/manifolds/response/eq022_S.dat";

    curve equilibrium(equilibriumFile);
    curve perturbed(perturbedFile);

    std::cout << "Equilibrium points: " << equilibrium.curvePoints.size() << std::endl;
    std::cout << "Perturbed points: " << perturbed.curvePoints.size() << std::endl;

    // Find the intersection points between the equilibrium and perturbed curves
    std::vector<point> intersection = perturbed.intersectionsWith(equilibrium);

    // Write the intersection points to a file with 15 decimal places
    std::ofstream intersectionFile("intersection.dat");
    for (size_t i = 0; i < intersection.size(); ++i) {
        intersectionFile << std::setprecision(15) << intersection[i].R << " " << intersection[i].Z << std::endl;
    }

    std::ofstream midpointFile("midpoint.dat");
    for (size_t i = 0; i < intersection.size() - 1; ++i) {
        lobe lobe_i(intersection[i], intersection[i + 1], equilibrium, perturbed);
        std::cout << "Lobe " << i << " perimeter: " << lobe_i.perimeter
                  << " area: " << lobe_i.area
                  << " midpoint: (" << lobe_i.midpoint.R << " " << lobe_i.midpoint.Z << ")"
                  << std::endl;
        // write midpoint of the lobe to a file
        midpointFile << std::setprecision(15) << lobe_i.midpoint.R << " " << lobe_i.midpoint.Z << std::endl;
    }

    return 0;
}