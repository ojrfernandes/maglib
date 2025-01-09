#include "lobes.h"
#include <iomanip>

int main() {
    std::string equilibriumFile = "/home/jfernandes/Software/results/manifolds/equilibrium/eq000_EQ.dat";
    std::string perturbedFile = "/home/jfernandes/Software/results/manifolds/response/eq000_S.dat";

    curve equilibrium(equilibriumFile);
    curve perturbed(perturbedFile);

    std::cout << "Equilibrium points: " << equilibrium.curvePoints.size() << std::endl;
    std::cout << "Perturbed points: " << perturbed.curvePoints.size() << std::endl;

    // Intersection points between the equilibrium and perturbed curves
    std::vector<point> intersection = perturbed.intersectionsWith(equilibrium);

    // Write the intersection points to a file with 15 decimal places
    std::ofstream intersectionFile("intersection.dat");
    for (size_t i = 0; i < intersection.size(); ++i) {
        intersectionFile << std::setprecision(15) << intersection[i].R << " " << intersection[i].Z << std::endl;
    }

    std::ofstream lobeFile("lobe.dat");

    // Write the header to the lobe file
    lobeFile << "# LobeIndex    Perimeter    Area    H-Parameter" << std::endl;

    for (size_t i = 0; i < intersection.size() - 1; ++i) {
        lobe lobe_i(intersection[i], intersection[i + 1], equilibrium, perturbed);
        std::cout << "Lobe " << i << " perimeter: " << lobe_i.perimeter
                  << " area: " << lobe_i.area << " h-parameter: " << lobe_i.hParameter
                  << std::endl;
        // write lobe data to a file
        lobeFile << std::setprecision(15) << i << " " << lobe_i.perimeter
                 << " " << lobe_i.area << " " << lobe_i.hParameter << std::endl;
    }

    return 0;
}