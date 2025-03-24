#include "hdf5.h"
#include "input_read.h"
#include "lobe.h"
#include <iomanip>

int main() {
    // read params from input file
    std::string pathsFile = "params.txt";

    // read input file
    input_read input(pathsFile);
    bool readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    // read xmag and zmag from the HDF5 file
    bool hdf5Status = input.readHDF5File();
    if (!hdf5Status) {
        std::cerr << "Error reading HDF5 file." << std::endl;
        return 1;
    }

    point magAxisCoord(input.xmag, input.zmag);
    std::cout << "magAxisCoord: " << magAxisCoord.R << " " << magAxisCoord.Z << std::endl;

    curve equilibrium(input.equilibriumFile, 0);
    curve perturbed(input.perturbedFile, 800);

    std::cout << "Equilibrium points: " << equilibrium.curvePoints.size() << std::endl;
    std::cout << "Perturbed points: " << perturbed.curvePoints.size() << std::endl;

    // Intersection points between the equilibrium and perturbed curves
    std::vector<point> intersection = perturbed.intersectionsWith(equilibrium);

    point xPoint = perturbed.curvePoints[0];

    std::cout << "xPoint: " << xPoint.R << " " << xPoint.Z << std::endl;

    // Write the intersection points to a file with 15 decimal places
    std::ofstream intersectionFile(input.intersectionFile);
    for (size_t i = 0; i < intersection.size(); ++i) {
        intersectionFile << std::setprecision(15) << intersection[i].R << " " << intersection[i].Z << std::endl;
    }

    std::ofstream lobeFile(input.lobeFile);

    // Write the header to the lobe file
    lobeFile << "#Rmid Zmid polAngle Perimeter Area H-Parameter" << std::endl;

    for (size_t i = 0; i < intersection.size() - 1; ++i) {
        lobe lobe_i(intersection[i], intersection[i + 1], equilibrium, perturbed, magAxisCoord, xPoint);
        double newAngle = lobe_i.referenceAngle;
        lobeFile << std::setprecision(5) << lobe_i.midpoint.R << " " << lobe_i.midpoint.Z << " " << std::setprecision(15) << newAngle << " " << lobe_i.perimeter << " " << lobe_i.area << " " << lobe_i.hParameter << std::endl;
    }

    return 0;
}