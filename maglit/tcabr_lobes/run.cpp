#include "hdf5.h"
#include "lobes.h"
#include <iomanip>

// read coordinates of xmag and zmag from hdf5 file
point magAxis(const std::string &hdf5File) {
    // paths to the xmag and zmag datasets in the HDF5 file
    const std::string xmagPath = "scalars/xmag";
    const std::string zmagPath = "scalars/zmag";

    // Open the HDF5 file
    hid_t file = H5Fopen(hdf5File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Open the datasets
    hid_t xmagDataset = H5Dopen(file, xmagPath.c_str(), H5P_DEFAULT);
    hid_t zmagDataset = H5Dopen(file, zmagPath.c_str(), H5P_DEFAULT);

    // Get the dataspace of the datasets
    hid_t xmagSpace = H5Dget_space(xmagDataset);
    hid_t zmagSpace = H5Dget_space(zmagDataset);

    // Get the number of elements in the datasets
    hsize_t xmagSize;
    hsize_t zmagSize;
    H5Sget_simple_extent_dims(xmagSpace, &xmagSize, nullptr);
    H5Sget_simple_extent_dims(zmagSpace, &zmagSize, nullptr);

    // Read the data into vectors
    std::vector<double> xmagData(xmagSize);
    std::vector<double> zmagData(zmagSize);

    H5Dread(xmagDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xmagData.data());
    H5Dread(zmagDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zmagData.data());

    // save the xmag and zmag values to R_xmag and Z_xmag
    point magAxis(xmagData[0], zmagData[0]);

    // Close datasets and file
    H5Dclose(xmagDataset);
    H5Dclose(zmagDataset);
    H5Fclose(file);

    return magAxis;
}

int main() {
    std::string hdf5File = "/home/jfernandes/Software/m3dc1_data/eq022/i_coils/n03/C1.h5";
    std::string equilibriumFile = "/home/jfernandes/Software/results/manifolds/equilibrium/eq022_EQ.dat";
    std::string perturbedFile = "/home/jfernandes/Software/results/manifolds/response/eq022_i03_S.dat";

    point magAxisCoord = magAxis(hdf5File);

    std::cout << "magAxisCoord: " << magAxisCoord.R << " " << magAxisCoord.Z << std::endl;

    curve equilibrium(equilibriumFile, 0);
    curve perturbed(perturbedFile, 800);

    std::cout << "Equilibrium points: " << equilibrium.curvePoints.size() << std::endl;
    std::cout << "Perturbed points: " << perturbed.curvePoints.size() << std::endl;

    // Intersection points between the equilibrium and perturbed curves
    std::vector<point> intersection = perturbed.intersectionsWith(equilibrium);

    point xPoint = perturbed.curvePoints[0];

    std::cout << "xPoint: " << xPoint.R << " " << xPoint.Z << std::endl;

    // Write the intersection points to a file with 15 decimal places
    std::ofstream intersectionFile("/home/jfernandes/Software/results/lobes/its_eq022_S.dat");
    for (size_t i = 0; i < intersection.size(); ++i) {
        intersectionFile << std::setprecision(15) << intersection[i].R << " " << intersection[i].Z << std::endl;
    }

    std::ofstream lobeFile("/home/jfernandes/Software/results/lobes/lobe_eq022_S.dat");

    // Write the header to the lobe file
    lobeFile << "#Rmid Zmid polAngle Perimeter Area H-Parameter" << std::endl;

    for (size_t i = 0; i < intersection.size() - 1; ++i) {
        lobe lobe_i(intersection[i], intersection[i + 1], equilibrium, perturbed, magAxisCoord, xPoint);
        double newAngle = lobe_i.referenceAngle;
        lobeFile << std::setprecision(5) << lobe_i.midpoint.R << " " << lobe_i.midpoint.Z << " " << std::setprecision(15) << newAngle << " " << lobe_i.perimeter << " " << lobe_i.area << " " << lobe_i.hParameter << std::endl;
    }

    return 0;
}