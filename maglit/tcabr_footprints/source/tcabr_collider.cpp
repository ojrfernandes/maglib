#include "tcabr_collider.h"

bool tcabr_shape::load_shape(const std::string path) {
    // open file from path and handle errors
    std::ifstream f0(path);
    if (!f0) {
        std::cerr << "Error loading vessel shape from file: " << path << std::endl;
        return false;
    }

    // the files first line contains it's integer number of rows
    if (!(f0 >> this->np)) {
        std::cerr << "Error reading 'np' from vessel shape file" << std::endl;
        return false;
    }

    // allocate memory for the read values
    this->R.resize(this->np);
    this->Z.resize(this->np);
    this->th.resize(this->np + 1);
    // initialize auxiliary variables
    double Rmin = 0, Rmax = 0, Zmin = 0, Zmax = 0;

    // search for the minimum and maximun R and Z values
    for (int i = 0; i < (this->np); i++) {
        if (!(f0 >> this->R[i])) {
            std::cerr << "Error reading 'R[" << i << "]' from file" << std::endl;
            return false;
        }

        if (!(f0 >> this->Z[i])) {
            std::cerr << "Error reading 'Z[" << i << "]' from file" << std::endl;
            return false;
        }

        if (this->R[i] < Rmin || i == 0)
            Rmin = this->R[i];
        if (this->R[i] > Rmax || i == 0)
            Rmax = this->R[i];
        if (this->Z[i] < Zmin || i == 0)
            Zmin = this->Z[i];
        if (this->Z[i] > Zmax || i == 0)
            Zmax = this->Z[i];
    }

    // compute the (R,Z) coordinate of the center point in the vessel
    this->Rc = 0.5 * (Rmin + Rmax);
    this->Zc = 0.5 * (Zmin + Zmax);

    // compute the angle between Zi-Zc and Ri-Rc for every i < np
    for (int i = 0; i < (this->np); i++) {
        this->th[i] = atan2(this->Z[i] - this->Zc, this->R[i] - this->Rc);
    }

    // set the last value as angularly periodic
    this->th[this->np] = this->th[0] + 2 * M_PI;

    // initialize index to the first sector
    this->idx = 0;

    f0.close();
    return true;
}

// search the index of a given angle
bool tcabr_shape::search_index(double ang) {
    // check if the angle is in the current sector
    if (ang >= this->th[this->idx] && ang < this->th[this->idx + 1]) {
        return true;
        // check if the angle is in the previous sector
    } else if (ang < this->th[this->idx]) {
        if (this->idx == 0)
            this->idx = this->np - 1;
        else
            this->idx--;
        return search_index(ang);
        // check if the angle is in the next sector
    } else if (ang > this->th[this->idx + 1]) {
        if (this->idx == this->np - 1)
            this->idx = 0;
        else
            this->idx++;
        return search_index(ang);
    } else {
        std::cerr << "tcabr::search_index undefined behavior" << std::endl;
        return false;
    }
}

// calculate cross product AB x AC
double tcabr_shape::cross(double RA, double ZA, double RB, double ZB, double RC, double ZC) {
    double dR1 = RB - RA;
    double dZ1 = ZB - ZA;
    double dR2 = RC - RA;
    double dZ2 = ZC - ZA;
    return dR1 * dZ2 - dR2 * dZ1;
}

// check if coordinate lies inside the vessel
bool tcabr_shape::tcabr_inside(double R, double Z, double phi, void *aux) {
    tcabr_shape *shape = (tcabr_shape *)aux;

    // compute the angle between Z-Zc and R-Rc
    double ang = atan2(Z - shape->Zc, R - shape->Rc);
    if (!search_index(ang))
        return false;

    // determine the next index with periodic wrap-around
    size_t next_idx = (shape->idx == shape->np - 1) ? 0 : shape->idx + 1;

    // cross-product calculation
    double axb = cross(shape->R[shape->idx],
                       shape->Z[shape->idx],
                       shape->R[next_idx],
                       shape->Z[next_idx], R, Z);
    return (axb >= 0);
}