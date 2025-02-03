#ifndef TCABR_COLLIDER
#define TCABR_COLLIDER
#define TCABR_COLLIDER_V 241119 // version (yy.mm.dd)

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class tcabr_shape {
  public:
    // class default constructor
    tcabr_shape() = default;
    // class constructor loading shape from file
    tcabr_shape(const std::string path) {
        load_shape(path);
    }

    // check if coordinate lies inside the vessel
    bool tcabr_inside(double R, double Z, double phi, void *aux);

  private:
    // load vessel shape from file
    bool load_shape(const std::string path);
    // search the index of a given angle
    bool search_index(double ang);
    // calculate AB x AC
    double cross(double RA, double ZA, double RB, double ZB, double RC, double ZC);

    int np;                 // number of nodes
    int idx;                // current index
    std::vector<double> R;  // node coordinates
    std::vector<double> Z;  //
    double Rc;              // center coordinates
    double Zc;              //
    std::vector<double> th; // sector angles
};

#endif