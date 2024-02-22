#include "tcabr_collider.h"

bool load_shape(const char path[], tcabr_shape *shape) {
    FILE *f0 = fopen(path, "r");
    if (!f0) {
        printf("load_shape: %s not found\n", path);
        return false;
    } else {
        printf("load_shape: %s loaded\n", path);
    }

    // Check the return value of fscanf and handle errors
    if (fscanf(f0, "%d", &(shape->np)) != 1) {
        fclose(f0);
        fprintf(stderr, "Error reading 'np' from file\n");
        return false;
    }

    shape->R = new double[shape->np];
    shape->Z = new double[shape->np];
    shape->th = new double[shape->np + 1];
    double Rmin = 0, Rmax = 0, Zmin = 0, Zmax = 0;

    for (int i = 0; i < (shape->np); i++) {
        if (fscanf(f0, "%lf", shape->R + i) != 1) {
            fclose(f0);
            fprintf(stderr, "Error reading 'R[%d]' from file\n", i);
            return false;
        }

        if (fscanf(f0, "%lf", shape->Z + i) != 1) {
            fclose(f0);
            fprintf(stderr, "Error reading 'Z[%d]' from file\n", i);
            return false;
        }

        if (shape->R[i] < Rmin || i == 0)
            Rmin = shape->R[i];
        if (shape->R[i] > Rmax || i == 0)
            Rmax = shape->R[i];
        if (shape->Z[i] < Zmin || i == 0)
            Zmin = shape->Z[i];
        if (shape->Z[i] > Zmax || i == 0)
            Zmax = shape->Z[i];
    }

    shape->Rc = 0.5 * (Rmin + Rmax);
    shape->Zc = 0.5 * (Zmin + Zmax);
    printf("shape Rc, Zc = %f, %f\n", shape->Rc, shape->Zc);

    for (int i = 0; i < (shape->np); i++) {
        shape->th[i] = atan2(shape->Z[i] - shape->Zc, shape->R[i] - shape->Rc);
    }

    shape->th[shape->np] = shape->th[0] + 2 * M_PI;
    fclose(f0);
    return true;
}

void free_shape(tcabr_shape *shape) {
    delete[] shape->R;
    delete[] shape->Z;
    delete[] shape->th;
    delete shape;
}

bool search_index(double ang, tcabr_shape *shape) {
    // printf("search_index: theta[%d] = %f, ang = %f\n", shape->idx, shape->th[shape->idx], ang);
    if (ang >= shape->th[shape->idx] && ang < shape->th[shape->idx + 1]) {
        return true;
    } else if (ang < shape->th[shape->idx]) {
        if (shape->idx == 0)
            shape->idx = shape->np - 1;
        else
            shape->idx--;
        return search_index(ang, shape);
    } else if (ang > shape->th[shape->idx + 1]) {
        if (shape->idx == shape->np - 1)
            shape->idx = 0;
        else
            shape->idx++;
        return search_index(ang, shape);
    } else {
        printf("search_index: undefined behavior\n");
        return false;
    }
}

// calculate AB x AC
double cross(double RA, double ZA, double RB, double ZB, double RC, double ZC) {
    double dR1 = RB - RA;
    double dZ1 = ZB - ZA;
    double dR2 = RC - RA;
    double dZ2 = ZC - ZA;
    return dR1 * dZ2 - dR2 * dZ1;
}

bool tcabr_inside(double R, double Z, double phi, void *aux) {
    tcabr_shape *shape = (tcabr_shape *)aux;
    double       ang = atan2(Z - shape->Zc, R - shape->Rc);
    if (!search_index(ang, shape))
        return false;

    double axb = cross(shape->R[shape->idx],
                       shape->Z[shape->idx],
                       shape->R[shape->idx + 1],
                       shape->Z[shape->idx + 1], R, Z);
    return (axb >= 0);
}