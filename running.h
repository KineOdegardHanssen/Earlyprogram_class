#ifndef RUNNING_H
#define RUNNING_H
#include <systems.h>
#include <diagonalization.h>
#include <quantities.h>
#include <armadillo>
#include <Eigen/Dense>

class Running
{
public:

    int systemsize, sectorsize, maxit;
    double h, J, tolerance;
    bool armadillobool;

    Systems system;


    // Initializers
    Running(bool armadillobool);
    Running(int systemsize, int maxit, double h, double J, double tolerance, bool armadillobool); // Do this some other way...


    // Calls to different kinds of systems. These will be called from some overlying function, or manually after the second initializer
    void homogenous_field(int systemsize, double h, double J, bool armadillobool);
    void alternating_field(int systemsize, double h, double J, bool armadillobool);
    void random_uniform_field(int systemsize, double h, double J, bool armadillobool);

};

#endif // RUNNING_H
