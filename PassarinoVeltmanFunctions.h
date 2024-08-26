#ifndef PassarinoVeltmanFunctions_h
#define PassarinoVeltmanFunctions_h
#include <iostream>
#include <math.h>
#include "ChiSquared.h"

double A0(double M);
double B0(double p, double m1, double m2);
double B1(double p, double M, double m);
double B22(double p, double m1, double m2);
double tB22(double p, double m1, double m2);
double F(double p, double m1, double m2);
double G(double p, double m1, double m2);
double H(double p, double m1, double m2);
double C0(double m1, double m2, double m3);
double D0(double m1, double m2, double m3, double m4);

#endif /* PassarinoVeltmanFunctions_h */
