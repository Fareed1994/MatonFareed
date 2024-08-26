#include "RungeKutta.h"

using namespace std;

// Function for a single Runge-Kutta step
void rungeKuttaStep(const VectorXcd &y, const VectorXcd &dydx, double x, double h, VectorXcd &yout, VectorXcd &yerr, VectorXcd (*RGEsystem)(double, const VectorXcd &)) {
    static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
                        b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0,
                        b41 = 0.3, b42 = -0.9, b43 = 1.2,
                        b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
                        b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                        b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
                        c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
                        dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
                        dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.0 / 14336.0, dc6 = c6 - 0.25;

    int n = 107;
    VectorXcd ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n);

    ytemp = y + b21 * h * dydx;
    ak2 = (*RGEsystem)(x + a2 * h, ytemp);

    for (int i = 0; i < n; i++)
        ytemp(i) = y(i) + h * (b31 * dydx(i) + b32 * ak2(i));
    ak3 = (*RGEsystem)(x + a3 * h, ytemp);

    for (int i = 0; i < n; i++)
        ytemp(i) = y(i) + h * (b41 * dydx(i) + b42 * ak2(i) + b43 * ak3(i));
    ak4 = (*RGEsystem)(x + a4 * h, ytemp);

    for (int i = 0; i < n; i++)
        ytemp(i) = y(i) + h * (b51 * dydx(i) + b52 * ak2(i) + b53 * ak3(i) + b54 * ak4(i));
    ak5 = (*RGEsystem)(x + a5 * h, ytemp);

    for (int i = 0; i < n; i++)
        ytemp(i) = y(i) + h * (b61 * dydx(i) + b62 * ak2(i) + b63 * ak3(i) + b64 * ak4(i) + b65 * ak5(i));
    ak6 = (*RGEsystem)(x + a6 * h, ytemp);

    for (int i = 0; i < n; i++)
        yout(i) = y(i) + h * (c1 * dydx(i) + c3 * ak3(i) + c4 * ak4(i) + c6 * ak6(i));
    for (int i = 0; i < n; i++)
        yerr(i) = h * (dc1 * dydx(i) + dc3 * ak3(i) + dc4 * ak4(i) + dc5 * ak5(i) + dc6 * ak6(i));
}

// Function for adaptive step size control using Runge-Kutta
int odeStepper(VectorXcd &y, const VectorXcd &dydx, double *x, double htry, double eps, VectorXcd &yscal, double *hdid, double *hnext, VectorXcd (*RGEsystem)(double, const VectorXcd &)) {
    const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

    int i, n = 107;
    double errmax, h, htemp, xnew;
    VectorXcd yerr(n), ytemp(n);  // Change VectorXd to VectorXcd here
    h = htry;
    for (;;) {
        rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, RGEsystem);
        errmax = 0.0;
        for (i = 0; i < n; i++)
            errmax = std::max(errmax, std::abs(yerr(i) / yscal(i)));
        errmax /= eps;

        if (errmax <= 1.0)
            break;

        htemp = SAFETY * h * pow(errmax, PSHRNK);
        h = (h >= 0.0 ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h));
        xnew = (*x) + h;
        if (xnew == *x) {
            std::cerr << "stepsize underflow in odeStepper\n";
            return 1;
        }
    }

    if (errmax > ERRCON)
        *hnext = SAFETY * h * pow(errmax, PGROW);
    else
        *hnext = 5.0 * h;

    *x += (*hdid = h);
    y = ytemp;
    return 0;
}

// Main integration function
int integrateOdes(VectorXcd &ystart, double from, double to, double eps, double h1, double hmin, VectorXcd (*RGEsystem)(double, const VectorXcd &)) {
    int nvar = 107;
    int nstp, i;
    double x, hnext, hdid, h;
    VectorXcd yscal(nvar), y(ystart), dydx(nvar);

    x = from;
    h = (to > from) ? fabs(h1) : -fabs(h1);

    const int MAXSTP = 5000;  // Increase the maximum number of steps
    const double TINY = 1.0e-20;  // Adjust the tiny value for more precision

    for (nstp = 1; nstp <= MAXSTP; nstp++) {
        dydx = (*RGEsystem)(x, y);
        for (i = 0; i < nvar; i++)
            yscal(i) = std::abs(y(i)) + std::abs(dydx(i) * h) + TINY;

        if ((x + h - to) * (x + h - from) > 0.0)
            h = to - x;

        if (odeStepper(y, dydx, &x, h, eps, yscal, &hdid, &hnext, RGEsystem))
            return 1;

        if ((x - to) * (to - from) >= 0.0) {
            ystart = y;
            return 0;
        }

        if (fabs(hnext) <= hmin) {
            std::cerr << "Step size too small in integrateOdes\n";
            return 1;
        }

        h = hnext;
    }

    std::cerr << "Too many steps in integrateOdes\n";
    return 1;
}




