#include "udf.h"
#include <math.h>

real get_thermal_conductivity(real Tf)
{
    return 0.6065 - 0.001484 * (Tf - 273.15);
}

real get_prandtl_number(real Tf)
{
    return 13.0 - 0.03 * (Tf - 273.15);
}

real get_viscosity(real Tf)
{
    real T_C = Tf - 273.15;
    return 2.414e-5 * pow(10, 247.8 / (T_C + 133.15));
}

real get_density(real Tf)
{
    real T_C = Tf - 273.15;
    return 1000 * (1 - ((T_C + 288.9414) / (508929.2 * (T_C + 68.12963))) * pow(T_C - 3.9863, 2));
}

DEFINE_PROFILE(heat_transfer_coefficient, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, rho, mu, k, Pr, Re, vel, D_h = 0.01; // Example hydraulic diameter
    real h;

    begin_f_loop(f, t)
    {
        Tw = F_T(f, t); // Wall temperature
        Tinf = F_PROFILE(f, t, i); // T∞ from BC profile
        Tf = 0.5 * (Tw + Tinf);

        k = get_thermal_conductivity(Tf);
        Pr = get_prandtl_number(Tf);
        mu = get_viscosity(Tf);
        rho = get_density(Tf);

        c0 = F_C0(f, t); // Adjacent fluid cell
        t0 = THREAD_T0(t);
        vel = sqrt(pow(C_U(c0, t0), 2) + pow(C_V(c0, t0), 2) + pow(C_W(c0, t0), 2));

        Re = rho * vel * D_h / mu;

        h = 0.023 * pow(Re, 0.8) * pow(Pr, 0.3) * k / D_h;

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}
