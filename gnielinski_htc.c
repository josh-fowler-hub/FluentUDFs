#include "udf.h"
#include <math.h>

// === Property Functions for Air ===
real get_thermal_conductivity(real Tf)
{
    return 0.024 + 7e-5 * (Tf - 273.15); // W/m·K
}

real get_prandtl_number(real Tf)
{
    return 0.71; // Approximate constant for air
}

real get_viscosity(real Tf)
{
    real T_C = Tf - 273.15;
    return 1.458e-6 * pow(Tf, 1.5) / (Tf + 110.4); // Pa·s
}

real get_density(real Tf)
{
    return 101325 / (287.058 * Tf); // kg/m³ using ideal gas law
}

// real lambert_W(real x)
// {
//     // Lambert W function using Taylor/Laurent series expansions
//     // Avoids iterative methods for better performance
    
//     if (x == 0.0) return 0.0; // Exact solution for W(0) = 0
    
//     if (x < -0.35) {
//         // Near the branch point x = -1/e ≈ -0.36788
//         // Use Laurent series expansion around the branch point
//         real p = sqrt(2.0 * (2.718281828459045 * x + 1.0));
//         if (p == 0.0) return -1.0; // Branch point exactly
        
//         // Laurent series: W(x) = -1 + p - p²/3 + 11p³/72 - 43p⁴/540 + 769p⁵/17280 + ...
//         return -1.0 + p * (1.0 - p * (1.0/3.0 - p * (11.0/72.0 - p * (43.0/540.0 - p * 769.0/17280.0))));
//     }
//     else if (x < 0.1) {
//         // For small x, use Taylor series around x = 0
//         // W(x) = x - x² + (3/2)x³ - (8/3)x⁴ + (125/24)x⁵ + ...
//         return x * (1.0 - x * (1.0 - x * (1.5 - x * (8.0/3.0 - x * 125.0/24.0))));
//     }
//     else if (x < 3.0) {
//         // For moderate x, use Padé approximation
//         // More accurate than simple log approximation
//         real lnx = log(x);
//         real lnx_minus_lnlnx = lnx - log(lnx);
//         real correction = log(lnx) / (1.0 + lnx);
//         return lnx_minus_lnlnx + correction;
//     }
//     else if (x < 100.0) {
//         // For large x, use asymptotic series
//         // W(x) = ln(x) - ln(ln(x)) + ln(ln(x))/ln(x) + ...
//         real lnx = log(x);
//         real lnlnx = log(lnx);
//         real ratio = lnlnx / lnx;
//         return lnx - lnlnx + ratio * (1.0 - 0.5 * ratio + ratio * ratio / 3.0);
//     }
//     else {
//         // For very large x, simpler asymptotic form
//         real lnx = log(x);
//         real lnlnx = log(lnx);
//         return lnx - lnlnx + lnlnx / lnx;
//     }
// }

// // === Gnielinski Correlation ===
// real get_friction_factor_blaisus(real Re)
// {
//     return 0.3164 * pow(Re, -0.25); // Blasius equation
// }

// real get_friction_factor_colebrook_white(real Re, real epsilon, real D_h)
// {
//     // Colebrook-White equation using Lambert W function
//     // Following the exact formulation from literature
//     // x = 1/sqrt(f), b = epsilon/(14.8*R_h), a = 2.51/Re
//     // f = 1/[2*W(ln(10)/(2a) * 10^(b/(2a))) / ln(10) - b/a]^2
    
//     real R_h = D_h / 4.0; // Hydraulic radius for circular pipes
//     real a = 2.51 / Re;
//     real b = epsilon / (14.8 * R_h);
    
//     // Calculate the argument for Lambert W function
//     real ln10 = log(10.0);
//     real arg1 = ln10 / (2.0 * a);
//     real arg2 = b / (2.0 * a);
//     real arg_exp = pow(10.0, arg2);
//     real W_argument = arg1 * arg_exp;
    
//     // Solve using Lambert W function
//     real W = lambert_W(W_argument);
    
//     // Calculate x = 1/sqrt(f)
//     real x = (2.0 * W) / ln10 - b / a;
    
//     // Calculate friction factor f = 1/x^2
//     real f = 1.0 / (x * x);
    
//     // Ensure reasonable bounds for friction factor
//     if (f < 0.008) f = 0.008;  // Lower bound
//     if (f > 0.1) f = 0.1;      // Upper bound
    
//     return f;
// }

// real get_friction_factor(real Re)
// {
//     real b = epsilon / (14.8 * (D_h / 4.0));
//     real a = 2.51 / Re;
//     real f = a + b;
//     return f;
// }

// // === Additional Friction Factor Correlations ===

// real get_friction_factor_moody(real Re, real epsilon, real D_h)
// {
//     // Moody (1947) correlation
//     real rel_rough = epsilon / D_h;
//     return 0.0055 * (1.0 + pow(2e4 * rel_rough + 1e6 / Re, 1.0/3.0));
// }

// real get_friction_factor_wood(real Re, real epsilon, real D_h)
// {
//     // Wood (1966) correlation
//     real rel_rough = epsilon / D_h;
//     real psi = 1.62 * pow(rel_rough, 0.134);
//     return 0.094 * pow(rel_rough, 0.225) + 0.53 * rel_rough + 88.0 * pow(rel_rough, 0.44) * pow(Re, -psi);
// }

// real get_friction_factor_eck(real Re, real epsilon, real D_h)
// {
//     // Eck (1973) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = -2.0 * log10(rel_rough / 3.715 + 15.0 / Re);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_swamee_jain(real Re, real epsilon, real D_h)
// {
//     // Swamee and Jain (1976) correlation
//     real rel_rough = epsilon / D_h;
//     real log_term = log10(rel_rough / 3.7 + 5.74 / pow(Re, 0.9));
//     return 0.25 / (log_term * log_term);
// }

// real get_friction_factor_churchill_1973(real Re, real epsilon, real D_h)
// {
//     // Churchill (1973) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = -2.0 * log10(rel_rough / 3.715 + pow(7.0 / Re, 0.9));
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_jain(real Re, real epsilon, real D_h)
// {
//     // Jain (1976) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = -2.0 * log10(rel_rough / 3.715 + pow(6.943 / Re, 0.9));
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_churchill_1977(real Re, real epsilon, real D_h)
// {
//     // Churchill (1977) correlation - all flow regimes
//     real rel_rough = epsilon / D_h;
//     real theta1 = pow(-2.457 * log(pow(7.0 / Re, 0.9) + 0.27 * rel_rough), 16.0);
//     real theta2 = pow(37530.0 / Re, 16.0);
//     return 8.0 * pow(pow(8.0 / Re, 12.0) + 1.0 / pow(theta1 + theta2, 1.5), 1.0/12.0);
// }

// real get_friction_factor_chen(real Re, real epsilon, real D_h)
// {
//     // Chen (1979) correlation
//     real rel_rough = epsilon / D_h;
//     real log_arg = rel_rough / 3.7065 - 5.0452 / Re * log10(1.0 / 2.8257 * pow(rel_rough, 1.1098) + 5.8506 / pow(Re, 0.8981));
//     real inv_sqrt_f = -2.0 * log10(log_arg);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_round(real Re, real epsilon, real D_h)
// {
//     // Round (1980) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = 1.8 * log10(Re / (0.135 * Re * rel_rough + 6.5));
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_haaland(real Re, real epsilon, real D_h)
// {
//     // Haaland (1983) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = -1.8 * log10(pow(rel_rough / 3.7, 1.11) + 6.9 / Re);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_serghides(real Re, real epsilon, real D_h)
// {
//     // Serghides (1984) correlation
//     real rel_rough = epsilon / D_h;
//     real psi1 = -2.0 * log10(rel_rough / 3.7 + 12.0 / Re);
//     real psi2 = -2.0 * log10(rel_rough / 3.7 + 2.51 * psi1 / Re);
//     real psi3 = -2.0 * log10(rel_rough / 3.7 + 2.51 * psi2 / Re);
//     real inv_sqrt_f = psi1 - (psi2 - psi1) * (psi2 - psi1) / (psi3 - 2.0 * psi2 + psi1);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_tsal(real Re, real epsilon, real D_h)
// {
//     // Tsal (1989) correlation
//     real rel_rough = epsilon / D_h;
//     real A = 0.11 * pow(68.0 / Re + rel_rough, 0.25);
//     if (A >= 0.018) {
//         return A;
//     } else {
//         return 0.0028 + 0.85 * A;
//     }
// }

// real get_friction_factor_manadilli(real Re, real epsilon, real D_h)
// {
//     // Manadilli (1997) correlation
//     real rel_rough = epsilon / D_h;
//     real inv_sqrt_f = -2.0 * log10(rel_rough / 3.7 + 95.0 / pow(Re, 0.983) - 96.82 / Re);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

// real get_friction_factor_buzzelli(real Re, real epsilon, real D_h)
// {
//     // Buzzelli (2008) correlation
//     real a = 1.0 / (1.0 + pow(Re / 2720.0, 9.0));
//     real b = 1.0 / (1.0 + pow(Re / (160.0 * D_h / epsilon), 2.0));
    
//     real term1 = pow(Re / 64.0, a);
//     real term2 = pow(1.8 * log10(Re / 6.8), 2.0 * (1.0 - a) * b);
//     real term3 = pow(2.0 * log10(3.7 * D_h / epsilon), 2.0 * (1.0 - a) * (1.0 - b));
    
//     return 1.0 / (term1 * term2 * term3);
// }

// real get_friction_factor_fang(real Re, real epsilon, real D_h)
// {
//     // Fang (2011) correlation
//     real rel_rough = epsilon / D_h;
//     real beta = log(Re / (1.816 * log(1.1 * Re / log(1.0 + 1.1 * Re))));
//     real inv_sqrt_f = -2.0 * log10(2.18 * beta / Re + rel_rough / 3.71);
//     return 1.0 / (inv_sqrt_f * inv_sqrt_f);
// }

real get_friction_factor_goudar_sonnad(real Re, real epsilon, real D_h)
{
    // Goudar-Sonnad (2006) correlation - exact formulation
    real a = 2.0 / log(10.0);
    real b = (epsilon / D_h) / 3.7;
    real d = log(10.0) * Re / 5.02;
    real s = b * d + log(d);
    real q = pow(s, s / (s + 1.0));
    real g = b * d + log(d / q);
    real z = log(q / g);
    real DLA = z * g / (g + 1.0);
    real DCFA = DLA * (1.0 + (z / 2.0) / ((g + 1.0) * (g + 1.0) + (z / 3.0) * (2.0 * g - 1.0)));
    
    real inv_sqrt_f = a * (log(d / q) + DCFA);
    return 1.0 / (inv_sqrt_f * inv_sqrt_f);
}

real get_nusselt_gnielinski(real Re, real Pr, real epsilon, real D_h)
{
    real f = get_friction_factor_goudar_sonnad(Re, epsilon, D_h);
    real numerator = (f / 8.0) * (Re - 1000.0) * Pr;
    real denominator = 1.0 + 12.7 * sqrt(f / 8.0) * (pow(Pr, 2.0/3.0) - 1.0);
    return numerator / denominator;
}

// === Main HTC UDF ===
DEFINE_PROFILE(heat_transfer_coefficient, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, rho, mu, k, Pr, Re, vel, Nu, h;
    const real D_h = 0.0127; // Hydraulic diameter in meters
    const real epsilon = 1.5e-6; // Surface roughness in meters (typical for smooth pipes)

    begin_f_loop(f, t)
    {
        Tw = F_T(f, t);            // Wall temperature
        Tinf = F_PROFILE(f, t, i); // Bulk temperature from BC
        Tf = 0.5 * (Tw + Tinf);    // Film temperature

        // Get fluid properties at film temperature
        k = get_thermal_conductivity(Tf);
        Pr = get_prandtl_number(Tf);
        mu = get_viscosity(Tf);
        rho = get_density(Tf);

        // Get velocity from adjacent fluid cell
        c0 = F_C0(f, t);
        t0 = THREAD_T0(t);
        vel = sqrt(pow(C_U(c0, t0), 2) + pow(C_V(c0, t0), 2) + pow(C_W(c0, t0), 2));

        // Reynolds number
        Re = rho * vel * D_h / mu;

        // Nusselt number via Gnielinski with Goudar-Sonnad friction factor
        Nu = get_nusselt_gnielinski(Re, Pr, epsilon, D_h);

        // Heat transfer coefficient
        h = Nu * k / D_h;

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}
