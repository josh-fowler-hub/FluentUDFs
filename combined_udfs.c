#include "udf.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// === Boundary Conditions ===
#define UDM_TINF 0                                                  // User Defined Memory location for Tinf
#define UDM_HTC 1                                                   // User Defined Memory location for HTC
#define UDM_MASS_FLOW 2                                       // UDM location for mass flow rate

// --- Dynamic Radiative and Stack Loss Tracking ---
#define UDM_RAD_WALL_LOSS 3      // UDM location for radiative loss
#define UDM_RAD_GAS_LOSS 4      // UDM location for gas radiative loss
#define UDM_STACK_LOSS_PERCENT 5    // UDM location for stack loss
#define UDM_STACK_LOSS 6    // UDM location for total stack loss
#define UDM_TOTAL_LOSS 7    // UDM location for total loss
#define UDM_CONDENSED 8 // UDM location for condensed water vapor
#define UDM_UNCONDENSED 9 // UDM location for uncondensed water vapor
#define UDM_CONDENSED_MASS 10 // UDM location for condensed mass

// === Stratified Temperature Profile ===
#define TBOT 294.261111                                             // Bottom temperature in Kelvin
#define TTOP 313.705556                                             // Top temperature in Kelvin
#define Y_MIN 0.00                                                  // Minimum Y-coordinate for this zone, meters
#define Y_MAX 1.4816836                                             // Maximum Y-coordinate for this zone, meters

// === Hydraulic Diameter ===
#define A1 0.194056*M_PI                                             // Cross-sectional area in m^2
#define P1 1.83495639                                               // Perimeter in meters
#define L1 1.2456127                                               // Length in meters
#define A11 0.143256*M_PI                                             // Cross-sectional area in m^2
#define P11 0.47969984                                             // Perimeter in meters
#define L11 0.01482392                                               // Length in meters
#define A12 ((0.182372+0.039116)/2.0)*M_PI                                             // Cross-sectional area in m^2
#define P12 0.48376839                                             // Perimeter in meters
#define L12 0.01143641                                               // Length in meters
#define A13 (((2*0.07497168)+(2*0.0947265))/2.0)*M_PI                                             // Cross-sectional area in m^2
#define P13 ((2*0.23553049)+(2*0.29759209))/2.0                                             // Perimeter in meters
#define L13 0.03532738                                               // Length in meters
#define A14 (((0.167132)+(0.026924))/2.0)*M_PI                                             // Cross-sectional area in m^2
#define P14 0.61839156                                             // Perimeter in meters
#define L14 0.00787182                                               // Length in meters
#define A2 0.00258806                                             // Cross-sectional area in m^2
#define P2 0.18033998                                            // Perimeter in meters
#define L2 1.1464569                                            // Length in meters
#define A3 0.00156958                                          // Cross-sectional area in m^2
#define P3 0.14044176                                         // Perimeter in meters
#define L3 1.4816836                                         // Length in meters
#define END_CAP_ROC 0.12701944                                // Measured radius of curvature in meters

// === Default HTC for Water (Used for initialization) ===
#define DEFAULT_HTC_WATER 1136.0                                    // Default heat transfer coefficient for water in W/(m^2*K)

// --- Empirical Radiation and Stack Loss Correlations ---
#define EMISSIVITY_WALL 0.85         // Wall surface emissivity (oxidized steel)
#define EMISSIVITY_GAS 0.3           // Effective flue gas emissivity (CO2/H2O)
#define SIGMA_SB 5.670374419e-8      // Stefan-Boltzmann constant (W/m^2/K^4)
#define T_AMB 298.15                 // Ambient temperature (K)
#define CP_GAS 1563.46 // J/kg/K
#define FLAME_TEMP 1366.48 // K
#define WALL_THICKNESS 0.007 // Example wall thickness in meters, set as needed

#define NUM_RAD_LOSS_BOUNDARIES 2 // Set to the number of boundaries
// WARNING: Ensure these boundary IDs correspond to existing wall boundary zones in your mesh
// TEMPORARILY DISABLED - These boundary IDs may not exist and cause segfaults during UDF loading
// int RAD_LOSS_BOUNDARY_IDS[NUM_RAD_LOSS_BOUNDARIES] = {137, 143}; // Update with correct boundary zone IDs
// static int RAD_LOSS_BOUNDARY_IDS[NUM_RAD_LOSS_BOUNDARIES] = {137, 143}; // Disabled - set to safe values

// --- Condensation detection and mass calculation ---
#define DEW_POINT 332.393408 // K
#define LATENT_HEAT_VAP 2257000 // J/kg
#define WATER_VAPOR_MASS_FRAC 0.123902 // mass fraction

// --- Mass Flow Rate Calculation ---
#define MASS_FLOW_RATE 0.05088578 // kg/s

// Global flag to control UDF behavior during initialization
static int g_udf_enabled = 0;  // Start disabled

// === Utility Functions ===

// Global UDM validation function
int validate_udm_setup(Domain *d)
{
    Thread *t;
    cell_t c;
    int udm_count = 0;
    
    if (d == NULL) {
        Message("Error: Domain is NULL in validate_udm_setup\n");
        return 0;
    }
    
    // Find a fluid thread to check UDM allocation
    thread_loop_c(t, d)
    {
        if (FLUID_THREAD_P(t) && THREAD_STORAGE(t, SV_UDM_I) != NULL)
        {
            begin_c_loop(c, t)
            {
                // Check if we can access our required UDM indices
                if (N_UDM > UDM_CONDENSED_MASS) {
                    udm_count = N_UDM;
                    break;
                }
            }
            end_c_loop(c, t)
            break;
        }
    }
    
    if (udm_count <= UDM_CONDENSED_MASS) {
        Message("Error: Insufficient UDM variables allocated. Need at least %d, found %d\n", 
                UDM_CONDENSED_MASS + 1, udm_count);
        Message("Please allocate at least %d UDM variables in Fluent before using this UDF\n", 
                UDM_CONDENSED_MASS + 1);
        return 0;
    }
    
    Message("UDM validation successful: %d UDM variables allocated\n", udm_count);
    return 1;
}

// Enhanced safe initialization state checker with maximum protection
int is_solution_initialized(face_t f, Thread *t)
{
    if (t == NULL) return 0;
    
    // Check if basic thread storage is available first
    if (!THREAD_STORAGE(t, SV_T)) return 0;
    
    // Check if face thread storage is available
    if (!FACE_THREAD_P(t)) return 0;
    
    // Try to safely get cell information with additional checks
    Thread *t0 = THREAD_T0(t);
    if (t0 == NULL) return 0;
    
    // Check if this is a valid fluid thread
    if (!FLUID_THREAD_P(t0)) return 0;
    
    cell_t c0 = F_C0(f, t);
    if (c0 < 0) return 0;
    
    // Check if cell thread storage is available
    if (!THREAD_STORAGE(t0, SV_T)) return 0;
    
    // Additional safety check - ensure we have valid data
    if (!NULLP(THREAD_STORAGE(t0, SV_T))) {
        real cell_temp = C_T(c0, t0);
        real face_temp = F_T(f, t);
        
        // Check if temperatures are reasonable (indicates solution is active)
        if (cell_temp <= 50.0 || cell_temp > 5000.0 || face_temp <= 50.0 || face_temp > 5000.0) {
            return 0;
        }
    } else {
        return 0;
    }
    
    return 1;
}

// Macro for safe profile initialization check
#define SAFE_PROFILE_CHECK(f, t, default_value) \\\n    do { \\\n        if (!is_solution_initialized(f, t)) { \\\n            F_PROFILE(f, t, i) = default_value; \\\n            continue; \\\n        } \\\n    } while(0)

real calculate_spherical_cap_height(real base_diameter, real sphere_radius)
{
    real base_radius = base_diameter / 2.0;
    real height_cap;
    
    if (base_radius < sphere_radius) {
        height_cap = sphere_radius - sqrt(sphere_radius * sphere_radius - base_radius * base_radius);
    } else {
        height_cap = sphere_radius;
    }
    
    return height_cap;
}

int check_flow_initialization(face_t f, Thread *t, cell_t *c0, Thread **t0)
{
    // Add input validation
    if (t == NULL || c0 == NULL || t0 == NULL) {
        Message("Error: NULL pointer in check_flow_initialization\n");
        return 0;
    }
    
    *c0 = F_C0(f, t);
    *t0 = THREAD_T0(t);
    
    // Basic checks for thread and cell validity
    if (*t0 == NULL || *c0 < 0) {
        return 0; // Invalid thread or cell
    }
    
    // Check if basic storage is available
    if (!THREAD_STORAGE(*t0, SV_T) || !THREAD_STORAGE(*t0, SV_UDM_I)) {
        return 0; // Required storage not available
    }
    
    // Check if temperatures are in reasonable range (indicates solution is initialized)
    real cell_temp = C_T(*c0, *t0);
    real face_temp = F_T(f, t);
    
    if (cell_temp <= 0.0 || cell_temp > 10000.0 || face_temp <= 0.0 || face_temp > 10000.0) {
        return 0; // Temperatures not initialized or unreasonable
    }
    
    return 1; // Flow is properly initialized
}

void set_default_profile_value(face_t f, Thread *t, int i, cell_t c0, Thread *t0, int udm_index)
{
    // Add input validation
    if (t == NULL) {
        Message("Error: Thread is NULL in set_default_profile_value\n");
        return;
    }
    
    // Set the profile value
    F_PROFILE(f, t, i) = DEFAULT_HTC_WATER;

    // Store default value in UDM if possible with bounds checking
    if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL && 
        udm_index >= 0 && udm_index < N_UDM) {
        C_UDMI(c0, t0, udm_index) = DEFAULT_HTC_WATER;
    }
}

real compute_stratified_tinf(real y)
{
    if (y < Y_MIN) y = Y_MIN;
    else if (y > Y_MAX) y = Y_MAX;
    
    real Tinf;
    if ((Y_MAX - Y_MIN) > 1e-12)
        Tinf = TBOT + (TTOP - TBOT) * ((y - Y_MIN) / (Y_MAX - Y_MIN));
    else
        Tinf = TBOT;

    if (Tinf < TBOT) Tinf = TBOT;
    if (Tinf > TTOP) Tinf = TTOP;

    return Tinf;
}

real get_stratified_tinf(face_t f, Thread *t, cell_t c0, Thread *t0)
{
    real Tinf;
    real xc[ND_ND];
    
    // Add input validation
    if (t == NULL) {
        Message("Warning: Thread is NULL in get_stratified_tinf, using default Tinf\n");
        return TBOT;
    }
    
    F_CENTROID(xc, f, t);
    
    // Try to get Tinf from UDM first (more efficient and consistent)
    if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL && 
        UDM_TINF >= 0 && UDM_TINF < N_UDM) {
        Tinf = C_UDMI(c0, t0, UDM_TINF);
        // Validate the UDM value - if it's unreasonable, recalculate
        if (Tinf < TBOT || Tinf > TTOP) {
            Tinf = compute_stratified_tinf(xc[1]); // Y-coordinate
        }
    } else {
        // Fallback to calculation if UDM not available
        Tinf = compute_stratified_tinf(xc[1]); // Y-coordinate
    }
    
    return Tinf;
}

real get_clamped_wall_temperature_and_film_temp(face_t f, Thread *t, real Tinf, real *Tf)
{
    real Tw;
    
    // Add input validation
    if (t == NULL || Tf == NULL) {
        Message("Error: Invalid inputs in get_clamped_wall_temperature_and_film_temp\n");
        if (Tf != NULL) *Tf = 300.0; // Default film temperature
        return 300.0; // Default wall temperature
    }
    
    Tw = F_T(f, t);
    
    if (Tw < TBOT) {
        Tw = TBOT;
    }
    else if (Tw > 1500.0) {
        Tw = 1500.0;
    }

    *Tf = 0.5 * (Tw + Tinf);
    
    return Tw;
}

real get_thermal_conductivity(real Tf)
{
    real T_C = Tf - 273.15;
    return 0.5718 + 0.00175 * T_C - 6.14e-6 * T_C * T_C;
}

real get_prandtl_number(real Tf)
{
    real T_C = Tf - 273.15;
    if (T_C < 0.0) T_C = 0.0;
    if (T_C > 100.0) T_C = 100.0;
    return 13.44 - 0.297 * T_C + 0.00243 * T_C * T_C - 7.14e-6 * T_C * T_C * T_C;
}

real get_viscosity(real Tf)
{
    real T_C = Tf - 273.15;
    if (T_C < 0.0) T_C = 0.0;
    if (T_C > 100.0) T_C = 100.0;
    return 0.001792 * exp(-0.0268 * T_C + 0.000102 * T_C * T_C);
}

real get_density(real Tf)
{
    real T_C = Tf - 273.15;
    if (T_C < 0.0) T_C = 0.0;
    if (T_C > 100.0) T_C = 100.0;
    return 999.84 - 0.0672 * T_C - 0.00796 * T_C * T_C + 2.53e-5 * T_C * T_C * T_C;
}

real get_thermal_expansion_coefficient(real Tf)
{
    real T_C = Tf - 273.15;
    if (T_C < 0.0) T_C = 0.0;
    if (T_C > 100.0) T_C = 100.0;
    return 2.14e-4 + 1.45e-6 * T_C + 6.79e-8 * T_C * T_C;
}

real get_friction_factor_goudar_sonnad(real Re, real epsilon, real D_h)
{
    if (Re < 4000.0) Re = 4000.0;
    if (epsilon <= 0.0) epsilon = 1e-6;
    if (D_h <= 0.0) D_h = 0.01;
    
    real a = 2.0 / log(10.0);
    real b = (epsilon / D_h) / 3.7;
    real d = log(10.0) * Re / 5.02;
    
    if (d <= 1.0) d = 2.0;
    
    real s = b * d + log(d);
    
    if (s <= 0.0) s = 0.1;
    
    real q = pow(s, s / (s + 1.0));
    
    if (q <= 0.0) q = 0.1;
    
    real g = b * d + log(d / q);
    
    if (q / g <= 0.0 || g <= 0.0) {
        return 0.3164 * pow(Re, -0.25);
    }
    
    real z = log(q / g);
    real DLA = z * g / (g + 1.0);
    real DCFA = DLA * (1.0 + (z / 2.0) / ((g + 1.0) * (g + 1.0) + (z / 3.0) * (2.0 * g - 1.0)));
    
    real inv_sqrt_f = a * (log(d / q) + DCFA);
    
    if (inv_sqrt_f <= 0.0) {
        return 0.3164 * pow(Re, -0.25);
    }
    
    real f = 1.0 / (inv_sqrt_f * inv_sqrt_f);
    
    if (f < 0.008) f = 0.008;
    if (f > 0.1) f = 0.1;
    
    return f;
}

real get_nusselt_churchill_chu_vertical(real Ra, real Pr)
{
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e12) Ra = 1e12;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu = pow(0.825 + (0.387 * pow(Ra, 1.0/6.0)) / pow(1.0 + pow(0.492/Pr, 9.0/16.0), 8.0/27.0), 2.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_nusselt_churchill_chu_horizontal(real Ra, real Pr)
{
    if (Ra < 1e-5) Ra = 1e-5;
    if (Ra > 1e12) Ra = 1e12;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu = pow(0.6 + (0.387 * pow(Ra, 1.0/6.0)) / pow(1.0 + pow(0.559/Pr, 9.0/16.0), 8.0/27.0), 2.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_nusselt_churchill_chu_elbow(real Ra, real Pr, real elbow_factor)
{
    if (elbow_factor < 0.0) elbow_factor = 0.0;
    if (elbow_factor > 1.0) elbow_factor = 1.0;
    
    real Nu_horizontal = get_nusselt_churchill_chu_horizontal(Ra, Pr);
    real Nu_vertical = get_nusselt_churchill_chu_vertical(Ra, Pr);
    real Nu = Nu_horizontal * (1.0 - elbow_factor) + Nu_vertical * elbow_factor;
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_nusselt_vertical_flat_plate(real Ra, real Pr)
{
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e9) Ra = 1e9;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu = pow(0.825 + (0.387 * pow(Ra, 1.0/6.0)) / pow(1.0 + pow(0.492/Pr, 9.0/16.0), 8.0/27.0), 2.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_nusselt_horizontal_flat_plate(real Ra, real Pr, int is_heated_up)
{
    if (Ra < 1e4) Ra = 1e4;
    if (Ra > 1e11) Ra = 1e11;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu;
    
    if (is_heated_up) {
        if (Ra < 1e7) {
            Nu = 0.54 * pow(Ra, 0.25);
        } else {
            Nu = 0.15 * pow(Ra, 1.0/3.0);
        }
    } else {
        Nu = 0.27 * pow(Ra, 0.25);
    }
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_nusselt_sphere(real Ra, real Pr)
{
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e11) Ra = 1e11;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu;
    
    Nu = 2.0 + (0.589 * pow(Ra, 0.25)) / pow(1.0 + pow(0.469/Pr, 9.0/16.0), 4.0/9.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;
    
    return Nu;
}

real get_finite_length_correction(real L_over_D)
{
    
    real correction_factor;
    
    if (L_over_D >= 35.0) {
        correction_factor = 1.0;
    } else if (L_over_D >= 20.0) {
        correction_factor = 0.85 + 0.15 * ((L_over_D - 20.0) / 15.0);
    } else if (L_over_D >= 10.0) {
        correction_factor = 0.75 + 0.10 * ((L_over_D - 10.0) / 10.0);
    } else {
        correction_factor = 0.75;
    }
    
    if (correction_factor < 0.5) correction_factor = 0.5;
    if (correction_factor > 1.0) correction_factor = 1.0;
    
    return correction_factor;
}

real get_nusselt_churchill_chu(real Ra, real Pr)
{
    return get_nusselt_churchill_chu_vertical(Ra, Pr);
}

real calculate_htc_vertical_cylinder_with_correction(real Tf, real Tw, real Tinf, real L_hydraulic, real L_actual)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real L_Ra = L_hydraulic;
    real L_h = L_hydraulic;
    
    real L_over_D = L_actual / L_hydraulic;
    
    real Gr = g_accel * beta * delta_T * pow(L_Ra, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e12) Ra = 1e12;

    real Nu = get_nusselt_churchill_chu_vertical(Ra, Pr);
    
    real correction_factor = get_finite_length_correction(L_over_D);
    Nu *= correction_factor;
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L_h;
    
    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_vertical_cylinder(real Tf, real Tw, real Tinf, real L)
{
    return calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L);
}

real calculate_htc_horizontal_cylinder(real Tf, real Tw, real Tinf, real D)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);

    real L = D;
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 1e-5) Ra = 1e-5;
    if (Ra > 1e12) Ra = 1e12;

    real Nu = get_nusselt_churchill_chu_horizontal(Ra, Pr);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L;

    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_elbow(real Tf, real Tw, real Tinf, real L, real elbow_factor)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);

    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 1e-5) Ra = 1e-5;
    if (Ra > 1e12) Ra = 1e12;

    real Nu = get_nusselt_churchill_chu_elbow(Ra, Pr, elbow_factor);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L;

    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_vertical_flat_plate(real Tf, real Tw, real Tinf, real L)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e9) Ra = 1e9;

    real Nu = get_nusselt_vertical_flat_plate(Ra, Pr);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_horizontal_flat_plate(real Tf, real Tw, real Tinf, real L, int is_heated_up)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 1e4) Ra = 1e4;
    if (Ra > 1e11) Ra = 1e11;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu;
    
    if (is_heated_up) {
        if (Ra < 1e7) {
            Nu = 0.54 * pow(Ra, 0.25);
        } else {
            Nu = 0.15 * pow(Ra, 1.0/3.0);
        }
    } else {
        Nu = 0.27 * pow(Ra, 0.25);
    }
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L;

    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_sphere(real Tf, real Tw, real Tinf, real D)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);

    real L = D;
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e11) Ra = 1e11;

    real Nu = get_nusselt_sphere(Ra, Pr);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 1000.0) Nu = 1000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 5000.0) h = 5000.0;

    return h;
}

real calculate_htc_spherical_cap(real Tf, real Tw, real Tinf, real D_hydraulic, real R_sphere, real height_cap, int cap_orientation)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real A_cap = 2.0 * M_PI * R_sphere * height_cap;
    
    real A_full_sphere = 4.0 * M_PI * R_sphere * R_sphere;

    real area_ratio = A_cap / A_full_sphere;
    
    if (area_ratio < 0.05) area_ratio = 0.05;
    if (area_ratio > 0.5) area_ratio = 0.5;
    
    real L = D_hydraulic;
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e11) Ra = 1e11;

    real Nu_base = get_nusselt_sphere(Ra, Pr);
    
    real geometric_factor = pow(area_ratio, 0.25);
    
    real orientation_factor;
    if (cap_orientation == 1) {
        orientation_factor = 0.8;
    } else {
        orientation_factor = 0.9;
    }
    
    real Nu = Nu_base * geometric_factor * orientation_factor;
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 500.0) Nu = 500.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 3000.0) h = 5000.0;

    return h;
}

real get_nusselt_vertical_plate_isothermal(real Ra, real Pr)
{
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e13) Ra = 1e13;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu;
    
    // Churchill and Chu correlation for vertical isothermal plates
    Nu = pow(0.825 + (0.387 * pow(Ra, 1.0/6.0)) / pow(1.0 + pow(0.492/Pr, 9.0/16.0), 8.0/27.0), 2.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;
    
    return Nu;
}

real get_nusselt_vertical_plate_uniform_heat_flux(real Ra, real Pr)
{
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e13) Ra = 1e13;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    
    real Nu;
    
    // Modified correlation for uniform heat flux boundary condition
    // Based on Sparrow and Gregg correlation with Churchill-Chu format
    Nu = 0.68 + (0.67 * pow(Ra, 0.25)) / pow(1.0 + pow(0.492/Pr, 9.0/16.0), 4.0/9.0);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;
    
    return Nu;
}

real get_nusselt_helical_coil(real Ra, real Pr, real D_coil, real D_tube, real pitch)
{
    if (Ra < 1.0) Ra = 1.0;
    if (Ra > 1e12) Ra = 1e12;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    if (D_coil <= 0.0) D_coil = 0.1;
    if (D_tube <= 0.0) D_tube = 0.01;
    if (pitch <= 0.0) pitch = D_tube;
    
    // Curvature ratio
    real delta = D_tube / D_coil;
    if (delta > 0.1) delta = 0.1; // Limit for validity
    
    // Pitch ratio
    real lambda = pitch / D_tube;
    if (lambda < 1.0) lambda = 1.0;
    if (lambda > 20.0) lambda = 20.0;
    
    // Base Nusselt number for straight tube (vertical orientation)
    real Nu_straight = get_nusselt_churchill_chu_vertical(Ra, Pr);
    
    // Curvature enhancement factor based on Dean number analogy
    real curvature_factor = 1.0 + 0.1 * pow(delta, 0.5) * pow(Ra, 0.1);
    
    // Pitch effect factor
    real pitch_factor = 1.0 + 0.05 * log(lambda);
    
    real Nu = Nu_straight * curvature_factor * pitch_factor;
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;
    
    return Nu;
}

real get_nusselt_spiral_coil(real Ra, real Pr, real D_coil, real D_tube)
{
    if (Ra < 1.0) Ra = 1.0;
    if (Ra > 1e12) Ra = 1e12;
    if (Pr < 0.1) Pr = 0.1;
    if (Pr > 100.0) Pr = 100.0;
    if (D_coil <= 0.0) D_coil = 0.1;
    if (D_tube <= 0.0) D_tube = 0.01;
    
    // Curvature ratio
    real delta = D_tube / D_coil;
    if (delta > 0.1) delta = 0.1; // Limit for validity
    
    // Base Nusselt number for straight tube
    real Nu_straight = get_nusselt_churchill_chu_horizontal(Ra, Pr);
    
    // Curvature enhancement factor for spiral coils
    real curvature_factor = 1.0 + 0.15 * pow(delta, 0.4) * pow(Ra, 0.08);
    
    real Nu = Nu_straight * curvature_factor;
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;
    
    return Nu;
}

real calculate_htc_vertical_plate_isothermal(real Tf, real Tw, real Tinf, real L)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e13) Ra = 1e13;

    real Nu = get_nusselt_vertical_plate_isothermal(Ra, Pr);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 10000.0) h = 10000.0;

    return h;
}

real calculate_htc_vertical_plate_uniform_heat_flux(real Tf, real Tw, real Tinf, real L)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 0.1) Ra = 0.1;
    if (Ra > 1e13) Ra = 1e13;

    real Nu = get_nusselt_vertical_plate_uniform_heat_flux(Ra, Pr);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 10000.0) h = 10000.0;

    return h;
}

real calculate_htc_helical_coil(real Tf, real Tw, real Tinf, real D_tube, real D_coil, real pitch)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real L = D_tube; // Characteristic length
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 1.0) Ra = 1.0;
    if (Ra > 1e12) Ra = 1e12;

    real Nu = get_nusselt_helical_coil(Ra, Pr, D_coil, D_tube, pitch);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 15000.0) h = 15000.0;

    return h;
}

real calculate_htc_spiral_coil(real Tf, real Tw, real Tinf, real D_tube, real D_coil)
{
    real k = get_thermal_conductivity(Tf);
    real Pr = get_prandtl_number(Tf);
    real mu = get_viscosity(Tf);
    real rho = get_density(Tf);

    real g_accel = 9.81;
    real beta = get_thermal_expansion_coefficient(Tf);
    real delta_T = fabs(Tw - Tinf);
    
    real L = D_tube; // Characteristic length
    
    real Gr = g_accel * beta * delta_T * pow(L, 3) / pow(mu / rho, 2);
    
    real Ra = Gr * Pr;
    
    if (Ra < 1.0) Ra = 1.0;
    if (Ra > 1e12) Ra = 1e12;

    real Nu = get_nusselt_spiral_coil(Ra, Pr, D_coil, D_tube);
    
    if (Nu < 2.0) Nu = 2.0;
    if (Nu > 2000.0) Nu = 2000.0;

    real h = Nu * k / L;
    
    if (h < 50.0) h = 50.0;
    if (h > 15000.0) h = 15000.0;

    return h;
}

// Template: Stratified temperature profile
// Usage: Apply to wall temperature boundary condition for bulk temperature
// Function: Sets linearly varying temperature based on Y-coordinate
DEFINE_PROFILE(Tinf_strat, t, i)
{
    // If UDF profiles are disabled, return immediately with default
    if (!g_udf_enabled) {
        face_t f;
        begin_f_loop(f, t) {
            F_PROFILE(f, t, i) = TBOT;
        }
        end_f_loop(f, t)
        return;
    }

    face_t f;
    real xc[ND_ND], y, Tinf;

    begin_f_loop(f, t)
    {
        // Use minimal checks - just set default if any basic operation fails
        F_CENTROID(xc, f, t);
        y = xc[1]; // Y-coordinate

        // Compute stratified Tinf using utility function
        Tinf = compute_stratified_tinf(y);

        F_PROFILE(f, t, i) = Tinf;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (first section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_first_vert, t, i)
{
    // If UDF profiles are disabled, return immediately with default
    if (!g_udf_enabled) {
        face_t f;
        begin_f_loop(f, t) {
            F_PROFILE(f, t, i) = DEFAULT_HTC_WATER;
        }
        end_f_loop(f, t)
        return;
    }

    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A1 / P1; // Hydraulic diameter as characteristic length
    const real L_actual = L1; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (first section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_first_vert_top, t, i)
{
    // If UDF profiles are disabled, return immediately with default
    if (!g_udf_enabled) {
        face_t f;
        begin_f_loop(f, t) {
            F_PROFILE(f, t, i) = DEFAULT_HTC_WATER;
        }
        end_f_loop(f, t)
        return;
    }

    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A11 / P11; // Hydraulic diameter as characteristic length
    const real L_actual = L11; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (first section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_first_vert_first_round, t, i)
{
    // If UDF profiles are disabled, return immediately with default
    if (!g_udf_enabled) {
        face_t f;
        begin_f_loop(f, t) {
            F_PROFILE(f, t, i) = DEFAULT_HTC_WATER;
        }
        end_f_loop(f, t)
        return;
    }

    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A12 / P12; // Hydraulic diameter as characteristic length
    const real L_actual = L12; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (first section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_first_vert_cone, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A13 / P13; // Hydraulic diameter as characteristic length
    const real L_actual = L13; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (first section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_first_vert_second_round, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A14 / P14; // Hydraulic diameter as characteristic length
    const real L_actual = L14; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (second section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_second_vert, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A2 / P2; // Hydraulic diameter as characteristic length
    const real L_actual = L2; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for horizontal cylinder
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for horizontal cylinder using Churchill-Bernstein correlation
DEFINE_PROFILE(htc_second_horz, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real D = 4 * A2 / P2; // Hydraulic diameter in meters

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using horizontal cylinder configuration
        h = calculate_htc_horizontal_cylinder(Tf, Tw, Tinf, D);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for elbow (mixed orientation)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for elbow using combined horizontal/vertical correlation
DEFINE_PROFILE(htc_second_elbow, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A2 / P2; // Characteristic length (hydraulic diameter) in meters
    const real elbow_factor = 0.5; // 0.5 = 45° elbow (half horizontal, half vertical)

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using elbow configuration (mixed horizontal/vertical)
        h = calculate_htc_elbow(Tf, Tw, Tinf, L, elbow_factor);

        // Store HTC in adjacent fluid cell for visualization (with safety comments)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for vertical cylinder (third section)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for vertical cylinder using Churchill-Chu correlation
DEFINE_PROFILE(htc_third_vert, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A3 / P3; // Hydraulic diameter as characteristic length
    const real L_actual = L3; // Actual tube length for finite length correction

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using vertical cylinder configuration with finite length correction
        h = calculate_htc_vertical_cylinder_with_correction(Tf, Tw, Tinf, L, L_actual);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for horizontal flat plate (diverter)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for heated plate facing down (stable configuration)
DEFINE_PROFILE(htc_diverter, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A3 / P3; // Hydraulic diameter in meters

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using horizontal flat plate configuration (heated down - stable)
        h = calculate_htc_horizontal_flat_plate(Tf, Tw, Tinf, L, 0);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Template: Heat transfer coefficient for horizontal flat plate (collector)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for heated plate facing up (unstable configuration)
DEFINE_PROFILE(htc_collector_plate, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real L = 4 * A3 / P3; // Hydraulic diameter in meters

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using horizontal flat plate configuration (heated up - unstable)
        h = calculate_htc_horizontal_flat_plate(Tf, Tw, Tinf, L, 1);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}


// Template: Heat transfer coefficient for horizontal flat plate (endcap)
// Usage: Apply to wall boundary condition for convective heat transfer coefficient
// Function: Calculates HTC for heated plate facing down (stable configuration)
DEFINE_PROFILE(htc_endcap, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real Tw, Tinf, Tf, h;
    const real D = 4 * A1 / P1; // Hydraulic diameter from first pass (A1/P1)
    const real R_sphere = END_CAP_ROC; // Use defined radius of curvature for endcap
    
    // Calculate cap height using utility function
    real height_cap = calculate_spherical_cap_height(D, R_sphere);
    
    // Endcap orientation: 1 = bottom cap, 0 = top cap
    // Set this based on your actual geometry
    const int cap_orientation = 1;  // Assuming bottom cap (change if needed)

    begin_f_loop(f, t)
    {
        // Check if flow is properly initialized
        if (!check_flow_initialization(f, t, &c0, &t0)) {
            set_default_profile_value(f, t, i, c0, t0, UDM_HTC);
            continue;
        }

        // Get bulk temperature using utility function
        Tinf = get_stratified_tinf(f, t, c0, t0);
        
        // Get clamped wall temperature and calculate film temperature
        Tw = get_clamped_wall_temperature_and_film_temp(f, t, Tinf, &Tf);
        
        // Calculate HTC using spherical cap correlation with geometric corrections
        h = calculate_htc_spherical_cap(Tf, Tw, Tinf, D, R_sphere, height_cap, cap_orientation);

        // Store HTC in adjacent fluid cell for visualization (with safety checks)
        if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
            C_UDMI(c0, t0, UDM_HTC) = h;
        }

        F_PROFILE(f, t, i) = h;
    }
    end_f_loop(f, t)
}

// Function to compute mass-averaged cp over a boundary
real compute_mass_averaged_cp(Thread *t)
{
    face_t f;
    cell_t c0;
    Thread *t0;
    real sum_mdot_cp = 0.0;
    real sum_mdot = 0.0;
    real mdot, cp_cell;
    
    begin_f_loop(f, t)
    {
        c0 = F_C0(f, t); // Adjacent cell
        t0 = THREAD_T0(t);
        mdot = MASS_FLOW_RATE;
        cp_cell = C_CP(c0, t0); // Cell cp
        sum_mdot_cp += mdot * cp_cell;
        sum_mdot += mdot;
    }
    end_f_loop(f, t)
    
    if (sum_mdot > 0.0)
        return sum_mdot_cp / sum_mdot;
    else
        return 0.0;
}

// Empirical radiative loss from wall
// DEFINE_PROFILE(empirical_radiative_wall_loss, t, i)
// {
//     face_t f;
//     real Tw, q_rad;
//     begin_f_loop(f, t)
//     {
//         Tw = F_T(f, t);
//         q_rad = EMISSIVITY_WALL * SIGMA_SB * (pow(Tw, 4) - pow(T_AMB, 4));
//         F_PROFILE(f, t, i) = q_rad;
//         F_UDMI(f, t, UDM_RAD_WALL_LOSS) = q_rad; // Store in UDM
//     }
//     end_f_loop(f, t)
// }

// // Empirical radiative loss from flue gas
// DEFINE_PROFILE(empirical_radiative_gas_loss, t, i)
// {
//     face_t f;
//     real Tg, q_rad_gas;
//     begin_f_loop(f, t)
//     {
//         Tg = F_T(f, t);
//         q_rad_gas = EMISSIVITY_GAS * SIGMA_SB * (pow(Tg, 4) - pow(T_AMB, 4));
//         F_PROFILE(f, t, i) = q_rad_gas;
//         F_UDMI(f, t, UDM_RAD_GAS_LOSS) = q_rad_gas; // Store in UDM
//     }
//     end_f_loop(f, t)
// }

// // ASME/Perry's percent stack loss percentage
// DEFINE_PROFILE(stack_loss_percent, t, i)
// {
//     face_t f;
//     real Tg, cp_massavg, percent_loss;
//     cp_massavg = compute_mass_averaged_cp(t); // Mass-averaged cp for boundary
//     begin_f_loop(f, t)
//     {
//         Tg = F_T(f, t);
//         // ASME/Perry's: percent loss = cp*(Tg-Tamb)/(cp*(Tflame-Tamb)) * 100
//         percent_loss = (Tg - T_AMB) / (FLAME_TEMP - T_AMB) * 100.0;
//         F_PROFILE(f, t, i) = percent_loss;
//         F_UDMI(f, t, UDM_STACK_LOSS_PERCENT) = percent_loss; // Store in UDM
//     }
//     end_f_loop(f, t)
// }

// DEFINE_PROFILE(stack_loss_actual, t, i)
// {
//     face_t f;
//     real Tg, mdot, cp_massavg, percent_loss, q_stack_actual;
//     cp_massavg = compute_mass_averaged_cp(t); // Mass-averaged cp for boundary
//     begin_f_loop(f, t)
//     {
//         Tg = F_T(f, t);
//         mdot = MASS_FLOW_RATE;
//         // ASME/Perry's percent loss
//         percent_loss = F_UDMI(f, t, UDM_STACK_LOSS_PERCENT);
//         // Actual stack loss (W)
//         q_stack_actual = percent_loss * mdot * cp_massavg * (FLAME_TEMP - T_AMB);
//         F_UDMI(f, t, UDM_STACK_LOSS) = q_stack_actual;
//         F_PROFILE(f, t, i) = q_stack_actual;
//     }
//     end_f_loop(f, t)
// }

// UDF to detect condensation and calculate condensed mass per cell
// TEMPORARILY DISABLED - MAY CAUSE SEGFAULT DURING INITIALIZATION
/*
DEFINE_EXECUTE_AT_END(condensation_detection)
{
    Domain *d = Get_Domain(1);
    Thread *t;
    cell_t c;
    real T_cell, m_cell, m_condensed;

    // Add domain validation
    if (d == NULL) {
        Message("Error: Domain 1 not found in condensation_detection\n");
        return;
    }
    
    // Validate UDM setup before proceeding
    if (!validate_udm_setup(d)) {
        Message("Error: UDM validation failed in condensation_detection\n");
        return;
    }

    thread_loop_c(t, d)
    {
        if (FLUID_THREAD_P(t)) // Correct check for fluid cell threads
        {
            begin_c_loop(c, t)
            {
                // Add safety checks for UDM storage
                if (THREAD_STORAGE(t, SV_UDM_I) == NULL) continue;
                
                T_cell = C_T(c, t);
                m_cell = C_VOLUME(c, t) * C_R(c, t); // Cell mass (kg)
                m_condensed = 0.0;
                if (T_cell < DEW_POINT)
                {
                    m_condensed = m_cell * WATER_VAPOR_MASS_FRAC;
                }
                C_UDMI(c, t, UDM_CONDENSED_MASS) = m_condensed;
            }
            end_c_loop(c, t)
        }
    }
}
*/

// --- Modified sum_total_heat_transfer_rate to include wall and gas radiative losses ---
// DEFINE_ON_DEMAND(sum_total_heat_transfer_rate)
// {
//     Domain *d = Get_Domain(1);
//     Thread *t;
//     face_t f;
//     cell_t c0;
//     Thread *t0;
//     real total_heat_transfer = 0.0;
//     real total_stack_loss = 0.0;
//     real total_rad_wall_loss = 0.0;
//     real total_rad_gas_loss = 0.0;

//     // Add domain validation
//     if (d == NULL) {
//         Message("Error: Domain 1 not found in sum_total_heat_transfer_rate\n");
//         return;
//     }

//     thread_loop_f(t, d)
//     {
//         if (THREAD_TYPE(t) == THREAD_F_WALL || THREAD_TYPE(t) == THREAD_F_OUTFLOW || THREAD_TYPE(t) == THREAD_F_POUTLET)
//         {
//             begin_f_loop(f, t)
//             {
//                 // Get adjacent cell for UDM access - face UDM is not typically allocated
//                 c0 = F_C0(f, t);
//                 t0 = THREAD_T0(t);
                
//                 // Safety checks for cell and UDM access
//                 if (t0 != NULL && c0 >= 0 && THREAD_STORAGE(t0, SV_UDM_I) != NULL) {
//                     total_stack_loss += C_UDMI(c0, t0, UDM_STACK_LOSS);
//                     total_rad_wall_loss += C_UDMI(c0, t0, UDM_RAD_WALL_LOSS);
//                     total_rad_gas_loss += C_UDMI(c0, t0, UDM_RAD_GAS_LOSS);
//                 }
//             }
//             end_f_loop(f, t)
//         }
//     }
//     total_heat_transfer = total_stack_loss + total_rad_wall_loss + total_rad_gas_loss;
//     Message("Total Stack Loss (W): %g\n", total_stack_loss);
//     Message("Total Wall Radiative Loss (W): %g\n", total_rad_wall_loss);
//     Message("Total Gas Radiative Loss (W): %g\n", total_rad_gas_loss);
//     Message("Total Heat Loss (Stack + Wall Rad + Gas Rad, W): %g\n", total_heat_transfer);
// }

// // --- On-demand function to sum total condensed mass in the domain ---
// DEFINE_ON_DEMAND(sum_total_condensed_mass)
// {
//     Domain *d = Get_Domain(1);
//     Thread *t;
//     cell_t c;
//     real total_condensed_mass = 0.0;

//     // Add domain validation
//     if (d == NULL) {
//         Message("Error: Domain 1 not found in sum_total_condensed_mass\n");
//         return;
//     }

//     thread_loop_c(t, d)
//     {
//         if (FLUID_THREAD_P(t))
//         {
//             begin_c_loop(c, t)
//             {
//                 // Add safety checks for UDM storage
//                 if (THREAD_STORAGE(t, SV_UDM_I) != NULL) {
//                     total_condensed_mass += C_UDMI(c, t, UDM_CONDENSED_MASS);
//                 }
//             }
//             end_c_loop(c, t)
//         }
//     }
//     Message("Total Condensed Mass (kg): %g\n", total_condensed_mass);
// }

// // --- On-demand function to compute latent heat lost from uncondensed vapor ---
// DEFINE_ON_DEMAND(sum_latent_heat_uncondensed)
// {
//     Domain *d = Get_Domain(1);
//     Thread *t;
//     cell_t c;
//     real total_water_vapor = 0.0;
//     real total_condensed_mass = 0.0;
//     real latent_heat_lost = 0.0;

//     // Add domain validation
//     if (d == NULL) {
//         Message("Error: Domain 1 not found in sum_latent_heat_uncondensed\n");
//         return;
//     }

//     thread_loop_c(t, d)
//     {
//         if (FLUID_THREAD_P(t))
//         {
//             begin_c_loop(c, t)
//             {
//                 total_water_vapor += C_VOLUME(c, t) * C_R(c, t) * WATER_VAPOR_MASS_FRAC;
//                 // Add safety checks for UDM storage
//                 if (THREAD_STORAGE(t, SV_UDM_I) != NULL) {
//                     total_condensed_mass += C_UDMI(c, t, UDM_CONDENSED_MASS);
//                 }
//             }
//             end_c_loop(c, t)
//         }
//     }
//     latent_heat_lost = (total_water_vapor - total_condensed_mass) * LATENT_HEAT_VAP;
//     Message("Latent Heat Lost from Uncondensed Vapor (J): %g\n", latent_heat_lost);
// }

// --- Area-weighted flue gas mass flow profile ---
DEFINE_PROFILE(flue_gas_mass_flow, t, i)
{
    face_t f;
    // real total_area = 0.0;
    // real face_area[ND_ND];

    // First loop: compute total inlet area
    // begin_f_loop(f, t)
    // {
    //     F_AREA(face_area, f, t);
    //     total_area += NV_MAG(face_area);
    // }
    // end_f_loop(f, t)

    // Second loop: assign area-weighted mass flow to each face
    begin_f_loop(f, t)
    {
        // F_AREA(face_area, f, t);
        // F_PROFILE(f, t, i) = MASS_FLOW_RATE * NV_MAG(face_area) / total_area;
        F_PROFILE(f, t, i) = MASS_FLOW_RATE;
    }
    end_f_loop(f, t)
    // Message("Flue gas mass flow rate set to %g kg/s on boundary %d\n", MASS_FLOW_RATE, THREAD_ID(t));
}

// // --- Convert radiative loss from multiple boundaries to negative volumetric heat generation rate ---
// DEFINE_PROFILE(radiative_loss_to_volumetric_source, d)
// {
//     Thread *t_boundary;
//     Thread *t_cell;
//     face_t f;
//     cell_t c0;
//     real total_rad_loss = 0.0;
//     real total_volume = 0.0;
//     real rad_loss_per_volume = 0.0;
//     int i;
    
//     // Loop over all radiative loss boundaries
//     for (i = 0; i < NUM_RAD_LOSS_BOUNDARIES; ++i)
//     {
//         t_boundary = Lookup_Thread(d, RAD_LOSS_BOUNDARY_IDS[i]);
//         begin_f_loop(f, t_boundary)
//         {
//             c0 = F_C0(f, t_boundary);
//             t_cell = THREAD_T0(t_boundary);
//             if (t_cell && c0 >= 0)
//             {
//                 // Assume radiative loss stored in UDM_RAD_WALL_LOSS per face
//                 total_rad_loss += C_UDMI(c0, t_cell, UDM_RAD_WALL_LOSS) * F_AREA(f, t_boundary);
//                 total_volume += C_VOLUME(c0, t_cell);
//             }
//         }
//         end_f_loop(f, t_boundary)
//     }
    
//     if (total_volume > 0.0)
//         rad_loss_per_volume = -total_rad_loss / total_volume;
//     else
//         rad_loss_per_volume = 0.0;
    
//     // Store as negative heat generation rate in a UDM for each cell
//     thread_loop_c(t_cell, d)
//     {
//         if (THREAD_TYPE(t_cell) == THREAD_CELL)
//         {
//             begin_c_loop(c0, t_cell)
//             {
//                 C_UDMI(c0, t_cell, UDM_RAD_WALL_LOSS) = rad_loss_per_volume;
//             }
//             end_c_loop(c0, t_cell)
//         }
//     }
//     Message("Radiative loss (all boundaries) converted to volumetric source: %g W/m^3\n", rad_loss_per_volume);
// }

DEFINE_PROFILE(radiative_wall_loss_volumetric, t, i)
{
    face_t f;
    real Tw, T_amb, q_rad, q_vol;
    T_amb = T_AMB; // Ambient temperature (K)
    
    begin_f_loop(f, t)
    {
        Tw = F_T(f, t); // Wall temperature
        // Radiative heat flux (W/m^2): q_rad = emissivity * sigma * (Tw^4 - T_amb^4)
        q_rad = EMISSIVITY_WALL * SIGMA_SB * (pow(Tw, 4) - pow(T_amb, 4));
        // Convert to volumetric source (W/m^3): q_vol = q_rad / wall_thickness
        q_vol = q_rad / WALL_THICKNESS;
        F_PROFILE(f, t, i) = -q_vol;
    }
    end_f_loop(f, t)
}


// Utility: Area average of a variable over a boundary
real area_average_on_boundary(Domain *d, int zone_id, real (*face_var)(face_t, Thread *)) {
    Thread *t = Lookup_Thread(d, zone_id);
    face_t f;
    real sum = 0.0, area = 0.0;
    real A[ND_ND];
    begin_f_loop(f, t) {
        F_AREA(A, f, t);
        sum += face_var(f, t) * NV_MAG(A);
        area += NV_MAG(A);
    } end_f_loop(f, t);
    return (area > 0.0) ? sum / area : 0.0;
}

real total_heat_transfer_on_boundary(Domain *d, int zone_id) {
    Thread *t = Lookup_Thread(d, zone_id);
    face_t f;
    real total = 0.0;
    real A[ND_ND];
    begin_f_loop(f, t) {
        F_AREA(A, f, t);
        total += F_FLUX(f, t) * NV_MAG(A);
    } end_f_loop(f, t);
    return total; // W
}

// Utility: Volume average of a cell variable over the domain
real volume_average_in_domain(Domain *d, real (*cell_var)(cell_t, Thread *)) {
    Thread *t;
    cell_t c;
    real sum = 0.0, vol = 0.0;
    thread_loop_c(t, d) {
        begin_c_loop(c, t) {
            sum += cell_var(c, t) * C_VOLUME(c, t);
            vol += C_VOLUME(c, t);
        } end_c_loop(c, t);
    }
    return (vol > 0.0) ? sum / vol : 0.0;
}

// Utility: Mass average of a cell variable over the domain
real mass_average_in_domain(Domain *d, real (*cell_var)(cell_t, Thread *)) {
    Thread *t;
    cell_t c;
    real sum = 0.0, mass = 0.0;
    thread_loop_c(t, d) {
        begin_c_loop(c, t) {
            sum += cell_var(c, t) * C_VOLUME(c, t) * C_R(c, t);
            mass += C_VOLUME(c, t) * C_R(c, t);
        } end_c_loop(c, t);
    }
    return (mass > 0.0) ? sum / mass : 0.0;
}

// Macros for cell/face variables
real cell_temp(cell_t c, Thread *t) { return C_T(c, t); }
real cell_press(cell_t c, Thread *t) { return C_P(c, t); }
real cell_dens(cell_t c, Thread *t) { return C_R(c, t); }
// real cell_yi0(cell_t c, Thread *t) { return C_YI(c, t, 0); } // Example: species 0
// real cell_yi1(cell_t c, Thread *t) { return C_YI(c, t, 1); } // Example: species 1
real face_temp(face_t f, Thread *t) { return F_T(f, t); }
real face_heat_flux(face_t f, Thread *t) { return F_FLUX(f, t); }

typedef struct {
    int id;
    const char *name;
} ZoneInfo;

ZoneInfo zone_map[] = {
    {137, "collector_plate_walls_air"},
    {138, "collector_plate_walls_steel"},
    {136, "collector_plate_walls_water"},
    {143, "diverter_walls_air"},
    {133, "diverter_walls_water"},
    {140, "first_pass_burner_walls_steel"},
    {139, "first_pass_end_cap_walls_steel"},
    {135, "first_pass_end_cap_walls_water"},
    {147, "first_pass_vertical_walls_cone_water"},
    {146, "first_pass_vertical_walls_first_round_water"},
    {148, "first_pass_vertical_walls_second_round_water"},
    {145, "first_pass_vertical_walls_top_water"},
    {134, "first_pass_vertical_walls_water"},
    {126, "internal_baffles.1"},
    {1, "internal_baffles.1-shadow"},
    {127, "internal_baffles.2"},
    {2, "internal_baffles.2-shadow"},
    {130, "second_pass_elbow_walls_water"},
    {131, "second_pass_horizontal_walls_water"},
    {129, "second_pass_vertical_walls_water"},
    {132, "second_weld_walls_water"},
    {128, "third_pass_vertical_walls_water"}
    // Add more as needed
};
int zone_map_size = sizeof(zone_map)/sizeof(zone_map[0]);

const char* get_zone_name(int id) {
    for (int i = 0; i < zone_map_size; ++i) {
        if (zone_map[i].id == id) return zone_map[i].name;
    }
    return "Unknown";
}

DEFINE_ON_DEMAND(write_markdown_summary)
{
    FILE *fp;
    Domain *d = Get_Domain(1);
    fp = fopen("fluent_summary_report.md", "w");
    if (fp == NULL) return;

    fprintf(fp, "# Fluent CFD Analysis Summary\n\n");
    fprintf(fp, "## Domain Averages\n");
    fprintf(fp, "- Volume Average Temperature: %.2f K\n", volume_average_in_domain(d, cell_temp));
    fprintf(fp, "- Mass Average Temperature: %.2f K\n", mass_average_in_domain(d, cell_temp));
    fprintf(fp, "- Volume Average Pressure: %.2f Pa\n", volume_average_in_domain(d, cell_press));
    fprintf(fp, "- Volume Average Density: %.4f kg/m³\n", volume_average_in_domain(d, cell_dens));
    // fprintf(fp, "- Species 0 Mass Fraction (avg): %.4f\n", volume_average_in_domain(d, cell_yi0));
    // fprintf(fp, "- Species 1 Mass Fraction (avg): %.4f\n", volume_average_in_domain(d, cell_yi1));
    // Add more species as needed

    // Boundary analysis
    fprintf(fp, "\n## Boundary Heat Transfer\n");
    int boundary_ids[] = {137, 136, 143, 133, 135, 147, 146, 148, 145, 134, 130, 131, 129, 132, 128}; // Example boundary zone IDs
    int n_boundaries = sizeof(boundary_ids)/sizeof(boundary_ids[0]);
    real total_heat_transfer = 0.0;
    real boundary_heat[n_boundaries];
    for (int i = 0; i < n_boundaries; ++i) {
        boundary_heat[i] = total_heat_transfer_on_boundary(d, boundary_ids[i]);
        total_heat_transfer += boundary_heat[i];
    }
    for (int i = 0; i < n_boundaries; ++i) {
        fprintf(fp, "### Boundary %d (%s)\n", boundary_ids[i], get_zone_name(boundary_ids[i]));
        fprintf(fp, "- Wall Temperature (area avg): %.2f K\n", area_average_on_boundary(d, boundary_ids[i], face_temp));
        fprintf(fp, "- Wall Heat Flux (area avg): %.2f W/m²\n", area_average_on_boundary(d, boundary_ids[i], face_heat_flux));
        fprintf(fp, "- Total Heat Transfer: %.2f W\n", boundary_heat[i]);
        fprintf(fp, "- %% of Total Heat Transfer: %.2f%%\n", (total_heat_transfer > 0.0) ? 100.0 * boundary_heat[i] / total_heat_transfer : 0.0);
    }
    fprintf(fp, "\n## Overall Results\n");
    fprintf(fp, "- Total Heat Transfer (all boundaries): %.2f W\n", total_heat_transfer);
    // Example: Losses and efficiency (replace with actual calculation)
    real input_power = MASS_FLOW_RATE * CP_GAS * (FLAME_TEMP - T_AMB); // Example input power
    real efficiency = (total_heat_transfer > 0.0 && input_power > 0.0) ? 100.0 * total_heat_transfer / input_power : 0.0;
    fprintf(fp, "- Input Power (approx): %.2f W\n", input_power);
    fprintf(fp, "- Efficiency: %.2f%%\n", efficiency);
    fprintf(fp, "- Losses: %.2f W\n", input_power - total_heat_transfer);

    fclose(fp);
}

// UDF initialization checker - call this first before other UDFs
DEFINE_ON_DEMAND(check_udf_initialization)
{
    Domain *d = Get_Domain(1);
    
    Message("=== UDF INITIALIZATION CHECK ===\n");
    
    if (d == NULL) {
        Message("ERROR: Domain 1 not found!\n");
        Message("Solution: Load a valid case file first\n");
        return;
    }
    
    if (!validate_udm_setup(d)) {
        Message("ERROR: UDM setup validation failed!\n");
        Message("Solution: Allocate at least %d UDM variables in Solution > Solution Methods > User-Defined Memory\n", 
                UDM_CONDENSED_MASS + 1);
        return;
    }
    
    Message("SUCCESS: UDF initialization checks passed\n");
    Message("You can now safely use the UDF profiles and execute functions\n");
    Message("================================\n");
}

// Enable UDF profiles after solution initialization
DEFINE_ON_DEMAND(enable_udf_profiles)
{
    g_udf_enabled = 1;
    Message("=== UDF PROFILES ENABLED ===\n");
    Message("UDF profiles are now active and will calculate HTC values\n");
    Message("Call this AFTER solution initialization is complete\n");
    Message("============================\n");
}

// Disable UDF profiles (use before solution initialization)
DEFINE_ON_DEMAND(disable_udf_profiles)
{
    g_udf_enabled = 0;
    Message("=== UDF PROFILES DISABLED ===\n");
    Message("UDF profiles will return default values only\n");
    Message("Use this BEFORE solution initialization\n");
    Message("==============================\n");
}