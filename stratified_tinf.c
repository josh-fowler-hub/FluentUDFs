#include "udf.h"

#define TBOT 294.26   // Bottom temperature in Kelvin
#define TTOP 302.59   // Top temperature in Kelvin
#define Y_MIN 0.25    // Minimum Y-coordinate for this zone, meters
#define Y_MAX 1.25    // Maximum Y-coordinate for this zone, meters

DEFINE_PROFILE(stratified_tinf, t, i)
{
    face_t f;
    real xc[ND_ND], y, Tinf;

    begin_f_loop(f, t)
    {
        F_CENTROID(xc, f, t);
        y = xc[1]; // Y-coordinate

        // Clamp y between Y_MIN and Y_MAX
        if (y < Y_MIN)
            y = Y_MIN;
        else if (y > Y_MAX)
            y = Y_MAX;

        // Compute normalized height and Tinf
        if ((Y_MAX - Y_MIN) > 1e-12)
            Tinf = TBOT + (TTOP - TBOT) * ((y - Y_MIN) / (Y_MAX - Y_MIN));
        else
            Tinf = TBOT;

        F_PROFILE(f, t, i) = Tinf;
    }
    end_f_loop(f, t)
}
