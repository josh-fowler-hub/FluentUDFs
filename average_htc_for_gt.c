#include "udf.h"
#include <stdio.h>
#include <math.h> // Added for fabs

#define NUM_BINS 100
#define WALL_ZONE_NAME "water_wall_group"  // Replace with your named wall zone
#define OUTPUT_FILE "htc_profile.csv"

DEFINE_EXECUTE_AT_END(log_htc_profile)
{
    Domain *d = Get_Domain(1);
    Thread *t;
    face_t f;

    real y_min = 1e9, y_max = -1e9;
    real bin_sum[NUM_BINS] = {0.0};
    real bin_area[NUM_BINS] = {0.0}; // For area-weighted averaging
    static int header_written = 0;

    // Find the thread by name
    thread_loop_f(t, d)
    {
        if (strcmp(THREAD_NAME(t), WALL_ZONE_NAME) == 0)
        {
            // First pass: find y_min and y_max
            begin_f_loop(f, t)
            {
                real y = F_Y(f, t);
                if (y < y_min) y_min = y;
                if (y > y_max) y_max = y;
            }
            end_f_loop(f, t)

            real bin_height = (y_max - y_min) / NUM_BINS;

            // Second pass: bin HTC values (area-weighted)
            begin_f_loop(f, t)
            {
                real y = F_Y(f, t);
                real area = F_AREA(f, t);
                real q = F_HEAT_FLUX(f, t);
                real Twall = F_T(f, t);
                real Tinf = F_T_INF(f, t);
                real htc = (Twall != Tinf) ? fabs(q) / fabs(Twall - Tinf) : 0.0;

                int bin = (int)((y - y_min) / bin_height);
                if (bin < 0) bin = 0;
                if (bin >= NUM_BINS) bin = NUM_BINS - 1; // Clamp to valid range

                bin_sum[bin] += htc * area;
                bin_area[bin] += area;
            }
            end_f_loop(f, t)

            // Write to CSV
            FILE *fp = fopen(OUTPUT_FILE, header_written ? "a" : "w");
            if (!fp) return; // File open failed, exit safely

            if (!header_written)
            {
                fprintf(fp, "t");
                for (int i = 0; i < NUM_BINS; i++)
                    fprintf(fp, ",y%d", i + 1);
                fprintf(fp, "\n");
                header_written = 1;
            }

            fprintf(fp, "%.6f", CURRENT_TIME);
            for (int i = 0; i < NUM_BINS; i++)
            {
                real avg_htc = (bin_area[i] > 0.0) ? bin_sum[i] / bin_area[i] : 0.0;
                fprintf(fp, ",%.4f", avg_htc);
            }
            fprintf(fp, "\n");

            fclose(fp);
        }
    }
}
