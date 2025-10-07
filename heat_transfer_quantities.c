#include "udf.h"
#include <stdio.h>
#include <string.h>

#define NUM_BINS 100
#define WALL_ZONE_NAME "water_wall_group"  // Replace with your wall zone name

// Utility function to write a profile to CSV
void write_profile_to_csv(const char *filename, real values[NUM_BINS], real time, int write_time)
{
    FILE *fp = fopen(filename, write_time ? "a" : "w");

    if (fp == NULL) return;

    static int header_written = 0;
    if (!header_written)
    {
        if (write_time) fprintf(fp, "t");
        for (int i = 0; i < NUM_BINS; i++)
            fprintf(fp, ",y%d", i + 1);
        fprintf(fp, "\n");
        header_written = 1;
    }

    if (write_time) fprintf(fp, "%.6f", time);
    for (int i = 0; i < NUM_BINS; i++)
        fprintf(fp, ",%.6f", values[i]);
    fprintf(fp, "\n");

    fclose(fp);
}

// Core function to compute and write a field
void compute_wall_profile(Domain *d, const char *field, const char *filename, int write_time)
{
    Thread *t;
    face_t f;

    real y_min = 1e9, y_max = -1e9;
    real bin_sum[NUM_BINS] = {0.0};
    real bin_area[NUM_BINS] = {0.0}; // For area-weighted averaging

    thread_loop_f(t, d)
    {
        if (strcmp(THREAD_NAME(t), WALL_ZONE_NAME) == 0)
        {
            // First pass: find Y range
            begin_f_loop(f, t)
            {
                real y = F_Y(f, t);
                if (y < y_min) y_min = y;
                if (y > y_max) y_max = y;
            }
            end_f_loop(f, t)

            real bin_height = (y_max - y_min) / NUM_BINS;

            // Second pass: compute field (area-weighted)
            begin_f_loop(f, t)
            {
                real y = F_Y(f, t);
                real area = F_AREA(f, t);
                int bin = (int)((y - y_min) / bin_height);
                if (bin < 0) bin = 0;
                if (bin >= NUM_BINS) bin = NUM_BINS - 1;

                real value = 0.0;

                if (strcmp(field, "heat_flux") == 0)
                    value = F_HEAT_FLUX(f, t);
                else if (strcmp(field, "wall_temp") == 0)
                    value = F_T(f, t);
                else if (strcmp(field, "adj_temp") == 0)
                    value = F_T_ADJ(f, t);
                else if (strcmp(field, "tinf") == 0)
                    value = F_T_INF(f, t);

                bin_sum[bin] += value * area;
                bin_area[bin] += area;
            }
            end_f_loop(f, t)
        }
    }

    // Average and write
    real avg[NUM_BINS];
    for (int i = 0; i < NUM_BINS; i++)
        avg[i] = (bin_area[i] > 0.0) ? bin_sum[i] / bin_area[i] : 0.0;

    write_profile_to_csv(filename, avg, CURRENT_TIME, write_time);
}

// One UDF per field
DEFINE_EXECUTE_AT_END(log_heat_flux_profile)
{
    compute_wall_profile(Get_Domain(1), "heat_flux", "heat_flux_profile.csv", 1);
}

DEFINE_EXECUTE_AT_END(log_wall_temp_profile)
{
    compute_wall_profile(Get_Domain(1), "wall_temp", "wall_temp_profile.csv", 1);
}

DEFINE_EXECUTE_AT_END(log_adj_temp_profile)
{
    compute_wall_profile(Get_Domain(1), "adj_temp", "adjacent_temp_profile.csv", 1);
}

DEFINE_EXECUTE_AT_END(log_tinf_profile)
{
    compute_wall_profile(Get_Domain(1), "tinf", "tinf_profile.csv", 1);
}
