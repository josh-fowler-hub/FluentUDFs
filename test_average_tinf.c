#include "udf.h"
#include <stdio.h>
#include <math.h> // Include for fabs

#define MAX_ROWS 1000
#define NUM_TC 6

static const real max_tc_height = 45.0;
static const real tc_heights[NUM_TC] = {3.797244094/max_tc_height, 11.39173228/max_tc_height, 18.98622047/max_tc_height, 26.58070866/max_tc_height, 34.17519685/max_tc_height, 41.76968504/max_tc_height};
static real avg_tc[NUM_TC];
static int data_loaded = 0;

void load_and_average_csv()
{
    FILE *fp = fopen("thermocouple_data.csv", "r");
    if (!fp) return;

    char line[256];
    int row = 0;
    fgets(line, sizeof(line), fp); // skip header

    for (int i = 0; i < NUM_TC; i++) avg_tc[i] = 0.0;

    while (fgets(line, sizeof(line), fp) && row < MAX_ROWS)
    {
        double t, vals[NUM_TC];
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &t, &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5]) == 7)
        {
            for (int i = 0; i < NUM_TC; i++)
                avg_tc[i] += vals[i];
            row++;
        }
    }

    if (row > 0)
        for (int i = 0; i < NUM_TC; i++)
            avg_tc[i] /= row;

    fclose(fp); // Close the file after reading
    data_loaded = 1;
}

DEFINE_PROFILE(static_tinf_profile, t, i)
{
    face_t f;
    real a = 0.0, b = 0.0;

    if (!data_loaded)
        load_and_average_csv();

    // Linear least squares fit: T = a*y + b
    real sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (int i = 0; i < NUM_TC; i++)
    {
        real x = tc_heights[i];
        real y = avg_tc[i];
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    real denom = NUM_TC * sum_x2 - sum_x * sum_x;
    if (fabs(denom) > 1e-12) // Avoid division by zero
    {
        a = (NUM_TC * sum_xy - sum_x * sum_y) / denom;
        b = (sum_y * sum_x2 - sum_x * sum_xy) / denom;
    }
    else
    {
        a = 0.0;
        b = sum_y / NUM_TC;
    }

    begin_f_loop(f, t)
    {
        real y = F_Y(f, t);
        if (y <= tc_heights[0])
            F_PROFILE(f, t, i) = avg_tc[0];
        else if (y >= tc_heights[NUM_TC-1])
            F_PROFILE(f, t, i) = avg_tc[NUM_TC-1];
        else
            F_PROFILE(f, t, i) = a * y + b;
    }
    end_f_loop(f, t)
}
