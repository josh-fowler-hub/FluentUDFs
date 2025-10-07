#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ROWS 1000
#define NUM_TC 6

static const real max_tc_height = 45.0;
static const real tc_heights[NUM_TC] = {3.797244094/max_tc_height, 11.39173228/max_tc_height, 18.98622047/max_tc_height, 26.58070866/max_tc_height, 34.17519685/max_tc_height, 41.76968504/max_tc_height};
static real time_data[MAX_ROWS];
static real tc_data[MAX_ROWS][NUM_TC];
static int num_rows = 0;
static int data_loaded = 0;

void load_csv_data()
{
    FILE *fp = fopen("thermocouple_data.csv", "r");
    if (!fp) return;

    char line[256];
    fgets(line, sizeof(line), fp); // skip header

    while (fgets(line, sizeof(line), fp) && num_rows < MAX_ROWS)
    {
        double t, vals[NUM_TC];
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &t, &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5]) == 7)
        {
            time_data[num_rows] = t;
            for (int i = 0; i < NUM_TC; i++)
                tc_data[num_rows][i] = vals[i];
            num_rows++;
        }
    }

    data_loaded = 1;
}

DEFINE_PROFILE(dynamic_tinf_profile, t, i)
{
    face_t f;
    real a = 0.0, b = 0.0;
    real time = CURRENT_TIME;

    if (!data_loaded)
        load_csv_data();

    // Find bracketing time steps
    int i1 = 0;
    while (i1 < num_rows - 1 && time_data[i1 + 1] < time)
        i1++;

    int i2 = (i1 < num_rows - 1) ? i1 + 1 : i1;

    real t1 = time_data[i1], t2 = time_data[i2];
    real interp_tc[NUM_TC];

    for (int j = 0; j < NUM_TC; j++)
    {
        real v1 = tc_data[i1][j];
        real v2 = tc_data[i2][j];
        interp_tc[j] = (t2 != t1) ? v1 + (v2 - v1) * (time - t1) / (t2 - t1) : v1;
    }

    // Linear fit: T = a*y + b
    real sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (int j = 0; j < NUM_TC; j++)
    {
        real x = tc_heights[j];
        real y = interp_tc[j];
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    real denom = NUM_TC * sum_x2 - sum_x * sum_x;
    if (denom != 0)
    {
        a = (NUM_TC * sum_xy - sum_x * sum_y) / denom;
        b = (sum_y * sum_x2 - sum_x * sum_xy) / denom;
    }

    begin_f_loop(f, t)
    {
        real y = F_Y(f, t);
        if (y <= tc_heights[0])
            F_PROFILE(f, t, i) = interp_tc[0];
        else if (y >= tc_heights[NUM_TC-1])
            F_PROFILE(f, t, i) = interp_tc[NUM_TC-1];
        else
            F_PROFILE(f, t, i) = a * y + b;
    }
    end_f_loop(f, t)
}
