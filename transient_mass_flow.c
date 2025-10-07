#include "udf.h"
#include <stdio.h>

#define MAX_DATA_POINTS 1000

static real time_data[MAX_DATA_POINTS];
static real massflow_data[MAX_DATA_POINTS];
static int num_points = 0;
static int data_loaded = 0;

void load_massflow_csv()
{
    FILE *fp = fopen("mass_flow_rate.csv", "r");
    if (!fp) return;

    char line[256];
    fgets(line, sizeof(line), fp); // skip header

    while (fgets(line, sizeof(line), fp) && num_points < MAX_DATA_POINTS)
    {
        sscanf(line, "%lf,%lf", &time_data[num_points], &massflow_data[num_points]);
        num_points++;
    }

    fclose(fp);
    data_loaded = 1;
}

real interpolate_massflow(real time)
{
    if (num_points < 2) return 0.0;

    if (time <= time_data[0])
        return massflow_data[0];

    for (int i = 0; i < num_points - 1; i++)
    {
        if (time >= time_data[i] && time <= time_data[i+1])
        {
            real t1 = time_data[i], t2 = time_data[i+1];
            real m1 = massflow_data[i], m2 = massflow_data[i+1];
            return m1 + (m2 - m1) * (time - t1) / (t2 - t1);
        }
    }

    return massflow_data[num_points - 1]; // fallback for time > last data point
}

DEFINE_PROFILE(transient_massflow_profile, t, i)
{
    face_t f;
    real time = CURRENT_TIME;

    if (!data_loaded)
        load_massflow_csv();

    real m_dot = interpolate_massflow(time);

    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = m_dot;
    }
    end_f_loop(f, t)
}
