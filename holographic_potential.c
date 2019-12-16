#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

#define T_SIMULATION 10000
#define X_RESOLUTION 200
#define MIN(a,b) ( (a < b) ? (a) : (b))

void bin_x(double, double*, double, double);

// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;
int seed;


int main(){

	//Init RNG
	seed = -1;
 	gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        if(seed==-1) seed = ((int) (((int) clock() ) % 100000));
        gsl_rng_set(r, seed);

	double flip_threshold = 0.25; // This is epsilon, the distance to the origin away of which the potentia gets flippd
	double beta = 1.0; /* Absolute value of strength of harmonic potential*/
	double offset = 5.0; // = a_0, the two potentials are +/- offset away from the origin
	double diffusivity_0 = 1.0;
	
	double delta_t = 0.001;
	double diffusivity = sqrt(2 * diffusivity_0 * delta_t);

	double x, x_min, x_max;

	double time;

	double write_out_time = 0.1;
	int write_out_index;

	int write_points = ((int) (T_SIMULATION/write_out_time));
	double x_variance_average;
	double x_variance;

	double *x_histogram = calloc(X_RESOLUTION, sizeof(double));

	int iterations = 1;
	double inv_iteration = (1.0/((double) iterations));


	int current_state; // 1-potential pulls to the right, 0 - pulls to the left
	int switch_counter;
	int iter = 1;
	long time_steps_counter = 0;

	do
	{
		offset = 5.0;
		// One realisation of experiment
		x = 0.0;
		time = 0.0;
		write_out_index = 1;
		current_state = 1;
		switch_counter = 0;
		x_variance= 0.0;
		time_steps_counter = 0;
		x_min = -2.0;
		x_max = 2.0;
		
		int i;

		do{
			if( current_state == 1)
			{
				// Potential pulls to right
				do
				{
					x += (- beta*(x - offset) *delta_t + gsl_ran_gaussian_ziggurat(r,diffusivity)); //Check 2's)
					x_variance += x*x;
					time_steps_counter++;
					bin_x(x, x_histogram, x_min, x_max);
					//printf("%g\t%g\t%i\n", time_steps_counter*delta_t, x, current_state);
				
				}
				while((x < flip_threshold) && (time_steps_counter*delta_t < T_SIMULATION));
			}
			else{
				// Potential pulls to the left
				do
				{
					x += (- beta*(x + offset) *delta_t + gsl_ran_gaussian_ziggurat(r,diffusivity)); //Check 2's)
				
					x_variance += x*x;
					time_steps_counter++;
					bin_x(x, x_histogram, x_min, x_max);
					//printf("%g\t%g\t%i\n", time_steps_counter*delta_t, x, current_state);
				}
				while((x >  -flip_threshold) && (time_steps_counter*delta_t < T_SIMULATION));
			}

			current_state = 1 - current_state;
			switch_counter++;
			
		}
		while(time_steps_counter*delta_t < T_SIMULATION);
		
		double delta_x = ((x_max-x_min)/(X_RESOLUTION));
		for(i = 0; i < X_RESOLUTION; i++)
		{
			printf("%g\t%g\n", (x_min + (((double)i)+0.5)*delta_x), x_histogram[i]/((double) time_steps_counter) );
		}
		//printf("%g\t%g\t%g\n", offset, x_variance/((double) time_steps_counter), (((double) switch_counter)/T_SIMULATION));
	iter++;
	}while(iter <= iterations);
	
	return 0;
}

void bin_x(double x, double* histogram, double x_min, double x_max)
{
	if((x_min < x) && (x < x_max))
	{
		double delta_x = ((x_max-x_min)/(X_RESOLUTION));
		int i = ceil((x-x_min)/delta_x);
		histogram[i] += 1.0;
	}
}
