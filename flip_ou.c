#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

#define T_SIMULATION 100
#define MIN(a,b) ( (a < b) ? (a) : (b))

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

	double mu = 1.0	/* flip rate from trapping to repellent*/, nu = 1.0 /* flip rate from R to T*/;
	double beta_abs = 1.0; /* Absolute value of strength of harmonic potential*/
	double beta;
	double diffusivity_0 = 1.0;
	double diffusivity;
	double sum_rates = mu + nu;
	double probability_trapping = nu/sum_rates;
	double probability_repellent = mu/sum_rates;

	double x;

	double time;
	double delta_t = 0.001;

	double last_flip_time;
	double next_flip_time;
	double write_out_time = 0.1;
	int write_out_index;

	int write_points = ((int) (T_SIMULATION/write_out_time));
	double *x_variance_average = calloc(write_points, sizeof(double));

	int iterations = 100000;
	double inv_iteration = (1.0/((double) iterations));

	int iter = 1;

	do
	{
		x = 0.0;
		time = 0.0;
		last_flip_time = 0.0;
		write_out_index = 1;
		int current_state = gsl_ran_binomial(r, probability_trapping, 1); // Draw initial state of trap
		do{
			if( current_state == 1)
			{
				// Trapping
				beta = beta_abs;
				next_flip_time = last_flip_time + gsl_ran_exponential(r, 1/mu);
			}
			else{
				// Repellent
				beta = -beta_abs;
				next_flip_time = last_flip_time + gsl_ran_exponential(r, 1/nu);
			}
			
			//delta_t = ((next_flip_time - last_flip_time)/1000.);
			diffusivity = sqrt(2 * diffusivity_0 * delta_t);
			for(time = last_flip_time; time <= MIN( next_flip_time, T_SIMULATION); time += delta_t)
			{

				x += (- beta*x*delta_t + gsl_ran_gaussian_ziggurat(r,diffusivity)); //Check 2's)
				if(time >= write_out_index*write_out_time)
				{ 
					x_variance_average[write_out_index] = ((((double)(iter-1))/((double)iter))*x_variance_average[write_out_index] + x*x/((double) iter));
					write_out_index++;
				}
			}

			current_state = (1 - current_state); //Flip from T <-> R
			last_flip_time = next_flip_time;
		}
		while(next_flip_time < T_SIMULATION);
	iter++;
	}while(iter <= iterations);

	//Output
	int i;
	for(i = 1; i <= write_points; i++)
	{
		printf("%g\t%g\n", i*write_out_time, x_variance_average[i]);
	}


	return 0;
}
