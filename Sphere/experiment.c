
// To use the math constants from math.h in a code compiled with the MinGW gcc.
#undef __STRICT_ANSI__  

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"

double 
volume(unsigned int layer) {

    double r = (double)layer * (RADIUS/LAYERS);

    return (4.0 * M_PI * (pow((r + (RADIUS/LAYERS)), 3.0) - pow(r, 3.0)))/3.0;

}

void 
sphere(const unsigned long his, const unsigned long rep, 
       result *const results) {

    medium sphere;
    neutron *particle;
    unsigned long n;
    short additional_n;

    double  n_collisions[LAYERS],   /* Stores the number of collisions that 
                                       occurred in each layer of the sphere. */
            n_escapes[LAYERS],      /* Stores the number of particles that 
                                       left the system from each region.     */
            d;                      /* Keep track of the distance between 
                                       a particle and the source.            */
    unsigned int curr_layer;        /* Keep track of the layer of the sphere 
                                       where the particle is located.        */
    unsigned int last_layer;        /* Keep track of the layer of the sphere 
                                       where the particle was located before 
                                       its next random random walk.          */
            
    for(unsigned long r = 0; r < rep; r++) {
        
        n = his;
        particle = (neutron*)malloc(sizeof(neutron) * n);
        for(unsigned int l = 0; l < LAYERS; l++) {
            n_collisions[l] = 0.00000;
            n_escapes[l] = 0.00000;
        }

        for(unsigned long h = 0; h < n; h++) {

            additional_n = 0;

            // Selects an emission by the source if not all of the initial 
            // particles have been followed yet:
            if(h < his)
                punctual_source(particle, &h);

            init_cross_sections(&particle[h], &sphere, 0.0);
            
            while(true) {
                
                last_layer = curr_layer;

                select_path(&particle[h], &sphere);

                d = DISTANCE(particle[h].x, particle[h].y, particle[h].z);
                curr_layer = (unsigned int)d * (LAYERS/RADIUS);  

                if(d > RADIUS) {
                    n_escapes[last_layer]++;
                    break;
                }

                init_cross_sections(&particle[h], &sphere, d);
                additional_n = select_interaction(&particle[h], &sphere, &n);

                n_collisions[curr_layer]++;
                
                if(additional_n != 0) {
                    if(additional_n > 0) {
                        particle = realloc(particle, sizeof(neutron) * n);
                        // Store the coordinates of the fission particles 
                        // to follow later:
                        for(unsigned long i = (n - 1); 
                            i >= (n - additional_n); i--)
                            init_fission_neutrons(&particle[h], &particle[i]);
                    }
                    break;
                }

                if(particle[h].energy <= THERMAL_ENERGY)
                    break;
                
            }
             
        }
        
        double total_flux, layer_r;

        if(rep > 1) {
            for(unsigned int l = 0; l < LAYERS; l++) {

                results->coll_mean[l] += n_collisions[l];
                results->esc_mean[l] += n_escapes[l];
                results->coll_dev[l] += SQUARE(n_collisions[l]);
                results->esc_dev[l] += SQUARE(n_escapes[l]);

                layer_r = l * (RADIUS/LAYERS);
                if(layer_r < (RADIUS/FUEL_PART))
                    total_flux = n_collisions[l]/volume(l)/n/FUEL_TCS;
                else if(layer_r <= RADIUS)
                    total_flux = n_collisions[l]/volume(l)/n/MOD_TCS;

                results->flux_mean[l] += total_flux;
                results->flux_dev[l] += SQUARE(total_flux);
                
            }
        }
        else {
            for(unsigned int l = 0; l < LAYERS; l++) {

                results->coll_mean[l] += n_collisions[l];
                results->esc_mean[l] += n_escapes[l];

                layer_r = l * (RADIUS/LAYERS);
                if(layer_r < (RADIUS/FUEL_PART))
                    total_flux = n_collisions[l]/volume(l)/n/FUEL_TCS;
                else if(layer_r <= RADIUS)
                    total_flux = n_collisions[l]/volume(l)/n/MOD_TCS;

                results->flux_mean[l] += total_flux;
                
            }
        }

        free(particle);

    }

    for(unsigned int l = 0; l < LAYERS; l++) {
        results->coll_mean[l] /= rep;
        results->esc_mean[l] /= rep;
        results->flux_mean[l] /= rep;
    }

    if(rep > 1) {
        for(unsigned int l = 0; l < LAYERS; l++) {
            results->coll_dev[l] = sqrt((results->coll_dev[l]/rep) 
                                 - SQUARE(results->coll_mean[l]))/(rep - 1);
            results->esc_dev[l] = sqrt((results->esc_dev[l]/rep) 
                                - SQUARE(results->esc_mean[l]))/(rep - 1);
            results->flux_dev[l] = sqrt((results->flux_dev[l]/rep) 
                                 - SQUARE(results->flux_mean[l]))/(rep - 1);
        }
    }

}