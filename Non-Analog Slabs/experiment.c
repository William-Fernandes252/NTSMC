#undef __STRICT_ANSI__

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"

void slabs(const unsigned long his, const unsigned long rep, result *const results) {

    medium slabs;
    neutron *particle;
    unsigned long n;

    double  moderations,
            reflections,
            escapes; 

    for(unsigned long r = 0; r < rep; r++) {
        
        n = his;
        particle = (neutron*)malloc(sizeof(neutron) * n);

        moderations = 0;
        reflections = 0;
        escapes     = 0;

        for(unsigned long h = 0; h < n; h++) {

            // Selects an emission by the source if not all of the initial particles have been followed yet:
            if(h < his)
                punctual_source(particle, &h);
            
            init_cross_sections(&particle[h], &slabs);

            while(0.0 <= particle[h].z && particle[h].z <= MODERATOR_SLAB_LIMIT) {

                select_path(&particle[h], &slabs);

                if(particle[h].z < 0.0) {
                    reflections += particle[h].weight;
                    break;
                }
                else if(MODERATOR_SLAB_LIMIT < particle[h].z) {
                    escapes += particle[h].weight;
                    break;
                }

                if(select_interaction(&particle[h], &slabs)) {
                    n += 1;
                    particle = (neutron*)realloc(particle, sizeof(neutron) * n);
                    // Store the coordinates of the fission particle to follow later:
                    init_fission_neutrons(&particle[h], &particle[n-1], &slabs);
                }

                if(particle[h].weight < MIN_WEIGHT) {
                    if(russian_roulette(&particle[h]))
                        break;
                }

                if(particle[h].energy <= THERMAL_ENERGY) {
                    moderations += particle[h].weight;
                    break;
                }

            }

        } 
        
        results->mod_m += moderations;
        results->ref_m += reflections;
        results->esc_m += escapes;
        results->mod_frac_m += moderations/n;
        results->ref_frac_m += reflections/n;
        results->esc_frac_m += escapes/n;
        results->hist_m += (double)n;

        if(rep > 1) {
            results->mod_d += SQUARE(moderations);
            results->ref_d += SQUARE(reflections);
            results->esc_d += SQUARE(escapes);
            results->mod_frac_d += SQUARE(moderations/n);
            results->ref_frac_d += SQUARE(reflections/n);
            results->esc_frac_d += SQUARE(escapes/n);
            results->hist_d += (double)n;
            results->hist_d += SQUARE((double)n);
        }

        free(particle);

    }

    // Mean of the number and the relative number of occurrences of each type of interaction and of the number of histories followed in each repetition of the experiment:
    results->mod_m /= rep;
    results->ref_m /= rep;
    results->esc_m /= rep;
    results->mod_frac_m /= rep;
    results->ref_frac_m /= rep;
    results->esc_frac_m /= rep;
    results->hist_m /= rep;

    if(rep > 1) {
        // Standard deviation:
        results->mod_d = sqrt((results->mod_d/rep) - SQUARE(results->mod_m))/(rep - 1);
        results->ref_d = sqrt((results->ref_d/rep) - SQUARE(results->ref_m))/(rep - 1);
        results->esc_d = sqrt((results->esc_d/rep) - SQUARE(results->esc_m))/(rep - 1);
        results->mod_frac_d = sqrt((results->mod_frac_d/rep) - SQUARE(results->mod_frac_m))/(rep - 1);
        results->ref_frac_d = sqrt((results->ref_frac_d/rep) - SQUARE(results->ref_frac_m))/(rep - 1);
        results->esc_frac_d = sqrt((results->esc_frac_d/rep) - SQUARE(results->esc_frac_m))/(rep - 1);
        results->hist_d = sqrt((results->hist_d/rep) - SQUARE(results->hist_m))/(rep - 1);
    }

}