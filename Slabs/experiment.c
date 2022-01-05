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
    short additional_n;

    double  moderations,
            absorptions,
            reflections,
            escapes; 

    for(unsigned long r = 0; r < rep; r++) {
        
        n = his;
        particle = (neutron*)malloc(sizeof(neutron) * n);

        moderations = 0;
        absorptions = 0;
        reflections = 0;
        escapes     = 0;

        for(unsigned long h = 0; h < n; h++) {
            
            additional_n = 0;

            // Selects an emission by the source if not all of the initial particles have been followed yet:
            if(h < his)
                punctual_source(particle, &h);
            
            init_cross_sections(&particle[h], &slabs);

            while(0.0 <= particle[h].z && particle[h].z <= MODERATOR_SLAB_LIMIT) {

                select_path(&particle[h], &slabs);

                if(particle[h].z < 0.0) {
                    reflections++;
                    break;
                }
                else if(MODERATOR_SLAB_LIMIT < particle[h].z) {
                    escapes++;
                    break;
                }

                init_cross_sections(&particle[h], &slabs);
                additional_n = select_interaction(&particle[h], &slabs, &n);

                if(additional_n != 0) {
                    if(additional_n > 0) {
                        particle = (neutron*)realloc(particle, sizeof(neutron) * n);
                        // Store the coordinates of the fission particles to follow later:
                        for(unsigned long i = (n - 1); i >= (n - additional_n); i--)
                            init_fission_neutrons(&particle[h], &particle[i]);
                    }
                    absorptions++;
                    break;
                }

                if(particle[h].energy <= THERMAL_ENERGY) {
                    moderations++;
                    break;
                }

            }

        } 
    
        results->mod_m += moderations;
        results->abs_m += absorptions;
        results->ref_m += reflections;
        results->esc_m += escapes;
        results->mod_frac_m += moderations/n;
        results->abs_frac_m += absorptions/n;
        results->ref_frac_m += reflections/n;
        results->esc_frac_m += escapes/n;
        results->hist_m += (double)n;

        if(rep > 1) {
            results->mod_d += SQUARE(moderations);
            results->abs_d += SQUARE(absorptions);
            results->ref_d += SQUARE(reflections);
            results->esc_d += SQUARE(escapes);
            results->mod_frac_d += SQUARE(moderations/n);
            results->abs_frac_d += SQUARE(absorptions/n);
            results->ref_frac_d += SQUARE(reflections/n);
            results->esc_frac_d += SQUARE(escapes/n);
            results->hist_d += (double)n;
            results->hist_d += SQUARE((double)n);
        }

        free(particle);

    }
    
    // Mean of the number and the relative number of occurrences of each type of interaction and of the number of histories followed in each repetition of the experiment:
    results->mod_m /= rep;
    results->abs_m /= rep;
    results->ref_m /= rep;
    results->esc_m /= rep;
    results->mod_frac_m /= rep;
    results->abs_frac_m /= rep;
    results->ref_frac_m /= rep;
    results->esc_frac_m /= rep;
    results->hist_m /= rep;

    if(rep > 1) {
        // Standard deviation:
        results->mod_d = sqrt((results->mod_d/rep) - SQUARE(results->mod_m))/(rep - 1);
        results->abs_d = sqrt((results->abs_d/rep) - SQUARE(results->abs_m))/(rep - 1);
        results->ref_d = sqrt((results->ref_d/rep) - SQUARE(results->ref_m))/(rep - 1);
        results->esc_d = sqrt((results->esc_d/rep) - SQUARE(results->esc_m))/(rep - 1);
        results->mod_frac_d = sqrt((results->mod_frac_d/rep) - SQUARE(results->mod_frac_m))/(rep - 1);
        results->abs_frac_d = sqrt((results->abs_frac_d/rep) - SQUARE(results->abs_frac_m))/(rep - 1);
        results->ref_frac_d = sqrt((results->ref_frac_d/rep) - SQUARE(results->ref_frac_m))/(rep - 1);
        results->esc_frac_d = sqrt((results->esc_frac_d/rep) - SQUARE(results->esc_frac_m))/(rep - 1);
        results->hist_d = sqrt((results->hist_d/rep) - SQUARE(results->hist_m))/(rep - 1);
    }

}