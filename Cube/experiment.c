#undef __STRICT_ANSI__

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"

void 
cube(const unsigned long his, const unsigned long rep, result *const results) {

    medium cube;
    neutron *particle;
    unsigned long n;
    short additional_n;

    double  moderations,
            absorptions,
            detections,
            reflections,
            escapes,
            distance,
            steps,
            free_path,
            sum_steps,
            sum_fp;

    for(unsigned long r = 0; r < rep; r++) {
        
        n = his;
        particle = (neutron*)malloc(sizeof(neutron) * n);

        moderations = 0;
        absorptions = 0;
        reflections = 0;
        detections  = 0;
        escapes     = 0;
        sum_steps   = 0;
        sum_fp      = 0;

        for(unsigned long h = 0; h < n; h++) {

            distance    = 0;
            steps       = 0;

            // Selects an emission by the source if not all of the initial 
            // particles have been followed yet:
            if(h < his)
                punctual_source(particle, &h);

            particle[h].detection = false;
            
            init_cross_sections(&particle[h], &cube);

            while((FUEL_X1 <= particle[h].x && particle[h].x <= MODERATOR_X2) 
               && (Y1 <= particle[h].y && particle[h].y <= Y2) 
               && (Z1 <= particle[h].z && particle[h].z <= Z2)) {

                distance += select_path(&particle[h], &cube);
                steps++;

                if(particle[h].x < FUEL_X1) {
                    reflections++;
                    break;
                }
                else if((particle[h].x < FUEL_X1 
                        || MODERATOR_X2 < particle[h].x) || (particle[h].y < Y1 
                        || Y2 < particle[h].y) || (particle[h].z < Z1 
                        || Z2 < particle[h].z)) {
                    escapes++;
                    break;
                }
                else if((DETECTOR_X1 <= particle[h].x 
                        && particle[h].x <= DETECTOR_X2) 
                        && (DETECTOR_Y1 <= particle[h].y 
                        && particle[h].y <= DETECTOR_Y2) 
                        && (DETECTOR_Z1 <= particle[h].z 
                        && particle[h].z <= DETECTOR_Z2)) {
                    if(!particle[h].detection) {
                        particle[h].detection = true;
                        detections++;
                    }
                }

                init_cross_sections(&particle[h], &cube);
                additional_n = select_interaction(&particle[h], &cube, &n);

                if(additional_n != 0) {
                    if(additional_n > 0) {
                        particle = realloc(particle, sizeof(neutron) * n);
                        // Store the coordinates of the fission particles to 
                        // follow later:
                        for(unsigned long i = (n - 1); 
                            i >= (n - additional_n); i--) {
                            init_fission_neutrons(&particle[h], &particle[i]);
                        }
                    }
                    absorptions++;
                    break;
                }

                if(particle[h].energy <= THERMAL_ENERGY) {
                    moderations++;
                    break;
                }

            }

            free_path = distance/steps;

            sum_steps += steps;
            sum_fp += free_path;

        }
    
        results->mod_m += moderations;
        results->abs_m += absorptions;
        results->ref_m += reflections;
        results->det_m += detections;
        results->esc_m += escapes;
        results->mod_frac_m += moderations/n;
        results->abs_frac_m += absorptions/n;
        results->ref_frac_m += reflections/n;
        results->det_frac_m += detections/n;
        results->esc_frac_m += escapes/n;
        results->hist_m += (double)n;
        results->stp_m += sum_steps/n;
        results->frp_m += sum_fp/n;

        if(rep > 1) {
            results->mod_d += SQUARE(moderations);
            results->abs_d += SQUARE(absorptions);
            results->ref_d += SQUARE(reflections);
            results->det_d += SQUARE(detections);
            results->esc_d += SQUARE(escapes);
            results->mod_frac_d += SQUARE(moderations/n);
            results->abs_frac_d += SQUARE(absorptions/n);
            results->ref_frac_d += SQUARE(reflections/n);
            results->det_frac_d += SQUARE(detections/n);
            results->esc_frac_d += SQUARE(escapes/n);
            results->hist_d += SQUARE((double)n);
            results->stp_d += SQUARE(sum_steps/n);
            results->frp_d += SQUARE(sum_fp/n);
        }

        free(particle);

    }
    
    // Mean of the number and the relative number of occurrences of each type 
    // of interaction and of the number of histories followed in each repetition 
    // of the experiment:
    results->mod_m /= rep;
    results->abs_m /= rep;
    results->ref_m /= rep;
    results->det_m /= rep;
    results->esc_m /= rep;
    results->mod_frac_m /= rep;
    results->abs_frac_m /= rep;
    results->ref_frac_m /= rep;
    results->det_frac_m /= rep;
    results->esc_frac_m /= rep;
    results->hist_m /= rep;
    results->stp_m /= rep;
    results->frp_m /= rep;

    if(rep > 1) {
        // Standard deviation:
        results->mod_d = sqrt((results->mod_d/rep) 
                            - SQUARE(results->mod_m))/(rep - 1);
        results->abs_d = sqrt((results->abs_d/rep) 
                            - SQUARE(results->abs_m))/(rep - 1);
        results->ref_d = sqrt((results->ref_d/rep) 
                            - SQUARE(results->ref_m))/(rep - 1);
        results->det_d = sqrt((results->det_d/rep) 
                            - SQUARE(results->det_m))/(rep - 1);
        results->esc_d = sqrt((results->esc_d/rep) 
                            - SQUARE(results->esc_m))/(rep - 1);
        results->mod_frac_d = sqrt((results->mod_frac_d/rep) 
                            - SQUARE(results->mod_frac_m))/(rep - 1);
        results->abs_frac_d = sqrt((results->abs_frac_d/rep) 
                                - SQUARE(results->abs_frac_m))/(rep - 1);
        results->ref_frac_d = sqrt((results->ref_frac_d/rep) 
                                - SQUARE(results->ref_frac_m))/(rep - 1);
        results->det_frac_d = sqrt((results->det_frac_d/rep) 
                                - SQUARE(results->det_frac_m))/(rep - 1);
        results->esc_frac_d = sqrt((results->esc_frac_d/rep) 
                                - SQUARE(results->esc_frac_m))/(rep - 1);
        results->hist_d = sqrt((results->hist_d/rep) 
                                - SQUARE(results->hist_m))/(rep - 1);
        results->stp_d = sqrt((results->stp_d/rep) 
                                - SQUARE(results->stp_m))/(rep - 1);
        results->frp_d = sqrt((results->frp_d/rep) 
                                - SQUARE(results->frp_m))/(rep - 1);
    }

}