#ifndef EXPERIMENT_H
#define EXPERIMENT_H

/*

    Slabs experiment:
    -----------------
    Punctual source of neutrons located in the origin under two slabs with an infinite area, 
    where the first one is filed with a fissionable material (fuel) and the outer
    one if filed with a non-fissionable material (moderator).

*/

#define MODERATOR_SLAB_LIMIT 3.0
#define FUEL_SLAB_LIMIT 1.0
#define CRITICALITY_LIMIT 1.0   /* Limit for the multiplication of particles in the system. */

typedef struct {
    double abs_m;       /* Average number of absorptions.                           */
    double abs_d;       /* Standard deviation of the number of absorptions.         */
    double ref_m;       /* Average number of reflections.                           */
    double ref_d;       /* Standard deviation of the number of reflections.         */
    double esc_m;       /* Average number of escapes.                               */
    double esc_d;       /* Standard deviation of the number of escapes.             */
    double mod_m;       /* Average number of moderations.                           */
    double mod_d;       /* Standard deviation of the number of moderations.         */

    double abs_frac_m;  /* Relative absorbtions average and standard deviation.     */
    double abs_frac_d;
    double ref_frac_m;  /* Relative reflecions average and standard deviation.      */
    double ref_frac_d;
    double esc_frac_m;  /* Relative escapes average and standard deviation.         */
    double esc_frac_d;
    double mod_frac_m;  /* Relative moderations average and standard deviation.     */
    double mod_frac_d;
    
    double hist_m;      /* Average number of particle histories followed.           */
    double hist_d;      /* Standard deviation of the number of histories followed.  */
} result;

/*
    slabs:      Performs r repetitions of the experiment, where 
                in each one a initial number of n histories of neutrons
                will be selected.
*/
void slabs(const unsigned long his, unsigned long rep, result *const results);

#endif

