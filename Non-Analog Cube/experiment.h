#ifndef EXPERIMENT_H
#define EXPERIMENT_H

/*

    Cube experiment:
    -----------------
    Punctual source of neutrons located in the center of a slab filled with
    a fissionable material (fuel), which is unitted with a cube filled with
    a non-fissionable material (moderator). There is no barriers between the
    regions, so the particles can travel from one to another with no limitations.
    in the center of the cube there is a smaller cube, which is the region
    monitored by a detector. Statistical weight avaliation and Russian Roulette
    are played, so the simulations are non-analog.

*/

#define MODERATOR_X1 10.0
#define MODERATOR_X2 40.0
#define FUEL_X1 0.0
#define FUEL_X2 10.0
#define Y1 0.0
#define Y2 30.0
#define Z1 0.0
#define Z2 30.0
#define DETECTOR_EDGE 30.0
#define DETECTOR_X1 (((MODERATOR_X2 + MODERATOR_X1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_X2 (((MODERATOR_X2 + MODERATOR_X1)/2) + (DETECTOR_EDGE/2))
#define DETECTOR_Y1 (((Y2 + Y1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_Y2 (((Y2 + Y1)/2) + (DETECTOR_EDGE/2))
#define DETECTOR_Z1 (((Z2 + Z1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_Z2 (((Z2 + Z1)/2) + (DETECTOR_EDGE/2))


typedef struct {
    double mod_m;       /* Average number of moderations.                   */
    double mod_d;       /* Standard deviation of the number of moderations. */
    double ref_m;       /* Average number of reflections.                   */
    double ref_d;       /* Standard deviation of the number of reflections. */
    double det_m;       /* Average number of detections.                    */
    double det_d;       /* Standard deviation of the number of detections.  */
    double esc_m;       /* Average number of escapes.                       */
    double esc_d;       /* Standard deviation of the number of escapes.     */
    double stp_m;       /* Average number of steps.                         */
    double stp_d;       /* Standard deviation of the number of steps.       */
    double frp_m;       /* Mean free path of the particles in the system.   */
    double frp_d;       /* Standard deviation of the mean free path.        */

    double mod_frac_m;  /* Relative moderations average and standard deviation.     */
    double mod_frac_d;
    double ref_frac_m;  /* Relative reflecions average and standard deviation.      */
    double ref_frac_d;
    double det_frac_m;  /* Relative detections average and standard deviation.      */
    double det_frac_d;
    double esc_frac_m;  /* Relative escapes average and standard deviation.         */
    double esc_frac_d;
    
    double hist_m;      /* Average number of particle histories followed.           */
    double hist_d;      /* Standard deviation of the number of histories followed.  */
} result;

/*
    cube:      Performs r repetitions of the experiment, where 
                in each one a initial number of n histories of neutrons
                will be selected.
*/
void cube(const unsigned long his, unsigned long rep, result *const results);

#endif

