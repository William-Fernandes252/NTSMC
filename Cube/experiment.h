#ifndef EXPERIMENT_H
#define EXPERIMENT_H

/*

    Cube experiment:
    -----------------
    Punctual source of neutrons located in the center of a slab filled with
    a fissionable material (fuel), which is unitted with a cube filled with
    a non-fissionable material (moderator). There is no barriers between the
    regions, so the particles can travel from one to another with no 
    limitations.
    In the center of the cube there is a smaller cube, which is the region
    monitored by a detector.

*/

#define MODERATOR_X1 0.00
#define MODERATOR_X2 10.0
#define FUEL_X1 -1.00
#define FUEL_X2 0.00
#define Y1 0.00
#define Y2 10.0
#define Z1 0.00
#define Z2 10.0
#define DETECTOR_EDGE 7.00
#define DETECTOR_X1 (((MODERATOR_X2 + MODERATOR_X1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_X2 (((MODERATOR_X2 + MODERATOR_X1)/2) + (DETECTOR_EDGE/2))
#define DETECTOR_Y1 (((Y2 + Y1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_Y2 (((Y2 + Y1)/2) + (DETECTOR_EDGE/2))
#define DETECTOR_Z1 (((Z2 + Z1)/2) - (DETECTOR_EDGE/2))
#define DETECTOR_Z2 (((Z2 + Z1)/2) + (DETECTOR_EDGE/2))


typedef struct {
    double mod_m;       /* Average number of moderations.                   */
    double mod_d;       /* Standard deviation of the number of moderations. */
    double abs_m;       /* Average number of absorptions.                   */
    double abs_d;       /* Standard deviation of the number of absorptions. */
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

    double mod_frac_m;  /* Relative moderations average and std. deviation. */
    double mod_frac_d;
    double abs_frac_m;  /* Relative absorbtions average and std. deviation. */
    double abs_frac_d;
    double ref_frac_m;  /* Relative reflecions average and std. deviation.  */
    double ref_frac_d;
    double det_frac_m;  /* Relative detections average and std. deviation.  */
    double det_frac_d;
    double esc_frac_m;  /* Relative escapes average and std. deviation.     */
    double esc_frac_d;
    
    double hist_m;      /* Average number of particle histories followed.     */
    double hist_d;      /* Std. deviation of the number of histories followed.*/
} result;

/*
    cube:   Performs r repetitions of the experiment, where 
            in each one a initial number of n histories of neutrons
            will be selected.
*/
void 
cube(const unsigned long his, unsigned long rep, result *const results);

#endif

