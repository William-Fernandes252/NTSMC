#ifndef EXPERIMENT_H
#define EXPERIMENT_H

/* 

    Sphere experiment: 
    ------------------
    Punctual source of neutrons located in the center of a sphere divided 
    into determined number of layers, where the inner ones are filed with 
    a fissionable material (fuel), and the outer ones are filed with the 
    a non-fissionable material (moderator).

*/

#define RADIUS 1.00000e2
#define LAYERS 100

/* Relative section of the sphere filed with the fissionable material. */
#define FUEL_PART 100

typedef struct {
    /* Average number of collisions in each cell of the sphere.     */
    double coll_mean[LAYERS];
    /* Standard deviation of the number of collisions in each cell. */        
    double coll_dev[LAYERS];
    /* Average number of escapes in each layer of the sphere.       */        
    double esc_mean[LAYERS];
    /* Standard deviation of the number of escapes in each layer.   */       
    double esc_dev[LAYERS];
    /* Average Flux per cell of the sphere.                         */       
    double flux_mean[LAYERS];
    /* Standard deviation of the flux per cell.                     */
    double flux_dev[LAYERS];
} result;

/*
    sphere:     Performs r repetitions of the experiment, where 
                in each one a initial number of n histories of neutrons
                will be selected.
*/
void 
sphere(const unsigned long his, const unsigned long rep, 
       result *const results);

/*
    volume:     Computes the volume of the i-th layer of the sphere.
*/
double 
volume(unsigned int i);

#endif

