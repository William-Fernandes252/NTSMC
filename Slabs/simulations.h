#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#define SQUARE(x) ((x)*(x))
#define ABS(x) ((x) >= 0.0 ? (x) : (-(x)))

// Particles mass number and thermal energy:
#define N_MASS_NUM 1.00000e0
#define THERMAL_ENERGY 5.00000e-1

// Fuel (U-235) mass number and cross-sections:
#define FUEL_MASS_NUM 2.35000e2
#define FUEL_TCS 2.00000e2
#define FUEL_FCS 5.00000e0
#define FUEL_ESCS 1.25000e2
#define FUEL_ISCS 5.00000e1
#define FUEL_ACS 2.00000e1
#define FISSION_NEUTRONS 2.43000e0 /* Average number of neutrons emitted after a fission event. */

// Moderator (H20-18) mass number and cross-sections:
#define MOD_MASS_NUM 1.80000e1
#define MOD_TCS 2.5000e2
#define MOD_ESCS 1.50000e2
#define MOD_ISCS 9.50000e1
#define MOD_ACS 5.00000e0

// Medium properties:
typedef struct {
    double mass_num;        /* Mass number of atoms that compounds the fuel.            */
    double t_cs;            /* Total cross-section of the fuel.                         */
    double el_cs;           /* Cross-section for elastic scattering.                    */
    double inl_cs;          /* Cross-section for inelastic scattering.                  */
    double fission_cs;      /* Fission cross-section of the fuel.                       */
    double abs_cs;          /* Absorption cross-section of the fuel.                    */
    double n_fission;       /* Average number of neutrons emitted in fission events.    */
} fuel;

typedef struct {
    double mass_num;        /* Mass number of atoms that compounds the moderator.       */
    double t_cs;            /* Total cross-section of the moderator.                    */
    double el_cs;           /* Cross-section for elastic scattering.                    */
    double inl_cs;          /* Cross-section for inelastic scattering.                  */
    double abs_cs;          /* Absorption cross-section of the moderator.               */
} moderator;

typedef enum {
    FUEL,
    MODERATOR
} region;     /* Tag to identify regions. */

typedef union {
    fuel u_diox;       /* Fuel simulated in the experiment (uranium dioxide).   */
    moderator water;   /* Moderator simulated in the experiment.                */
} material;

typedef struct {
    region tag;
    material comp; 
} medium;

// Particles states:
typedef struct {
    double x, y, z; /* Cartesian coordinates of the particle position.                              */
    double u, v, w; /* Direction vector of the particle where.                                      */
    double in, out; /* Variables to keep check if the particle is getting out of a region bound.    */
    double energy;  /* Energy of the particle.                                                      */
} neutron;

// Functions for the simulation of particles histories:

/*
    ini_cross_sections:     Given the particle coordinates in the phase space, compute 
                            the cross sections of the material in the surroundings.
*/
void init_cross_sections(const neutron *const particle, medium *const slabs);

/*
    init_fission_neutrons:  Given the coordinates of the incident particle, compute the 
                            coordites of the particles emitted by the fission event and
                            stores then in the memory to be followed later.
*/
void init_fission_neutrons(const neutron *const incident_n, neutron *const fission_n);

/*
    check_path:     Given the location of the particle in the medium and its
                    current direction, checks if the particle is entering or 
                    leaving a region.
*/
bool check_path(neutron *const particle, const medium *const slabs);

/*
    select_path:     Given the location of the particle in the medium, selects 
                            the walk of the particle to the site where its next
                            interaction with the medium will happen.
*/
void select_path(neutron *const particle, const medium *const slabs);

/*
    select_interaction:     Given the location of the particle in the medium, selects a
                            a interaction with some nuclide in the surroundings. Returns the 
                            number of additional particles in the medium after the event.
*/
short select_interaction(neutron *particle, const medium *const slabs, unsigned long *const n);

/*
    direction_after_scattering:     Selects the new direction of a neutron 
                                    after a collision event using a rejection
                                    method.
                                    
*/
void direction_after_colission(neutron *const particle, double cos_theta);

/*
    iso_elastic_scattering:     Simulates the scattering of particles, that takes 
                                place after collisions events, with the assumption that the 
                                distribution of the new direction after a non-absorption 
                                interaction is isotropic. For simplicity, the cosine of the 
                                angle of the collision is selected from de elastic scattering 
                                model by default. If a inelastic collision was selected, the
                                energy of the particle is reduced by a factor depending on the
                                target element.  
*/
void iso_scattering(neutron *const particle, const double target_mass, bool inelastic);

/*
    iso_direction:  Selects a direction for a particle from a isotropic distribution. 
*/
void iso_direction(neutron *const particle);

/*
    punctual_source:    Simulates the emission of particles by a punctual font 
                        located in a selected point of the sistem. Returns the
                        adress of the coordinates of the particles to be followed. 
*/
void punctual_source(neutron *particle, const unsigned long *const h);

#endif