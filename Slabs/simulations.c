#undef __STRICT_ANSI__

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"

void init_cross_sections(const neutron *const particle, medium *const slabs) {

    // Compute the cross-sections from the particle atual position in the medium:
    if(particle->z <= FUEL_SLAB_LIMIT) {
        slabs->tag = FUEL;
        slabs->comp.u_diox.mass_num = FUEL_MASS_NUM;
        slabs->comp.u_diox.t_cs = FUEL_TCS;
        slabs->comp.u_diox.el_cs = FUEL_ESCS;
        slabs->comp.u_diox.inl_cs = FUEL_ISCS;
        slabs->comp.u_diox.fission_cs = FUEL_FCS;
        slabs->comp.u_diox.abs_cs = FUEL_ACS;
        slabs->comp.u_diox.n_fission = FISSION_NEUTRONS;
    }
    else if(particle->z <= MODERATOR_SLAB_LIMIT) {
        slabs->tag = MODERATOR;
        slabs->comp.water.mass_num = MOD_MASS_NUM;
        slabs->comp.water.t_cs = MOD_TCS;
        slabs->comp.water.el_cs = MOD_ESCS;
        slabs->comp.water.inl_cs = MOD_ISCS;
        slabs->comp.water.abs_cs = MOD_ACS;
    }

}

void init_fission_neutrons(const neutron *const incident_n, neutron *const fission_n) {

    // Set the position of the emitted particles to be the same of the incident particle:
    fission_n->x = incident_n->x;
    fission_n->y = incident_n->y;
    fission_n->z = incident_n->z;

    //Select the direction of the emitted particles from a isotropic distribution:
    iso_direction(fission_n);

    //Set the energy and of the emitted particles:
    fission_n->energy = 1.50000e6; /* 1.50 MeV */

}

bool check_path(neutron *const particle, const medium *const slabs) {

    double R_O, R_R, disc;

    // Compute the dot product between the direction vector of the particle and the vector perpendicular to the slab that it is located:
    R_O = (particle->x * particle->u) + (particle->y * particle->v) + (particle->z * particle->w);
    switch (slabs->tag) {
        case FUEL:
            R_R = SQUARE(particle->x) + SQUARE(particle->y) + SQUARE(particle->z) - SQUARE(FUEL_SLAB_LIMIT/2);
            break;
        case MODERATOR:
            R_R = SQUARE(particle->x) + SQUARE(particle->y) + SQUARE(particle->z) - SQUARE((MODERATOR_SLAB_LIMIT - FUEL_SLAB_LIMIT)/2);
            break;
    }    
    disc = SQUARE(R_O) - R_R;

    // Check if the particle still inside the same slab:
    if(disc < 0) 
        return false;

    particle->in = -R_O - sqrt(disc);
    particle->out = -R_O + sqrt(disc);

    // Check if the particle is leaving or entering in the region:
    if(particle->out < 0) 
        return false;

    return true;
    
}

void select_path(neutron *const particle, const medium *const slabs) {

    double path_length;

    // Select the path length of the particle from the total cross-section of the region that it is located:  
    switch(slabs->tag) {
        case FUEL:
            path_length = -log(mersenne())/slabs->comp.u_diox.t_cs;
            break;
        case MODERATOR:
            path_length = -log(mersenne())/slabs->comp.water.t_cs;
            break;
    }

    // Check if the particle is leaving or re-entering the region:
    if(check_path(particle, slabs))
        path_length += (particle->out - particle->in); 

    // Compute the new position after the ramdom walk:
    particle->x += particle->u * path_length;
    particle->y += particle->v * path_length;
    particle->z += particle->w * path_length;

}

short select_interaction(neutron *const particle, const medium *const slabs, unsigned long *const n) {

    short additional_n = 0;

    double t, mass_num;
    double els_p, inls_p, abs_p, fission_p;
    double R;

    // Compute the probability of occurrence of each possible event from the partial cross-sections:
    switch(slabs->tag) {
        case FUEL:
            t = slabs->comp.u_diox.t_cs;
            els_p = slabs->comp.u_diox.el_cs/t;
            inls_p = slabs->comp.u_diox.inl_cs/t;
            abs_p = slabs->comp.u_diox.abs_cs/t;               
            fission_p = slabs->comp.u_diox.fission_cs/t;
            mass_num = slabs->comp.u_diox.mass_num;
            break;
        case MODERATOR:
            t = slabs->comp.water.t_cs;
            els_p = slabs->comp.water.el_cs/t;
            inls_p = slabs->comp.water.inl_cs/t;
            abs_p = slabs->comp.water.abs_cs/t;
            fission_p = 0.00000;
            mass_num = slabs->comp.water.mass_num;
            break;
    }

    // Select the interaction an its outcome:
    R = mersenne();
    if(R < els_p) {
        iso_scattering(particle, mass_num, false);
    }
    else if(R < els_p + inls_p) {
        iso_scattering(particle, mass_num, true);
    }
    else if(R < els_p + inls_p + abs_p) {
        additional_n = -1;
    }
    else if(R <= els_p + inls_p + abs_p + fission_p) {  /* Never is selected if the particle is in a region filled with the moderator. */
        R = mersenne();
        if(R < ABS(2 - slabs->comp.u_diox.n_fission)) {
            additional_n = 3;
            *n += additional_n;
        }
        else if(R <= ABS(2 - slabs->comp.u_diox.n_fission) + ABS(3 - slabs->comp.u_diox.n_fission)) {
            additional_n = 2;
            *n += additional_n;
        }
    }

    return additional_n;

}

void direction_after_colission(neutron *const particle, const double cos_theta) {

    double c1, c2, c3, c4, c5, sin_delta, cos_delta;
    double w_mod = ABS(particle->w);
    
    do {
        c1 = 2 * mersenne() - 1;
        c2 = 2 * mersenne() - 1;
        c3 = SQUARE(c1) + SQUARE(c2);
    } while(c3 - 1 > 0);

    c4 = sqrt(c3);
    sin_delta = c1/c4;
    cos_delta = c2/c4;
    c5 = sqrt(1 - SQUARE(cos_theta));
    
    if(1 - w_mod < 0.0) {
        particle->u = c5 * cos_delta;
        particle->v = c5 * sin_delta;
        particle->w = cos_theta * (particle->w/w_mod);
    }
    else {
        double c6 = sqrt(1 - SQUARE(particle->w));
        particle->u = particle->u * cos_theta + ((c5/c6) * (particle->u * particle->v * cos_delta - particle->v * sin_delta));
        particle->v = particle->v * cos_theta + ((c5/c6) * (particle->u * particle->v * cos_delta + particle->u * sin_delta));
        particle->w = (particle->w * cos_theta) - (c5 * c6 *cos_delta);
    }

}

void iso_scattering(neutron *const particle, const double target_mass, bool inelastic) {

    // Selection of the cosine of the angle of collision in the laboratory ref. from the angle in the center-of-mass ref.
    double cos_theta_cm = 1 - 2 * mersenne();
    double A = target_mass/N_MASS_NUM;
    double h = A * cos_theta_cm;
    double alfa2 = SQUARE((A-1)/(A+1));
    double cos_theta = (1 + h)/sqrt(1 + SQUARE(A) + 2 * h);

    // Compute the new direction vector of the particle:
    direction_after_colission(particle, cos_theta);

    // Selection of the new energy.
    particle->energy = (particle->energy * ((1 - alfa2) * cos_theta_cm + 1 + alfa2))/2;

    if(inelastic) {
        if(target_mass == FUEL_MASS_NUM) particle->energy *= 0.65;
        if(target_mass == MOD_MASS_NUM) particle->energy *= 0.40;
    }

}

void iso_direction(neutron *const particle) {

    // Select of the azimuth scattering angle:
    double cos_theta = 2.0 * mersenne() - 1.0;
    double phi       = 2 * M_PI * mersenne();
	double sin_theta = sqrt(1.0 - SQUARE(cos_theta));
    
    // Compute the new direction vector:
    particle->u = cos(phi) * sin_theta;  
	particle->v = sin(phi) * sin_theta;  
	particle->w = cos_theta;

}

void punctual_source(neutron *particle, const unsigned long *const h) {

    // Initial position and energy:
    particle[*h].x = 0.0;
    particle[*h].y = 0.0;
    particle[*h].z = 0.0;
    particle[*h].energy = 5.00000e5; /* 0.500 MeV */

    // Selection of the initial direction:
    iso_direction(&particle[*h]);

}