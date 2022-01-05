#undef __STRICT_ANSI__

#include <math.h>
#include <stdbool.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"

void init_cross_sections(const neutron *const particle, medium *const cube) {

    // Compute the cross-sections from the particle atual position in the medium:
    if(particle->x <= FUEL_X2) {
        cube->tag = FUEL;
        cube->comp.u_diox.mass_num = FUEL_MASS_NUM;
        cube->comp.u_diox.t_cs = FUEL_TCS;
        cube->comp.u_diox.el_cs = FUEL_ESCS;
        cube->comp.u_diox.inl_cs = FUEL_ISCS;
        cube->comp.u_diox.fission_cs = FUEL_FCS;
        cube->comp.u_diox.abs_cs = FUEL_ACS;
        cube->comp.u_diox.n_fission = FISSION_NEUTRONS;
    }
    else if(particle->x <= MODERATOR_X2) {
        cube->tag = MODERATOR;
        cube->comp.water.mass_num = MOD_MASS_NUM;
        cube->comp.water.t_cs = MOD_TCS;
        cube->comp.water.el_cs = MOD_ESCS;
        cube->comp.water.inl_cs = MOD_ISCS;
        cube->comp.water.abs_cs = MOD_ACS;
    }

}

void init_fission_neutrons(neutron *const incident_n, neutron *const fission_n, const medium *const cube) {

    // Set the position of the emitted particles to be the same of the incident particle:
    fission_n->x = incident_n->x;
    fission_n->y = incident_n->y;
    fission_n->z = incident_n->z;

    // Select the direction of the emitted particles from a isotropic distribution:
    iso_direction(fission_n);

    // Set the energy and the weight of the emitted particles:
    fission_n->energy = 1.50000e6; /* 1.50 MeV */
    fission_n->weight = incident_n->weight * cube->comp.u_diox.n_fission;

    // Reduce the weight of the incident particle:
    incident_n->weight = (cube->comp.u_diox.fission_cs - cube->comp.u_diox.fission_cs)/cube->comp.u_diox.fission_cs;

}

double select_path(neutron *const particle, const medium *const cube) {

    double path_length;
    
    // Select the path length of the particle from the total cross-section of the region that it is located:  
    switch(cube->tag) {
        case FUEL:
            path_length = -log(mersenne())/cube->comp.u_diox.t_cs;
            break;
        case MODERATOR:
            path_length = -log(mersenne())/cube->comp.water.t_cs;
            break;
    }
    
    // Compute the new position after the ramdom walk:
    particle->x += particle->u * path_length;
    particle->y += particle->v * path_length;
    particle->z += particle->w * path_length;
    
    return path_length;

}

bool select_interaction(neutron *const particle, const medium *const cube) {

    double t, mass_num, wr_factor;
    double els_p, inls_p, abs_p, fission_p;
    double R;

    // Compute the probability of occurrence of each possible event from the partial cross-sections:
    switch(cube->tag) {
        case FUEL:
            t = cube->comp.u_diox.t_cs;
            els_p = cube->comp.u_diox.el_cs/t;
            inls_p = cube->comp.u_diox.inl_cs/t;
            abs_p = cube->comp.u_diox.abs_cs/t;               
            fission_p = cube->comp.u_diox.fission_cs/t;
            mass_num = cube->comp.u_diox.mass_num;
            wr_factor = (t - cube->comp.u_diox.abs_cs)/t;
            break;
        case MODERATOR:
            t = cube->comp.water.t_cs;
            els_p = cube->comp.water.el_cs/t;
            inls_p = cube->comp.water.inl_cs/t;
            abs_p = cube->comp.water.abs_cs/t;
            fission_p = 0.00000;
            mass_num = cube->comp.u_diox.mass_num;
            wr_factor = (t - cube->comp.water.abs_cs)/t;
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
        particle->weight *= wr_factor;
    }
    else if(R <= els_p + inls_p + abs_p + fission_p) { /* Never is selected if the particle is in region filled with the moderator material. */
        return true;
    }

    return false;

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

bool russian_roulette(neutron *const particle) {

    // Compute the survival probability of the particle:
    double survival_p = particle->weight/INITIAL_WEIGHT;

    // Play the Russian Roulette:
    double R = mersenne();
    if(R < survival_p) {
        particle->weight = INITIAL_WEIGHT;
        return false;
    }
    return true;

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
    particle[*h].x = (FUEL_X2 + FUEL_X1)/2.0;
    particle[*h].y = (Y2 + Y1)/2.0;
    particle[*h].z = (Z2 + Z1)/2.0;
    particle[*h].energy = 5.00000e5; /* 0.500 MeV */
    particle[*h].weight = INITIAL_WEIGHT;

    // Selection of the initial direction:
    iso_direction(&particle[*h]);
    
}