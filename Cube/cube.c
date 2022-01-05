#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>
#include "mersenne.h"
#include "simulations.h"
#include "experiment.h"
#include "report.h"

char usage[] = "\
\nNeutron Transport Simulation with Monte Carlo\n\
---------------------------------------------\n\n\
Cube experiment:\n\
Punctual source of neutrons located in the center of a slab filled with\n\
a fissionable material (fuel), which is unitted with a cube filled with\n\
a non-fissionable material (moderator). There is no barriers between the\n\
regions, so the particles can travel from one to another with no limitations.\n\
In the center of the cube there is a smaller cube, which is the region\n\
monitored by a detector.\n\n\
Usage: experiment n r o\n\
n - number of neutrons emitted bt the source;\n\
r - number of repetitions of the experiment;\n\
o - report format option:\n\
\t'p' prints a report in the terminal;\n\
\t'w' writes a report into the file report.txt;\n\
\t'c' writes the results in the csv file report.csv;\n\
\t'a' both writes the results in csv format in the file report.csv and prints \
the report in the terminal.\n\n\
";

int main(int argc, char **argv) {

    /* Seed for the generation of the pseudo-random numbers sequence. */
    seed_mersenne(1363891763);

    if (argc < 4) {
        fprintf(stderr, "%s", usage);
        exit(1);
    }

    /* Number of particles emitted by the font. */
    unsigned long n = atol(argv[1]);
    /* Number of repetitions of the experiment. */
    unsigned long r = atol(argv[2]);
    char output = argv[3][0];

    result results = {
        .mod_m = 0.00000,
        .abs_m = 0.00000,
        .ref_m = 0.00000,
        .det_m = 0.00000,
        .esc_m = 0.00000,
        .stp_m = 0.00000,
        .frp_m = 0.00000,
        .mod_d = 0.00000,
        .abs_d = 0.00000,
        .ref_d = 0.00000,
        .det_d = 0.00000,
        .esc_d = 0.00000,
        .stp_d = 0.00000,
        .frp_d = 0.00000,
        .mod_frac_m = 0.00000,
        .abs_frac_m = 0.00000,
        .ref_frac_m = 0.00000,
        .det_frac_m = 0.00000,
        .esc_frac_m = 0.00000,
        .mod_frac_d = 0.00000,
        .abs_frac_d = 0.00000,
        .ref_frac_d = 0.00000,
        .det_frac_d = 0.00000,
        .esc_frac_d = 0.00000,
    };

    cube(n, r, &results);

    switch(output) {
        case 'p': print_report(n, r, &results); break;
        case 'w': write_report(n, r, &results); break;
        case 'c': csv_report(n, r, &results); break;
        case 'a': 
            print_report(n, r, &results);
            csv_report(n, r, &results);
            break;
        default: puts("Invalid output."); exit(0);
    }

    return EXIT_SUCCESS;

}