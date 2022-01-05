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
Slabs experiment:\n\
Punctual source of neutrons located in the origin under two slabs\n\
with an infinite area, where the first one is filed with a fissionable\n\
material (fuel) and the outer one if filed with a non-fissionable\n\
material (moderator).\n\n\
Usage: experiment n r o\n\
n - number of neutrons emitted bt the source;\n\
r - number of repetitions of the experiment;\n\
o - report format option:\n\
\t'p' prints a report in the terminal;\n\
\t'w' writes a report into the file report.txt;\n\
\t'c' writes the results in the semi-colon or comma separeted csv file report.csv\n\
\t'a' both writes the results in csv format in the file report.csv and prints the report in the terminal.\n\n\
";

int main(int argc, char **argv) {

    seed_mersenne(1363891763); /* Seed for the generation of the pseudo-random numbers sequence. */

    if (argc < 4) {
        fprintf(stderr, "%s", usage);
        exit(1);
    }

    unsigned long n = atol(argv[1]);   /* number of particles emitted by the font. */
    unsigned long r = atol(argv[2]);   /* number of repetitions of the experiment. */
    char output = argv[3][0];

    result results = {
        .mod_m = 0.00,
        .abs_m = 0.00,
        .ref_m = 0.00,
        .esc_m = 0.00,
        .mod_d = 0.00,
        .abs_d = 0.00,
        .ref_d = 0.00,
        .esc_d = 0.00,
        .mod_frac_m = 0.00,
        .abs_frac_m = 0.00,
        .ref_frac_m = 0.00,
        .esc_frac_m = 0.00,
        .mod_frac_d = 0.00,
        .abs_frac_d = 0.00,
        .ref_frac_d = 0.00,
        .esc_frac_d = 0.00,
    };

    slabs(n, r, &results);

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