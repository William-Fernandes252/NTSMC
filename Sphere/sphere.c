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
Sphere experiment:\n\
Punctual source of neutrons located in the center of a sphere\n\
divided into determined number of layers, where the inner ones\n\
are filed with a fissionable material (fuel), and the outer\n\
ones are filed with the a non-fissionable material (moderator).\n\n\
Usage: experiment n r o\n\
n - number of neutrons emitted bt the source;\n\
r - number of repetitions of the experiment;\n\
o - report format option:\n\
\t'p' prints a report in the terminal;\n\
\t'w' writes a report into the file report.txt;\n\
\t'c' writes the results in the csv file report.csv\n\
\t'a' both writes the results in csv format in the file and prints the report \
in the terminal.\n\n\
";

int main(int argc, char **argv) {

    /* Seed for the generation of the pseudo-random numbers sequence. */
    seed_mersenne(1363891763); 

    if (argc < 4) {
        fprintf(stderr, "%s",usage);
        exit(1);
    }
    
    /* Number of particles emitted by the font. */    
    unsigned long n = atol(argv[1]);
    /* number of repetitions of the experiment. */   
    unsigned long r = atol(argv[2]);   
    char output = argv[3][0];

    result results;
    for(unsigned int l = 0; l < LAYERS; l++) {
        results.coll_mean[l] = 0.00000;
        results.coll_dev[l] = 0.00000;
        results.esc_mean[l] = 0.00000;
        results.esc_dev[l] = 0.00000;
        results.flux_mean[l] = 0.00000;
        results.flux_dev[l] = 0.00000;
    }

    sphere(n, r, &results);

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