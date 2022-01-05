#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "experiment.h"
#include "report.h"

void print_report(const unsigned long n, const unsigned long r, const result *const results) {

    printf("\n\nResults for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n\n", r, n);

    puts("Average number of collisions per layer:\n");
	printf("%-12s%-s\n", "Layer", "Collisions");	
	for(int l = 0; l != LAYERS; l++)
	    printf("%-11d%13.5e %c %-13.5e\n", (l+1), results->coll_mean[l], (char)241, results->coll_dev[l]);

	puts("\n\nAverage number of escapes per layer:\n");
	printf("%-12s%-s\n", "Layer", "Escapes");	
	for(int l = 0; l != LAYERS; l++)
	    printf("%-11d%13.5e %c %-13.5e\n", (l+1), results->esc_mean[l], (char)241, results->esc_dev[l]);

	puts("\n\nAverage total flux vs. radius:\n");
    printf("%-12s%-s\n", "Radius", "Total flux");
	for(int l = 0; l != LAYERS; l++)
	    printf("%-11.2f%13.5e %c %-13.5e\n", ((double)(l+1)*(RADIUS/LAYERS)), results->flux_mean[l], (char)241, results->flux_dev[l]);

	puts("\n");

}

void write_report(const unsigned long n, const unsigned long r, const result *const results) {

	FILE *report = fopen("report.txt", "a");
	if(report == NULL) {
		fprintf(stderr, "Failed to create report file.");
		exit(1);
	}
 
	fprintf(report, "Results for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n", r, n);

    fputs("\nAverage number of collisions per layer:\n", report);
	fprintf(report, "%-9s%-s\n", "Layer", "Collisions");	
	for(int l = 0; l != LAYERS; l++)
	    fprintf(report, "%-8d%13.5e +/- %-13.5e\n", (l+1), results->coll_mean[l],results->coll_dev[l]);

	fputs("\nAverage number of collisions per layer:\n", report);
	fprintf(report, "%-9s%-s\n", "Layer", "Escapes");	
	for(int l = 0; l != LAYERS; l++)
	    fprintf(report, "%-8d%13.5e +/- %-13.5e\n", (l+1), results->esc_mean[l],results->esc_dev[l]);

	fputs("\nAverage total flux vs. radius:\n", report);
    fprintf(report, "%-9s%-s\n", "Radius", "Total flux");
	for(int l = 0; l != LAYERS; l++)
	    fprintf(report, "%-8.2f%13.5e +/- %-13.5e\n", ((double)(l+1)*(RADIUS/LAYERS)), results->flux_mean[l], results->flux_dev[l]);

	fputs("\n\n", report);

}

void csv_report(const unsigned long n, const unsigned long r, const result *const results) {

	char option, sep;
	puts("Create new report? (s/n)");
	do {
	option = tolower(getchar());
	} while (option != 's' && option != 'n');

	if(option == 's') {

		puts("Enter the separator character: (comma ',' or semicolon ';')");
		do {
			sep = getchar();
		} while(sep != ',' && sep != ';');

		FILE *report = fopen("report.csv", "w");
		if(report == NULL) {
			fprintf(stderr, "Failed to create report file.");
			exit(1);
		}

		puts("Write average number of collisions per layer (c), average number of escapes per layer (e) or total flux vs. radius (f)?");
		do {
		option = tolower(getchar());
		} while (option != 'c' && option != 'f' && option != 'e');

		write_csv(report, option, sep, true, results);

		close_report(report);

	}

	else if(option == 'n') {

		FILE *report = fopen("report.csv", "r");
		if(report == NULL) {
			fprintf(stderr, "Failed to open report file.");
			exit(1);
		}

		sep = check_separator(report);
		assert(sep != 0);

		report = freopen("report.csv", "a", report);
		if(report == NULL) {
			fprintf(stderr, "Failed to open report file.");
			exit(1);
		}

		puts("Write average number of collisions per layer (c), average number of escapes per layer (e) or total flux vs. radius (f)?");
		do {
		option = tolower(getchar());
		} while (option != 'c' && option != 'f');

		write_csv(report, option, sep, false, results);

		close_report(report);

	}

}

char check_separator(FILE *const report) {

	int c;
	while ((c = fgetc(report)) != '\n') {
		if(ispunct(c)) 
			return c;
	}

	return 0;

}

void write_csv(FILE *const report, char option, char sep, bool new_report, const result *const results) {

	switch(option) {

		case 'c':
			if(new_report) 
				fprintf(report, "Layer%cCollisions mean%cStd. deviation\n", sep, sep);
			for(int l = 0; l != LAYERS; l++)
    			fprintf(report, "%d%c%.5e%c%.5e\n", (l+1), sep, results->coll_mean[l], sep, results->coll_dev[l]);
			break;

		case 'f':
			if(new_report)
				fprintf(report, "Radius%cTotal flux mean%cStd. deviation\n", sep, sep);
			for(int l = 0; l != LAYERS; l++)
				fprintf(report, "%.5f%c%.5e%c%.5e\n", ((double)(l+1)*(RADIUS/LAYERS)), sep, results->flux_mean[l], sep, results->flux_dev[l]);
			break;
		case 'e':
			if(new_report) 
				fprintf(report, "Layer%cEscapes mean%cStd. deviation\n", sep, sep);
			for(int l = 0; l != LAYERS; l++)
    			fprintf(report, "%d%c%.5e%c%.5e\n", (l+1), sep, results->esc_mean[l], sep, results->esc_dev[l]);
			break;

	}

}

void close_report(FILE *const report) {

	int check = fclose(report);
	if(check != 0) {
		fprintf(stderr, "Failed to close the report file, data may be corrupted.\n");
		exit(1);
	}

}