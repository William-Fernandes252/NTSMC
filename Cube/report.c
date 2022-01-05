#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "experiment.h"
#include "report.h"

void 
print_report(const unsigned long n, const unsigned long r, 
			 const result *const results) {

    printf("\n\nResults for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n", r, n);

	printf("\n\t\t%16s%17s\n\n", "Mean", "Std. dev.");
	printf("Moderations:\t%16.5e%17.5e\n\n", results->mod_m, results->mod_d);
	printf("Absorptions:\t%16.5e%17.5e\n\n", results->abs_m, results->abs_d);
	printf("Reflections:\t%16.5e%17.5e\n\n", results->ref_m, results->ref_d);
	printf("Detections:\t%16.5e%17.5e\n\n", results->det_m, results->det_d);
	printf("Escapes:\t%16.5e%17.5e\n\n", results->esc_m, results->esc_d);
	printf("Rel. mod.:\t%16.5e%17.5e\n\n", 
		   results->mod_frac_m, results->mod_frac_d);
	printf("Rel. abs.:\t%16.5e%17.5e\n\n", 
		   results->abs_frac_m, results->abs_frac_d);
	printf("Rel. ref.:\t%16.5e%17.5e\n\n", 
	   	   results->ref_frac_m, results->ref_frac_d);
	printf("Rel. det.:\t%16.5e%17.5e\n\n", 
		   results->det_frac_m, results->det_frac_d);
	printf("Rel. esc.:\t%16.5e%17.5e\n\n", 
	  	   results->esc_frac_m, results->esc_frac_d);
	printf("Histories:\t%16.5e%17.5e\n\n", results->hist_m, results->hist_d);
	printf("Steps:\t\t%16.5e%17.5e\n\n", results->stp_m, results->stp_d);
	printf("Mean free path:\t%16.5e%17.5e\n\n\n", results->frp_m, results->frp_d);

}

void 
write_report(const unsigned long n, const unsigned long r, 
			 const result *const results) {

	FILE *report = fopen("report.txt", "a");
	if(report == NULL) {
		fprintf(stderr, "Failed to create report file.");
		exit(1);
	}

	fprintf(report, "Results for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n", r, n);

	fprintf(report, "\n\t\t\t\t%16s%17s\n\n", "Mean", "Std. dev.");
	fprintf(report, "Moderations:\t%16.5e%17.5e\n\n", 
			results->mod_m, results->mod_d);
	fprintf(report, "Absorptions:\t%16.5e%17.5e\n\n", 
			results->abs_m, results->abs_d);
	fprintf(report, "Reflections:\t%16.5e%17.5e\n\n", 
			results->ref_m, results->ref_d);
	fprintf(report, "Detections:\t\t%16.5e%17.5e\n\n", 
			results->det_m, results->det_d);
	fprintf(report, "Escapes:\t\t%16.5e%17.5e\n\n", 
			results->esc_m, results->esc_d);
	fprintf(report, "Rel. mod.:\t\t%16.5e%17.5e\n\n", 
			results->mod_frac_m, results->mod_frac_d);
	fprintf(report, "Rel. abs.:\t\t%16.5e%17.5e\n\n", 
			results->abs_frac_m, results->abs_frac_d);
	fprintf(report, "Rel. ref.:\t\t%16.5e%17.5e\n\n", 
			results->ref_frac_m, results->ref_frac_d);
	fprintf(report, "Rel. det.:\t\t%16.5e%17.5e\n\n", 
			results->det_frac_m, results->det_frac_d);
	fprintf(report, "Rel. esc.:\t\t%16.5e%17.5e\n\n", 
			results->esc_frac_m, results->esc_frac_d);
	fprintf(report, "Histories:\t\t%16.5e%17.5e\n\n", 
			results->hist_m, results->hist_d);
	fprintf(report, "Steps:\t\t\t%16.5e%17.5e\n\n", 
			results->stp_m, results->stp_d);
	fprintf(report, "Mean free path:\t%16.5e%17.5e\n\n\n", 
			results->frp_m, results->frp_d);

	close_report(report);

}

void 
csv_report(const unsigned long n, const unsigned long r, 
		   const result *const results) {

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

		puts("Write absolute (a) or relative (r) quantities?");
		do {
		option = tolower(getchar());
		} while (option != 'a' && option != 'r');

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

		puts("Write absolute (a) or relative (r) quantities?");
		do {
		option = tolower(getchar());
		} while (option != 'a' && option != 'r');

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

void 
write_csv(FILE *const report, char option, char sep, bool new_report, 
		  const result *const results) {

	switch(option) {

		case 'a':
			if(new_report) 
				fprintf(report,"Moderations%cAbsorptions%cReflections%cDetections%cEscapes%cHistories%cSteps%cMean free path\n", sep, sep, sep, sep, sep, sep, sep);
			fprintf(report, "%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", results->mod_m, sep, results->abs_m, sep, results->ref_m, sep, results->det_m, sep, results->esc_m, sep, results->hist_m, sep, results->stp_m, sep, results->frp_m);
			fprintf(report, "%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", results->mod_d, sep, results->abs_d, sep, results->ref_d, sep, results->det_d, sep, results->esc_d, sep, results->hist_d, sep, results->stp_d, sep, results->frp_d);
			break;

		case 'r':
			if(new_report)
				fprintf(report,"Moderations%cAbsorptions%cReflections%cDetections%cEscapes", sep, sep, sep, sep);
			fprintf(report, "%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", results->mod_m, sep, results->abs_m, sep, results->ref_m, sep, results->det_m, sep, results->esc_m, sep, results->stp_m, sep, results->frp_m);
			fprintf(report, "%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", results->mod_d, sep, results->abs_d, sep, results->ref_d, sep, results->det_d, sep, results->esc_d);
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