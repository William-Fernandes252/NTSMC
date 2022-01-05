#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "experiment.h"
#include "report.h"

void print_report(const unsigned long n, const unsigned long r, const result *const results) {

    printf("\n\nResults for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n", r, n);

	printf("\n\t\t%13s%17s\n\n", "Mean", "Std. dev.");
	printf("Moderations:\t%13.5e%17.5e\n\n", results->mod_m, results->mod_d);
	printf("Absorptions:\t%13.5e%17.5e\n\n", results->abs_m, results->ref_d);
	printf("Reflections:\t%13.5e%17.5e\n\n", results->ref_m, results->ref_d);
	printf("Escapes:\t%13.5e%17.5e\n\n", results->esc_m, results->esc_d);
	printf("Rel. mod.:\t%13.5e%17.5e\n\n", results->mod_frac_m, results->mod_frac_d);
	printf("Rel. abs.:\t%13.5e%17.5e\n\n", results->abs_frac_m, results->abs_frac_d);
	printf("Rel. ref.:\t%13.5e%17.5e\n\n", results->ref_frac_m, results->ref_frac_d);
	printf("Rel. esc.:\t%13.5e%17.5e\n\n", results->esc_frac_m, results->esc_frac_d);
	printf("Histories:\t%13.5e%17.5e\n\n\n", results->hist_m, results->hist_d);

}

void write_report(const unsigned long n, const unsigned long r, const result *const results) {

	FILE *report = fopen("report.txt", "a");
	if(report == NULL) {
		fprintf(stderr, "Failed to create report file.");
		exit(1);
	}

	fprintf(report, "Results for %lu repetitions of the experiment with %lu initial particles from the source:\n------------------------------------------------------------------------------------------\n", r, n);

	fprintf(report, "\n\t\t\t\t%13s%17s\n\n", "Mean", "Std. dev.");
	fprintf(report, "Moderations:\t%13.5e%17.5e\n\n", results->mod_m, results->mod_d);
	fprintf(report, "Absorptions:\t%13.5e%17.5e\n\n", results->abs_m, results->abs_d);
	fprintf(report, "Reflections:\t%13.5e%17.5e\n\n", results->ref_m, results->ref_d);
	fprintf(report, "Escapes:\t\t%13.5e%17.5e\n\n", results->esc_m, results->esc_d);
	fprintf(report, "Rel. mod.:\t\t%13.5e%17.5e\n\n", results->mod_frac_m, results->mod_frac_d);
	fprintf(report, "Rel. abs.:\t\t%13.5e%17.5e\n\n", results->abs_frac_m, results->abs_frac_d);
	fprintf(report, "Rel. ref.:\t\t%13.5e%17.5e\n\n", results->ref_frac_m, results->ref_frac_d);
	fprintf(report, "Rel. esc.:\t\t%13.5e%17.5e\n\n", results->esc_frac_m, results->esc_frac_d);
	fprintf(report, "Histories:\t\t%13.5e%17.5e\n\n\n", results->hist_m, results->hist_d);

	close_report(report);

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

		fprintf(report, "Histories%cRepetitions%cMod. mean%cMod. std.%cAbs. mean%cAbs. std.%cRef. mean%cRef. std.%cEsc. mean%cEsc. std.%cRelative mod. mean%cRelative mod. std.%cRelative abs. mean%cRelative abs. std.%cRelative ref. mean%cRelative ref. std.%cRelative esc. mean%cRelative esc. std.\n", sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep);
		fprintf(report, "%lu%c%lu%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", n, sep, r, sep, results->mod_m, sep, results->mod_d, sep, results->abs_m, sep, results->abs_d, sep, results->ref_m, sep, results->ref_d, sep, results->esc_m, sep, results->esc_d, sep, results->mod_frac_m, sep, results->mod_frac_d, sep, results->abs_frac_m, sep, results->abs_frac_d, sep, results->ref_frac_m, sep, results->ref_frac_d, sep, results->esc_frac_m, sep, results->esc_frac_d);

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

		fprintf(report, "%lu%c%lu%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e%c%.5e\n", n, sep, r, sep, results->mod_m, sep, results->mod_d, sep, results->abs_m, sep, results->abs_d, sep, results->ref_m, sep, results->ref_d, sep, results->esc_m, sep, results->esc_d, sep, results->mod_frac_m, sep, results->mod_frac_d, sep, results->abs_frac_m, sep, results->abs_frac_d, sep, results->ref_frac_m, sep, results->ref_frac_d, sep, results->esc_frac_m, sep, results->esc_frac_d);

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

void close_report(FILE *const report) {

	int check = fclose(report);
	if(check != 0) {
		fprintf(stderr, "Failed to close the report file, data may be corrupted.\n");
		exit(1);
	}

}