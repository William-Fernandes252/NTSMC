#ifndef REPORT_H
#define REPORT_H

/*
    print_report:   Prints a formatted report in the terminal.
*/
void print_report(const unsigned long n, const unsigned long r, const result *const results);

/*
    write_report:   Writes a formatted report in the file specified.
*/
void write_report(const unsigned long n, const unsigned long r, const result *const results);

/*
    csv_report:     Writes the results in a semi-colon or comma separated csv file.
*/
void csv_report(const unsigned long n, const unsigned long r, const result *const results);

/*
    check_separator:    Look for a separator for a csv file in the first line of of the
                        input file.      
*/
char check_separator(FILE *const report);

/*
    close_report:   Closes the specified file and check if any error happened in this process.
*/
void close_report(FILE *const report);

#endif