#ifndef REPORT_H
#define REPORT_H

/*
    print_report:   Prints a formatted report in the terminal.
*/
void 
print_report(const unsigned long n, const unsigned long r, 
             const result *const results);

/*
    write_report:   Writes a formatted report in the file specified.
*/
void 
write_report(const unsigned long n, const unsigned long r, 
             const result *const results);

/*
    csv_report:     Writes the results in a semi-colon or comma separated csv 
                    file.
*/
void 
csv_report(const unsigned long n, const unsigned long r, 
           const result *const results);

/*
    check_separator:    Look for a separator for a csv file in the first line of 
                        the input file.      
*/
char 
check_separator(FILE *const report);

/*
    write_data:     Depending on the chosen option, either writes the absolute 
                    or relative  quantities estimated in the given stream. The 
                    in the first line after the collumn labels the mean of the 
                    quantities is written, and in the second line, the standard
                    deviation. If it is a new file, writes the collumns labels.
*/
void 
write_csv(FILE *const report, char option, char sep, bool new_report, 
          const result *const results);

/*
    close_report:   Closes the specified file and check if any error happened 
                    in this process.
*/
void 
close_report(FILE *const report);

#endif