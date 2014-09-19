#ifndef OUTPUT_H
#define OUTPUT_H

#include "traceback.h"

/* maybe more args if we get other options which prohibit certain formats */
int output_check_format_string(const char *, int die); 
void output_print_format_strings(FILE *,int indent);
void output_tb(tb_t *, const char *);

extern char * east_formats[];

#endif /* OUTPUT_H */
