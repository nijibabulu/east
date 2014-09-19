#include "output.h"

char * east_formats[] = {
  "blast",
  "stitch",
  NULL
};

int 
output_check_format_string(const char *s, int die)
{
  int i;

  for(i = 0; east_formats[i] != NULL; i ++) 
    if(east_formats[i][0] == s[0] || !strcmp(east_formats[i],s))
      return 1;
  
  if(die) {
    fprintf(stderr, "Unrecognized format %s\n", s);
    exit(1);
  }
  else
    return 0;
}

void
output_print_format_strings(FILE *f, int indent)
{
  int i;

  for(i = 0; east_formats[i] != NULL; i ++) 
    fprintf(f,"%*c[%c]%s\n", indent, ' ', east_formats[i][0], &east_formats[i][1]);
}
