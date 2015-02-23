//**************************************************************
//  (C) written and directecd by Robert Kloefkorn
//**************************************************************

#ifndef __PRINTHELP_CC__
#define __PRINTHELP_CC__

static void print_help(const char *funcName)
{
  printf("Alternative 1:\n");
  printf("usage: %s paramfile:paramfile <i_start> <i_end>", funcName);
  printf("%s reads a sequence of grids with discretefunctions\n",
   funcName);
  printf("      and displays all data with GRAPE\n");
  printf("      The parameters fem.prefix and fem.io.datafileprefix \n");
  printf("      are read from the parameterfile. \n");

  printf("Alternative 2:\n");
  printf("usage: %s <i_start> <i_end> -i increment", funcName);
  printf(" [-h] [-help] [-p path] ");
  printf("[-m grid] [-s df] [-v vecdf] [[-s drv] [-v drdv]] [-replay manager.replay file]\n");

  printf("%s reads a sequence of grids with discretefunctions\n",
   funcName);
  printf("      and displays all data with GRAPE\n");
  printf("options:\n");
  printf("   -h or -help: display this help\n");
  printf("   -f: use one fixed grid (from non adaptive simulations)\n");
  printf("   -p path: read data from path \n");
  printf("   -m grid: basename of data files \n");
  printf("   -s df:  read discrete function with basename 'df'\n");

  printf("Example\n");
  printf("%s 0 10 -i 5 -s u_h \n",funcName);
  printf("  reads grid grid0000000000 with scalar function u_h0000000000 and\n");
  printf("  then grid0000000005 with u_h0000000005 and\n");
  printf("  and finally grid0000000010 with u_h0000000010\n");
  exit(0);
}
#endif
