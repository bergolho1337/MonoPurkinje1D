//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"
#include "../utils/logfile_utils.h"
#include "../utils/utils.h"
#include <float.h>
#include <math.h>
#include <time.h>

#ifdef _MSC_VER
#include <process.h>
    #define getpid _getpid
#else
#include <unistd.h>
#endif

// TO DO: Implement this function
// Set a a custom Purkinje network from a file that stores its graph structure
void set_custom_purkinje_network (struct grid *the_grid, const char *file_name) 
{

    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen (file_name, "r");

    int N, E;
    fscanf(file,"%d %d",&N,&E);
    print_to_stdout_and_file("|| N = %d || E = %d ||\n",N,E);

    if (!file) 
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", file_name);
        exit (0);
    }

    fclose(file);

    // deallocate memory

}





