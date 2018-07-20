//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../libraries_common/config_helpers.h"
#include "../monodomain/config/purkinje_config.h"
#include "../utils/logfile_utils.h"
#include <assert.h>
#include <time.h>

#ifdef _MSC_VER
    #include <process.h>
    #define getpid _getpid
#else
    #include <unistd.h>
#endif

SET_SPATIAL_PURKINJE (initialize_purkinje_with_custom_mesh) 
{

    char *network_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (network_file,config->config_data.config,"network_file");

    print_to_stdout_and_file("Loading Purkinje Network: %s\n",config->domain_name);
    set_custom_purkinje_network(the_grid, network_file);

    free (network_file);

    print_to_stdout_and_file("Leaving program ...\n");
    exit(EXIT_FAILURE);
}

// TO DO: Build some benchmark Purkinje network for quick tests ...