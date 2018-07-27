//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"
#include "../utils/logfile_utils.h"
#include "../utils/utils.h"
#include <float.h>
#include <math.h>
#include <string.h>
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

    struct graph *purkinje = the_grid->the_purkinje_network;

    set_purkinje_network_from_file(purkinje,file_name);
    
    /*
    FILE *file = fopen (file_name, "r");

    int N, E;
    fscanf(file,"%d %d",&N,&E);

    if (!file) 
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", file_name);
        exit (0);
    }

    fclose(file);
    */

    // deallocate memory

}

void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name) 
{
    struct graph *skeleton_network = new_graph();

    // TO DO: Build the graph from the file
    //read_purkinje_network_from_file(file_name,&points,&branches,&N,&E);
    build_skeleton_purkinje(file_name,skeleton_network);
    
    // TO DO: Finish implementing this function
    //build_mesh_purkinje(the_purkinje_network,skeleton_network);
    
    //print_graph(skeleton_network);

    //double pos1[3] = {0.0, 0.0, 0.0};
    //double pos2[3] = {1.0, 1.0, 1.0};
    //double pos3[3] = {2.0, 2.0, 2.0};

    //insert_node_graph(the_purkinje_network,pos1);
    //insert_node_graph(the_purkinje_network,pos2);
    //insert_node_graph(the_purkinje_network,pos3);

    //print_graph(the_purkinje_network);
    
}

// TO DO: Find a way to build the purkinje network mesh directly without constructing the skeleton graph
void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network)
{

    FILE *file = fopen(filename,"r");
    if (!file)
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", filename);
        exit (EXIT_FAILURE);
    }

    int N;
    char str[100];

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    if (!fscanf(file,"%d",&N))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    if (!fscanf(file,"%s",str))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }

    // Read points
    for (int i = 0; i < N; i++)
    {
        double pos[3];
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        } 
        insert_node_graph(skeleton_network,pos);
    }

    // Read edges
    int trash, E;
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(file,"%d %d",&E,&trash))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < E; i++)
    {
        int e[2];
        if (!fscanf(file,"%d %d %d",&trash,&e[0],&e[1]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        insert_edge_graph(skeleton_network,e[0],e[1]);
    }

    fclose(file);
}

void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network)
{
    uint32_t n = skeleton_network->total_nodes;
    uint32_t *map_skeleton_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));
    
    // Construct the first node
    struct node *tmp = skeleton_network->list_nodes;
    double pos[3]; pos[0] = tmp->x; pos[1] = tmp->y; pos[2] = tmp->z;
    insert_node_graph(the_purkinje_network,pos);
    
    // Make a Depth-First-Search to build the mesh of the Purkije network
    depth_first_search(the_purkinje_network,tmp,0);

    free(map_skeleton_to_mesh);
}

void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level)
{
    // TO DO: Include the diameter of the Purkinje branch here ...

    struct edge *v = u->list_edges;
    while (v != NULL)
    {
        // TO DO: Implement this function
        //grow_segment(the_purkinje_network,u,v);
        depth_first_search(the_purkinje_network,v->dest,level+1);
        v = v->next;
    }
}

void read_purkinje_network_from_file (const char *filename, struct point **points, struct branch **branches, int *N, int *E)
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", filename);
        exit (EXIT_FAILURE);
    }

    char str[100];

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    if (!fscanf(file,"%d",N))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    if (!fscanf(file,"%s",str))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }

    // Read points
    *points = (struct point*)malloc(sizeof(struct point)*(*N));
    for (int i = 0; i < (*N); i++)
    {
        double pos[3];
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        } 
        (*points)[i].x = pos[0];
        (*points)[i].y = pos[1];
        (*points)[i].z = pos[2];
    }

    // Read edges
    int trash;
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(file,"%d %d",E,&trash))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    *branches = (struct branch*)malloc(sizeof(struct branch)*(*E));
    for (int i = 0; i < (*E); i++)
    {
        int e[2];
        if (!fscanf(file,"%d %d %d",&trash,&e[0],&e[1]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        (*branches)[i].source = e[0];
        (*branches)[i].destination = e[1];
    }

    fclose(file);

}





