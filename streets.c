#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"
#include <limits.h>

#define MAX_WAY_NAME_LENGTH 100

//-------------------------------------------Task1-------------------------------------//
struct node {
    int id; 
    double lat; 
    double lon; 
    int num_ways; 
    int *way_ids; 
};

struct way {
    int id; 
    char name[MAX_WAY_NAME_LENGTH]; 
    double maxspeed; 
    bool oneway; 
    int num_nodes; 
    int *node_ids; 
};

struct ssmap {
    struct node *nodes; 
    int num_nodes; 
    struct way *ways; 
    int num_ways; 
};
//-------------------------------------------Task2-------------------------------------//
struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    if (nr_nodes == 0 || nr_ways == 0) {
        return NULL;
    }

    struct ssmap *map = malloc(sizeof(struct ssmap));
    if (map == NULL) {
        return NULL;
    }

    map->nodes = malloc(sizeof(struct node) * nr_nodes);
    if (map->nodes == NULL) {
        free(map);
        return NULL;
    }
    
    map->ways = malloc(sizeof(struct way) * nr_ways);
    if (map->ways == NULL) {
        free(map->nodes);
        free(map);
        return NULL;
    }

    map->num_nodes = nr_nodes;
    map->num_ways = nr_ways;

    return map;
}

bool
ssmap_initialize(struct ssmap * m)
{
    /* TODO: task 2
     * additional initialization code can be added here */
    return true;
}

void
ssmap_destroy(struct ssmap * m)
{
   for (int i = 0; i < m->num_nodes; i++) {
        free(m->nodes[i].way_ids);
    }
    free(m->nodes);
    for (int i = 0; i < m->num_ways; i++) {
        free(m->ways[i].node_ids);
    }
    free(m->ways);
    free(m);
}

struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    if (id < 0 || id >= m->num_ways) {
        return NULL;
    }
    struct way *new_way = &m->ways[id];

    new_way->id = id;
    strncpy(new_way->name, name, MAX_WAY_NAME_LENGTH); 
    new_way->name[MAX_WAY_NAME_LENGTH - 1] = '\0'; 
    new_way->maxspeed = maxspeed;
    new_way->oneway = oneway;
    new_way->num_nodes = num_nodes;

    new_way->node_ids = malloc(sizeof(int) * num_nodes);
    if (new_way->node_ids == NULL) {
        return NULL;
    }

    for (int i = 0; i < num_nodes; i++) {
        new_way->node_ids[i] = node_ids[i];
    }

    return new_way;

}

struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    if (id < 0 || id >= m->num_nodes) {
        return NULL;
    }
    
    struct node *new_node = &m->nodes[id];
    new_node->id = id;
    new_node->lat = lat;
    new_node->lon = lon;
    new_node->num_ways = num_ways;

    new_node->way_ids = malloc(sizeof(int) * num_ways);
    if (new_node->way_ids == NULL) {

        return NULL;
    }

    for (int i = 0; i < num_ways; i++) {
        new_node->way_ids[i] = way_ids[i];
    }

    return new_node;
}
//--------------------------------------------Task3--------------------------------------------//
void
ssmap_print_way(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->num_ways) {
        printf("error: way %d does not exist.\n", id);
        return;
    }
    const struct way *way = &m->ways[id];
    if (way != NULL) {
        printf("Way %d: %s\n", id, way->name);
    } else {
        printf("error: way %d does not exist.\n", id);
    }
}

void
ssmap_print_node(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->num_nodes) {
        printf("error: node %d does not exist.\n", id);
        return;
    }

    const struct node *node = &m->nodes[id];
    if (node != NULL) {
        printf("Node %d: (%.7f, %.7f)\n", id, node->lat, node->lon);
    }else {
        printf("error: node %d does not exist.\n", id);
    }
}

void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    int found = 0;
    
    for (int i = 0; i < m->num_ways; i++) {
        if (strstr(m->ways[i].name, name) != NULL) {
            printf("%d ", m->ways[i].id);
            found = 1; 
        }
    }
    
    if (!found) {
        printf("No ways found with the name containing '%s'\n", name);
    } else {
        printf("\n");
    }
}

void 
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    if (name2 == NULL)
    {
        int nod_num = m->num_nodes;
        
        for (int i = 0;i<nod_num;i++){
            struct node node = m->nodes[i];
            for (int j = 0; j < node.num_ways; j++)
            {
                if(strstr(m->ways[node.way_ids[j]].name, name1) != NULL){
                    printf("%d ",i);
                    break;
                }
            }
        }
    }else
    {   
        int nod_num = m->num_nodes;
        for (int i = 0; i < nod_num; i++)
        {
            struct node node = m->nodes[i];
            bool match_name1 = false;
            bool match_name2 = false;
            bool *match_name1_pt = &match_name1;
            bool *match_name2_pt = &match_name2;
            for (int j=0;j<node.num_ways;j++){
                if ((*match_name1_pt == false) && strstr(m->ways[node.way_ids[j]].name,name1)!=NULL)
                {
                    *match_name1_pt = true;
                }
                else if (strstr(m->ways[node.way_ids[j]].name,name2) != NULL)
                {
                    *match_name2_pt = true;
                }
                if (*match_name1_pt && *match_name2_pt){
                    printf("%d ",i);
                    break;
                }
            }
        }
        printf("\n");
        return;
    }
}

/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}
//----------------------------------Helper Function---------------------------//
bool nodes_are_connected(const struct ssmap *m, int node_a_id, int node_b_id) {
    const struct node *node_a = &m->nodes[node_a_id];
    // const struct node *node_b = &m->nodes[node_b_id];
    if(node_a_id == node_b_id){
        return true;
    }

    // Check all ways associated with node_a
    for (int i = 0; i < node_a->num_ways; i++) {
        const struct way *way = &m->ways[node_a->way_ids[i]];

        // Check if node_b is also in this way
        for (int j = 0; j < way->num_nodes - 1; j++) {
            if ((way->node_ids[j] == node_a_id && way->node_ids[j + 1] == node_b_id) ||
                (!way->oneway && way->node_ids[j] == node_b_id && way->node_ids[j + 1] == node_a_id)) {
                return true;
            }
        }
    }
    // No connecting way found
    return false;
}

bool are_nodes_adjacent_in_way(const struct way *way, int node_a_id, int node_b_id) {
    for (int i = 0; i < way->num_nodes - 1; i++) {
        if (way->node_ids[i] == node_a_id && way->node_ids[i + 1] == node_b_id) {
            return true; 
        }
        // If the way is not one-way, check for adjacency in the reverse direction (b -> a)
        if (!way->oneway && way->node_ids[i + 1] == node_a_id && way->node_ids[i] == node_b_id) {
            return true; 
        }
    }
    return false; 
}

const struct way *find_way_with_nodes(const struct ssmap *m, int node_a_id, int node_b_id) {
    // Retrieve the nodes from the map
    const struct node *node_a = &m->nodes[node_a_id];

    // Check all ways associated with node_a to see if node_b is also part of the way
    for (int i = 0; i < node_a->num_ways; i++) {
        const struct way *way = &m->ways[node_a->way_ids[i]];
        // Check if node_b is in this way
        for (int j = 0; j < way->num_nodes; j++) {
            if (way->node_ids[j] == node_b_id) {
                // Found the way that contains both node_a and node_b
                return way;
            }
        }
    }

    // If no way is found that contains both nodes, return NULL
    return NULL;
}

//--------------------------------Task4-----------------------------------//
double 
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    //Error1
    for (int i = 0; i < size; i++) {
        if (node_ids[i] < 0 || node_ids[i] >= m->num_nodes) {
            printf("error: node %d does not exist.\n", node_ids[i]);
            return -1.0;
        }
    }
    
    //Error2 3 4
    for (int i = 0; i < size - 1; i++) {
        const struct way *connecting_way = find_way_with_nodes(m, node_ids[i], node_ids[i+1]);
        if (connecting_way == NULL) {
            printf("error: there are no roads between node %d and node %d.\n", node_ids[i], node_ids[i + 1]);
            return -1.0;
        }

        if (!are_nodes_adjacent_in_way(connecting_way, node_ids[i], node_ids[i+1])) {
            printf("error: cannot go directly from node %d to node %d.\n", node_ids[i], node_ids[i + 1]);
            return -1.0;
        }

        // If it's a one-way street, ensure the nodes are in the correct order
        if (connecting_way->oneway) {
            bool correct_order = false;
            for (int j = 0; j < connecting_way->num_nodes - 1; j++) {
                if (connecting_way->node_ids[j] == node_ids[i] && connecting_way->node_ids[j + 1] == node_ids[i + 1]) {
                    correct_order = true;
                    break;
                }
            }
            if (!correct_order) {
                printf("error: cannot go in reverse from node %d to node %d.\n", node_ids[i], node_ids[i + 1]);
                return -1.0;
            }
        }
    }
    //Error 5
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (node_ids[i] == node_ids[j]) {
                printf("error: node %d appeared more than once.\n", node_ids[i]);
                return -1.0;
            }
        }
    }
    //calculate time
    double total_travel_time = 0.0;
    for (int i = 0; i < size - 1; i++) {
        struct node *current_node = &m->nodes[node_ids[i]];
        struct node *next_node = &m->nodes[node_ids[i + 1]];
        const struct way *connecting_way = find_way_with_nodes(m,node_ids[i],node_ids[i+1]);
        
        // Calculate distance and travel time for this segment
        double distance = distance_between_nodes(current_node, next_node);
        double segment_time = distance / connecting_way->maxspeed;
        total_travel_time += segment_time;
    }

    return total_travel_time*60;
    
}

//--------------------------------Helper Functions--------------------------//
// Define the priority queue node
typedef struct PriorityQueueNode {
    int node_id;
    double priority; // This is the distance in Dijkstra's algorithm
    struct PriorityQueueNode *next;
} PriorityQueueNode;

// Define the priority queue (min-heap)
typedef struct PriorityQueue {
    PriorityQueueNode *head;
} PriorityQueue;

// Function to create a new Priority Queue
PriorityQueue *create_priority_queue() {
    PriorityQueue *pq = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    if (!pq) {
        fprintf(stderr, "error: out of memory.\n");
        return NULL;
    }
    pq->head = NULL;
    return pq;
}

// Function to insert a new node into the priority queue
void pq_insert(PriorityQueue *pq, int node_id, double priority) {
    PriorityQueueNode *new_node = (PriorityQueueNode *)malloc(sizeof(PriorityQueueNode));
    if (!new_node) {
        fprintf(stderr, "error: out of memory.\n");
        return;
    }
    new_node->node_id = node_id;
    new_node->priority = priority;
    new_node->next = NULL;

    if (pq->head == NULL || pq->head->priority > priority) {
        // Insert at the beginning
        new_node->next = pq->head;
        pq->head = new_node;
    } else {
        // Find the correct position and insert
        PriorityQueueNode *current = pq->head;
        while (current->next != NULL && current->next->priority <= priority) {
            current = current->next;
        }
        new_node->next = current->next;
        current->next = new_node;
    }
}

// Function to pop the node with the minimum priority from the priority queue
int pq_pop(PriorityQueue *pq) {
    if (pq->head == NULL) {
        return -1; // Priority Queue is empty
    }

    int node_id = pq->head->node_id;
    PriorityQueueNode *temp = pq->head;
    pq->head = pq->head->next;
    free(temp);

    return node_id;
}

// Function to check if the priority queue is empty
int pq_is_empty(PriorityQueue *pq) {
    return pq->head == NULL;
}

// Function to free the priority queue
void pq_free(PriorityQueue *pq) {
    while (!pq_is_empty(pq)) {
        pq_pop(pq);
    }
    free(pq);
}


void ssmap_path_create(const struct ssmap *m, int start_id, int end_id) {
    if (start_id == end_id) {
        printf("%d\n", start_id);
        return;
    }

    // Create and initialize the priority queue
    PriorityQueue *pq = create_priority_queue();
    if (!pq) return; // If failed to allocate memory for priority queue

    // Create arrays to store distances and previous node indices
    double *distances = (double *)malloc(sizeof(double) * m->num_nodes);
    int *previous = (int *)malloc(sizeof(int) * m->num_nodes);
    if (!distances || !previous) {
        fprintf(stderr, "error: out of memory.\n");
        pq_free(pq);
        free(distances);
        free(previous);
        return;
    }
    for (int i = 0; i < m->num_nodes; i++) {
        distances[i] = INFINITY;
        previous[i] = -1;
    }

    // Set the distance for the start node
    distances[start_id] = 0.0;
    pq_insert(pq, start_id, 0.0);

    // Dijkstra's algorithm
    while (!pq_is_empty(pq)) {
        int current_id = pq_pop(pq);
        struct node current_node = m->nodes[current_id];

        if (current_id == end_id) {
            break; // We've reached the destination node
        }

        for (int i = 0; i < current_node.num_ways; i++) {
            struct way way = m->ways[current_node.way_ids[i]];

            // For each node in the way
            for (int j = 0; j < way.num_nodes; j++) {
                if (way.node_ids[j] == current_id && j < way.num_nodes - 1){
                    // Get the neighbor's id
                    int neighbor_id = way.node_ids[j + 1];

                    // Calculate the distance to the neighbor
                    double alt = distances[current_id] + distance_between_nodes(&current_node, &m->nodes[neighbor_id]);

                    // Check if we found a shorter path to neighbor
                    if (alt < distances[neighbor_id]) {
                        distances[neighbor_id] = alt;
                        previous[neighbor_id] = current_id;
                        pq_insert(pq, neighbor_id, alt);
                    }
                }
            }
        }
    }

    // Reconstruct the shortest path
    int u = end_id;
    if (previous[u] == -1) {
        // No path was found
        printf("error: could not find a path from node %d to node %d.\n", start_id, end_id);
    } else {
        // We have found a path, let's reconstruct it backwards
        int *path = (int *)malloc(sizeof(int) * m->num_nodes);
        if (!path) {
            fprintf(stderr, "error: out of memory.\n");
            free(distances);
            free(previous);
            pq_free(pq);
            return;
        }

        int path_length = 0;
        while (u != -1) {
            path[path_length++] = u;
            u = previous[u];
        }

        // Print the path in the correct order
        for (int i = path_length - 1; i >= 0; i--) {
            printf("%d ", path[i]);
        }
        printf("\n");

        free(path);
    }

    // Free allocated resources
    free(distances);
    free(previous);
    pq_free(pq);
}
