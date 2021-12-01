This program runs various flow analysis functions on individual flows.
To run this program execute the following command in a terminal on a machine with CUDA
enabled:

    ./proj3 config.txt


Configuration File Options

- load
    Loads a *.csv file with flow data into the program(each flow file is assigned an index starting with 0)
    Syntax:
    load `filename` `source_x_column` `source_y_column` `destination_x_column` `destination_y_column` `flow_length_column`

- filter
    Filters a flow file accourding to a set of lower and upper bounds in an attribute column of the *.csv
    Assumes that column has a floating point type
    (NOT APPLIED UNTIL THE `apply` command is used,filters can be concatenated)
    Syntax:
    filter `flow_file_index` `column` `lower_bound` `upper_bound`

- apply
    Applies filters previously set filters to a flow file
    Syntax:
    apply  `flow_file_index`

- knn
    Runs k-nearest neighbors function on a flow file for a particular index
    writes nearest neigbors to an output file
    Syntax:
    knn `flow_file_index` `number_neighbors` `flow_index` `func_type*` `alpha` `output_file`

- flow_k
    Finds clusters in flow file using flow k function.
    Radius is a spatial radius for calculating the k-function 
    Syntax:
    flow_k `flow_file_index` `number_iterations` `func_type` `alpha` `radius` `output_file`

- cross_flow_k
    Find joint clusters for two flows using the cross k function
    Syntax:
    cross_flow_k `flow_file_index_1` `flow_file_index_2` `number_interations` `func_type` `alpha` `radius` `output_file`

- colocation
    Finds spatial colocation patterns in multivariate flows
    Runs on all loaded flows (need at least two loaded in program)
    Syntax:
    extract_colocation_patterns `frequency_threshold` `spatial_threshold` `output_file`


*func_type: Determines which distance function is used to calculate distance between flows
    0: Flow distance
    1: Flow disimilarity
    2: Maximum distance

*alpha: Determines weight between source and destination ends of flow in calculatin flow distance
    Ranges [0,2]
    2 - Only use source end points
    1 - Balanced between source and destination
    0 - Only use destination endpoints


Example Configuration File
---------------------------------

load tmp1.csv pickup_latitude pickup_longitude dropoff_latitude dropoff_longitude trip_distance
load tmp2.csv pickup_latitude pickup_longitude dropoff_latitude dropoff_longitude trip_distance
load tmp3.csv pickup_latitude pickup_longitude dropoff_latitude dropoff_longitude trip_distance
filter 0 passenger_count 0 2
apply 0
filter 1 passenger_count 1 3
apply 1
filter 2 passenger_count 2 4
apply 2
knn 0 5 1 0 1.0 knn_out.csv
flow_k 0 1000 0 1 .1 flow_k_out.csv
cross_flow_k 0 1 1000 0 1.0 .1 cross_flow_k.csv
extract_colocation_patterns .1 .01 colocation_out.csv