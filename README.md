# Modularity_based_Hypergraph_Clustering_code

To run the code, please run exp/Panel.java then input the commands based on the tasks. For example,

------------------------------------------------------------------------------------------------------------
Sample Task-1: To do the clustering (using PIC without optimization technique) on data set citeseer_cociting, the commands will be

citeseer_cociting
Clustering PIC_noPrune 0.3 0.3 20
0

Explainantion for the commands: 
data_set_name 
Clustering method_to_run min_theta# max_theta# trial_No# 
0

(NOTE: Please remember to type "0" at then end, which is for termination of the program, then click Enter)  To save the time, the program will be run in parallel on multi-cores. The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-2: To do the evaluation on the results (of PIC) on data sets citeseer_cociting, the commands will be

citeseer_cociting
evaluation PIC 0.3 0.3 20
0

Explainantion for the commands:
data_set_name
evaluation method min_theta# max_theta# trial_No#
0

(NOTE: Please remember to type "0" at then end, which is for termination of the program, then click Enter)  To save the time, the program will be run in parallel on multi-cores. The metric scores obtained would be stored under the folder data/citeseer_cociting/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-3: To test the running time of PIC without optimization technique on data sets citeseer_cociting, the commands will be

citeseer_cociting
runningTime PIC_noPrune 0.3 20
0

Explainantion for the commands:
data_set_name
evaluation method_to_run theta# trial_No#
0

Different from above Task-1 and Task-2, the program of Task-3 to Task-6 will be run in one-core only. The running time obtained would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-4: To test the running time of PIC1 on data sets citeseer_cociting, the commands will be

citeseer_cociting
runningTime PIC_prune1 0.3 20
0

The program will be run in one-core only. The running time obtained would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-5: To test the running time of PIC2 on data sets citeseer_cociting, the commands will be

citeseer_cociting
runningTime PIC_prune2 0.3 20
0

The program will be run in one-core only. The running time obtained would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-6: To test the running time of PIC12 on data sets citeseer_cociting, the commands will be

citeseer_cociting
runningTime PIC_prune12 0.3 20
0

The program will be run in one-core only. The running time obtained would be saved under the folder data/citeseer_cociting/clustering/pic/.
------------------------------------------------------------------------------------------------------------
