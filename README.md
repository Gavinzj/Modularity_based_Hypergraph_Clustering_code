<meta name="robots" content="noindex">
# Modularity_based_Hypergraph_Clustering_code

To run the code, please run exp/Panel.java then input the commands based on the tasks. Below lists several sample tasks.

------------------------------------------------------------------------------------------------------------
Sample Task-1. Do the clustering (using PIC without optimization technique) on data set citeseer_cociting.

Commands:

citeseer_cociting
Clustering PIC_noPrune 0.3 0.3 20
0

Explainantion for the commands: 
data_set_name 
Clustering method_to_run min_theta# max_theta# trial_No# 
0

(NOTE: Please remember to type "0" at then end, which is for termination of the program, then click Enter)  

Output:
The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-2. Evaluate the clustering results (of PIC) on data sets citeseer_cociting.

Commands:

citeseer_cociting
evaluation PIC 0.3 0.3 20
0

Explainantion for the commands:
data_set_name
evaluation method min_theta# max_theta# trial_No#
0

(NOTE: Please remember to type "0" at then end, which is for termination of the program, then click Enter)  

Output:
The metric scores obtained would be stored under the folder data/citeseer_cociting/measures/.

------------------------------------------------------------------------------------------------------------

Sample Task-3. Test the running time of PIC without optimization technique on data sets citeseer_cociting.

Commands:

citeseer_cociting
runningTime PIC_noPrune 0.3 20
0

Explainantion for the commands:
data_set_name
evaluation method_to_run theta# trial_No#
0

Output:
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-4. Test the running time of PIC1 on data sets citeseer_cociting.

Commands:

citeseer_cociting
runningTime PIC_prune1 0.3 20
0

Output:
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-5. Test the running time of PIC2 on data sets citeseer_cociting.

Commands:

citeseer_cociting
runningTime PIC_prune2 0.3 20
0

Output:
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-6. Test the running time of PIC12 on data sets citeseer_cociting.

Commands:

citeseer_cociting
runningTime PIC_prune12 0.3 20
0

Output:
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.
------------------------------------------------------------------------------------------------------------
