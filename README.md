<meta name="robots" content="noindex">


Code: code is available under the folder code/src/. <br/>
Data: data is available under the folder code/data/. <br/>

Instructions on running the code: <br/>

To run the code, please run exp/Panel.java and then input the commands based on the tasks. Below are the details and several sample tasks are provided.

We included the parameters in a soft-coded way. The program exp/Panel.java provides a command line interface for a user to input the commands in the console. Specifically, when the program exp/Panel.java run, a user is required to input the information in the console including: <br/>
1) data_set_name,  <br/>
2) which_method_to_run,  <br/>
3) ratio (the theta value), and  <br/>
4) number of trials, etc.

Note that based on the user input, the values of global variables like file_paths and parameters will be modified by the program correspondingly. In other words, a user does not need to modify the code by hand.

------------------------------------------------------------------------------------------------------------
Sample Task-1. Do the clustering using method PIC with no optimization technique (PIC_noPrune) on data set citeseer_cociting.

Commands:<br/>
citeseer_cociting<br/>
Clustering PIC_noPrune 0.31 0.31 20<br/>
0<br/>
<br/>

The explanation for the commands: <br/>
Line 1: data_set_name <br/>
Line 2: Clustering which_method_to_run min_theta# max_theta# trial_No# <br/>
Line 3: 0

Line 1 specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".  <br/>
Line 2 specifies the task and related parameters. The task is to do the "Clustering" using the method "PIC with no optimization technique". The theta value varies from "0.31" to "0.31". We will do the clustering "20" times.  <br/>
Line 3 terminates the program.  <br/>
(NOTE: Please remember to type "0" at the end, which is for terminating the program, then click Enter)  

Output:<br/>
The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-2. Evaluate the clustering results of PIC on data set citeseer_cociting.

Commands:<br/>
citeseer_cociting<br/>
evaluation PIC 0.31 0.31 20<br/>
0<br/>
<br/>

The explanation for the commands:<br/>
Line 1: data_set_name<br/>
Line 2: evaluation the_results_of_which_method min_theta# max_theta# trial_No#<br/>
Line 3: 0

Line 1 specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".  <br/>
Line 2 specifies the task and related parameters. The task is to "evaluate" the clustering results discovered by the method "PIC". The theta value varies from "0.31" to "0.31". We will evaluate "20" sets of clustering results to calculate the average metric scores.  <br/>
Line 3 terminates the program.  <br/>
(NOTE: Please remember to type "0" at the end, which is for termination of the program, then click Enter)  

Output:<br/>
The metric scores obtained would be stored under the folder data/citeseer_cociting/measures/.

------------------------------------------------------------------------------------------------------------

Sample Task-3. Calculate the running time of the method PIC with no optimization technique (PIC_noPrune) on data set citeseer_cociting.

Commands:<br/>
citeseer_cociting<br/>
runningTime PIC_noPrune 0.31 20<br/>
0<br/>
<br/>

The explanation for the commands:<br/>
Line 1: data_set_name<br/>
Line 2: runningTime which_method_to_run theta# trial_No#<br/>
Line 3: 0

Line 1 specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".  <br/>
Line 2 specifies the task and related parameters. The task is to calculate the "running time" of method "PIC with no optimization technique". The theta value is "0.31". We will run the method by "20" times to calculate the average.  <br/>
Line 3 terminates the program.  <br/>

Output:<br/>
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-4. Do the clustering using method PIC with optimization techniques 1 and 2 (PIC_prune12) on data set citeseer_cociting.

Commands:<br/>
citeseer_cociting<br/>
Clustering PIC_prune12 0.31 0.31 20<br/>
0<br/>
<br/>

Output:<br/>
The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

------------------------------------------------------------------------------------------------------------

Sample Task-5. Calculate the running time of method PIC with optimization techniques 1 and 2 (PIC_prune12) on data set citeseer_cociting.

Commands:<br/>
citeseer_cociting<br/>
runningTime PIC_prune12 0.31 20<br/>
0<br/>
<br/>

Output:<br/>
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.
