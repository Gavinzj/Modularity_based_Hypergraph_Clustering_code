<meta name="robots" content="noindex">

<h2>Repository Overview</h2>
Code: code is available under the folder code/src/. <br/>
Data: data is available under the folder code/data/. 
Note: If a dataset contains multiple disconnected components, only the largest component is retained. <br/>

<h2>Instructions on running the code:</h2>

To run the code, please run exp/Panel.java and then input the commands based on the tasks. Below are the details and several sample tasks are provided.

We included the parameters in a soft-coded way. The program exp/Panel.java provides a command line interface for a user to input the commands in the console. Specifically, when the program exp/Panel.java run, a user is required to input the information in the console including: <br/>
<ol>
  <li>data_set_name</li>
  <li>which_method_to_run</li>
  <li>ratio (the theta value)</li>
  <li>number of trials, and</li>
  <li>other relevant parameters</li>
</ol>

Note that based on the user input, the values of global variables like file_paths and parameters will be modified by the program correspondingly. In other words, a user does not need to modify the code by hand.

<hr>

<h3>Sample Task 1</h3>
<p>Task: Do the clustering using method PIC with no optimization technique (PIC_noPrune) on data set citeseer_cociting.</p>


<p>Commands:</p>
<pre>
citeseer_cociting
Clustering PIC_noPrune 0.31 0.31 20
0
</pre>

<!-- Commands:<br/>
citeseer_cociting<br/>
Clustering PIC_noPrune 0.31 0.31 20<br/>
0<br/>
<br/> -->

Explanation of commands: <br/>
<ul>
  <li>Line 1:<code>data_set_name</code> → Specifies the dataset "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".</li>
  <li>Line 2:<code>Clustering which_method_to_run min_theta# max_theta# trial_No#</code> → Specifies the task and related parameters. The task is to do the "Clustering" using the method "PIC with no optimization technique". The theta value varies from "0.31" to "0.31". We will do the clustering "20" times.</li>
  <li>Line 3:<code>0</code> → Terminates the program.</li>
</ul>
<p>(NOTE: Please remember to type "0" at the end, which is for terminating the program, then click Enter)</p>

<!-- Line 1: data_set_name <br/>
Line 2: Clustering which_method_to_run min_theta# max_theta# trial_No# <br/>
Line 3: 0

Line 1 specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".  <br/>
Line 2 specifies the task and related parameters. The task is to do the "Clustering" using the method "PIC with no optimization technique". The theta value varies from "0.31" to "0.31". We will do the clustering "20" times.  <br/>
Line 3 terminates the program.  <br/>
(NOTE: Please remember to type "0" at the end, which is for terminating the program, then click Enter)   -->

Output:<br/>
The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

<hr>

<h3>Sample Task 2</h3>
<p>Task: Evaluate the clustering results of PIC on data set citeseer_cociting.</p>


<p>Commands:</p>
<pre>
citeseer_cociting
evaluation PIC 0.31 0.31 20
0
</pre>

<!-- Commands:<br/>
citeseer_cociting<br/>
evaluation PIC 0.31 0.31 20<br/>
0<br/>
<br/> -->

Explanation of commands: <br/>
<ul>
  <li>Line 1:<code>data_set_name</code> → Specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".</li>
  <li>Line 2:<code>evaluation the_results_of_which_method min_theta# max_theta# trial_No#</code> → Specifies the task and related parameters. The task is to "evaluate" the clustering results discovered by the method "PIC". The theta value varies from "0.31" to "0.31". We will evaluate "20" sets of clustering results to calculate the average metric scores.</li>
  <li>Line 3:<code>0</code> → Terminates the program.</li>
</ul>
<p>(NOTE: Please remember to type "0" at the end, which is for termination of the program, then click Enter)</p>


<!-- The explanation for the commands:<br/>
Line 1: data_set_name<br/>
Line 2: evaluation the_results_of_which_method min_theta# max_theta# trial_No#<br/>
Line 3: 0 -->

<!-- Line 1 specifies the data set "citeseer_cociting". The file paths will then point to the folder corresponding to the data set "citeseer_cociting".  <br/>
Line 2 specifies the task and related parameters. The task is to "evaluate" the clustering results discovered by the method "PIC". The theta value varies from "0.31" to "0.31". We will evaluate "20" sets of clustering results to calculate the average metric scores.  <br/>
Line 3 terminates the program.  <br/>
(NOTE: Please remember to type "0" at the end, which is for termination of the program, then click Enter)   -->

Output:<br/>
The metric scores obtained would be stored under the folder data/citeseer_cociting/measures/.

<hr>

<h3>Sample Task 3</h3>
<p>Task: Calculate the running time of the method PIC with no optimization technique (PIC_noPrune) on data set citeseer_cociting.</p>

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

<hr>

<h3>Sample Task 4</h3>
<p>Task: Do the clustering using method PIC with optimization techniques 1 and 2 (PIC_prune12) on data set citeseer_cociting.</p>

Commands:<br/>
citeseer_cociting<br/>
Clustering PIC_prune12 0.31 0.31 20<br/>
0<br/>
<br/>

Output:<br/>
The clustering results obtained would be stored under the folder data/citeseer_cociting/clustering/pic/.

<hr>

<h3>Sample Task 5</h3>
<p>Task: Calculate the running time of method PIC with optimization techniques 1 and 2 (PIC_prune12) on data set citeseer_cociting.</p>

Commands:<br/>
citeseer_cociting<br/>
runningTime PIC_prune12 0.31 20<br/>
0<br/>
<br/>

Output:<br/>
The running time would be saved under the folder data/citeseer_cociting/clustering/pic/.

<hr>

<h2>Instructions on adding Jama library (IntelliJ IDEA as an example):<h2>

<p>Download Jama:</p>
You can download the Jama JAR file from the official site (https://math.nist.gov/javanumerics/jama/).

    Create/ Open Your Project:
        Open IntelliJ IDEA and either create a new project or open your existing project.

    Add Jama JAR to Your Project:
        Right-click on your project folder in the Project pane.
        Select "New" > "Directory" and name it something like libs (or another name of your choice).
        Drag the Jama-1.0.3.jar file into the libs directory.

    Include Jama in Your Project:
        Right-click on Jama-1.0.3.jar in the libs directory.
        Select "Add as Library".
        In the dialog, choose "Module SDK" (if not automatically selected) and click OK.

    Verify Library Addition:
        In Project Structure, navigate to File > Project Structure > Libraries.
        Make sure Jama-1.0.3.jar is listed as one of the libraries.

    Use Jama in Your Code:
        Now you can use the Jama library in your Java classes