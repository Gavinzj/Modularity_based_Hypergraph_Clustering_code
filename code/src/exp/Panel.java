package exp;

import java.util.Scanner;

import community.clustering.Clustering_PIC;
import evaluation.SingleGTruth;
import graph.Graph;
import hyperGraph.Hypergraph;
import utilities.FilePath_Mon;

// To process the user input, modify the global variables correspondingly, then run the program
public class Panel {
	
	public static void main(String arg[]) throws Exception {
		
		// Create a Scanner object to read the user input
		Scanner myObj = new Scanner(System.in); 
		
		// read and process the user input
		while (true) {
			
			// Read user input
			String input = myObj.nextLine(); 
			
			// if user input is "0", terminate the program
			if (input.equals("0")) break;
			String[] graphs = input.split(",");
			
			// process user input
			input = myObj.nextLine();
			String[] strs = input.split(" ");
			
			// for each data set that user input
			for (String graph : graphs) {
				
				// Modify the values of global variable filePathPre based on the user input
				FilePath_Mon.filePathPre = System.getProperty("user.dir") + "/data/" + graph;
				
				// load the graph
				Hypergraph.loadGraph();
				Graph.loaded = false;
				
				// tasks include: Clustering, runningTime, evaluation, etc.
				String PACKAGE = strs[0];
				// methods include: PIC, PIC's variants with optimization techniques, and baselines
				String FUNCTION = strs[1];
				
				switch (PACKAGE) {
				// to do clustering
				case "Clustering":
					{
						switch (FUNCTION) {
							
							case "PIC_noPrune":
							{
								// hypergraph clustering using PIC with no optimization technique
								// varying on ratio -- theta
								
								String ordering = "randomOrder";	// by default (we fix it here), the order of node move is random
								double minratio = Double.parseDouble(strs[2]);	// minimum ratio (theta)
								double maxRatio = Double.parseDouble(strs[3]);	// maximum ratio (theta)
								int trials = Integer.parseInt(strs[4]);	// number of trials that the algorithm runs
								
								// for each ratio (theta) vary from minratio to maxratio with step size 0.1
								for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
									Clustering_PIC clustering = new Clustering_PIC();
									clustering.doClustering_PIC_noPrune(ordering, ratio, trials);
									Hypergraph.garbbageCollector.gc();
								}
							}
								break;
							
							case "PIC_prune12":
							{
								// hypergraph clustering using PIC with optimization techniques 1 and 2
								// varying on ratio -- theta
								
								String ordering = "randomOrder";	// by default (we fix it here), the order of node move is random
								double minratio = Double.parseDouble(strs[2]);	// minimum ratio (theta)
								double maxRatio = Double.parseDouble(strs[3]);	// maximum ratio (theta)
								int trials = Integer.parseInt(strs[4]);	// number of trials that the algorithm runs
								
								// for each ratio (theta) vary from minratio to maxratio with step size 0.1
								for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
									Clustering_PIC clustering = new Clustering_PIC();
									clustering.doClustering_PIC_prune12(ordering, ratio, trials);
									Hypergraph.garbbageCollector.gc();
								}
							}
								break;
								
							default:
							{
								System.out.println("Invalid Function");
							}
						}
					}
					
					break;
					
				// to calculate the running time (clustering results will not be saved)
				case "runningTime":
					{
						switch (FUNCTION) {
							
							case "PIC_noPrune":
							{
								// hypergraph clustering using PIC with no optimization technique
								
								String ordering = "randomOrder";	// by default (we fix it here), the order of node move is random
								double ratio = Double.parseDouble(strs[2]);	// the value of ratio (theta)
								int trials = Integer.parseInt(strs[3]);	// calculate the average running time over several number of trials
								
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.runningTime_PIC_noPrune(ordering, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
								break;
								
							case "PIC_prune12":
							{
								// hypergraph clustering using PIC with optimization techniques 1 and 2
								
								String ordering = "randomOrder";	// by default (we fix it here), the order of node move is random
								double ratio = Double.parseDouble(strs[2]);	// the value of ratio (theta)
								int trials = Integer.parseInt(strs[3]);	// calculate the average running time over several number of trials
								
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.runningTime_PIC_prune12(ordering, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
								break;
								
							default:
							{
								System.out.println("Invalid Function");
							}
						}
					}
				
				break;
				
				// to evaluate the clustering results
				case "evaluation":
				{
					
					switch (FUNCTION) {
						case "PIC":
						{
							// evaluate the hypergraph clustering results using PIC algorithm
							// varying on ratio -- theta
							
							String ordering = "randomOrder";	// by default (we fix it here), the order of node move is random
							double minRatio = Double.parseDouble(strs[2]);	// minimum ratio (theta)
							double maxRatio = Double.parseDouble(strs[3]);	// maximum ratio (theta)
							int trials = Integer.parseInt(strs[4]);	// calculate the average metric scores over trials number of clustering results
							
							// for each ratio (theta) vary from minratio to maxratio with step size 0.1
							for (double ratio = minRatio; ratio <= maxRatio; ratio += 0.1) {
								SingleGTruth calculator = new SingleGTruth(false, 0, ordering, ratio);
								calculator.DiscoverLoaded = false;
								calculator.ComponentLoaded = false;
								calculator.PRF_Prime(trials);
								calculator.NMI(trials);
								calculator.ARI_Prime(trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
						break;
						
						default:
						{
							System.out.println("Invalid Function");
						}
					}
				}
				break;
				
				// default
				default:
					{
						System.out.println("Invalid Package");
					}
				}
			}
		}
		
		myObj.close();
	}
}
