package exp;

import java.util.Scanner;

import community.clustering.Clustering_PIC;
import evaluation.SingleGTruth;
import exp.scalability.Scalability;
import exp.sensitivity.ProcessDataSets;
import graph.Graph;
import hyperGraph.Hypergraph;
import statistic.Statistic_Hypergraph;
import utilities.FilePath_Mon;


public class Panel {
	
	public static void main(String arg[]) throws Exception {
		
		Scanner myObj = new Scanner(System.in); // Create a Scanner object
		
		while (true) {
			String input = myObj.nextLine(); // Read user input
			if (input.equals("0")) break;
			String[] graphs = input.split(",");
			
			input = myObj.nextLine(); // Read user input
			String[] strs = input.split(" ");
			
			for (String graph : graphs) {
				
				FilePath_Mon.filePathPre = System.getProperty("user.dir") + "/data/" + graph;
				Hypergraph.loadGraph();
				Graph.loaded = false;
				
				String PACKAGE = strs[0];
				String FUNCTION = strs[1];
				
				switch (PACKAGE) {
				case "Clustering":
					{
						switch (FUNCTION) {
						case "PIC":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double minratio = Double.parseDouble(strs[2]);
							double maxRatio = Double.parseDouble(strs[3]);
							int trials = Integer.parseInt(strs[4]);
							
							for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.doClustering_PIC(toHigherOrder, ordering, moveStrategy, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
							break;
							
						case "PIC_noPrune":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double minratio = Double.parseDouble(strs[2]);
							double maxRatio = Double.parseDouble(strs[3]);
							int trials = Integer.parseInt(strs[4]);
							
							for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.doClustering_PIC_noPrune(toHigherOrder, ordering, moveStrategy, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
							break;
							
						case "PIC_prune1":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double minratio = Double.parseDouble(strs[2]);
							double maxRatio = Double.parseDouble(strs[3]);
							int trials = Integer.parseInt(strs[4]);
							
							for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.doClustering_PIC_prune1(toHigherOrder, ordering, moveStrategy, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
							break;
							
						case "PIC_prune2":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double minratio = Double.parseDouble(strs[2]);
							double maxRatio = Double.parseDouble(strs[3]);
							int trials = Integer.parseInt(strs[4]);
							
							for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.doClustering_PIC_prune2(toHigherOrder, ordering, moveStrategy, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
							break;
							
						case "PIC_prune12":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double minratio = Double.parseDouble(strs[2]);
							double maxRatio = Double.parseDouble(strs[3]);
							int trials = Integer.parseInt(strs[4]);
							
							for (double ratio = minratio; ratio <= maxRatio; ratio += 0.1) {
								Clustering_PIC clustering = new Clustering_PIC();
								clustering.doClustering_PIC_prune12(toHigherOrder, ordering, moveStrategy, ratio, trials);
								Hypergraph.garbbageCollector.gc();
							}
						}
							break;
						}
					}
					
					break;
					
				case "runningTime":
					{
						switch (FUNCTION) {
						case "PIC":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double ratio = Double.parseDouble(strs[2]);
							int trials = Integer.parseInt(strs[3]);
							
							Clustering_PIC clustering = new Clustering_PIC();
							clustering.runningTime_PIC(toHigherOrder, ordering, moveStrategy, ratio, trials);
							Hypergraph.garbbageCollector.gc();
						}
							break;
							
						case "PIC_noPrune":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double ratio = Double.parseDouble(strs[2]);
							int trials = Integer.parseInt(strs[3]);
							
							Clustering_PIC clustering = new Clustering_PIC();
							clustering.runningTime_PIC_noPrune(toHigherOrder, ordering, moveStrategy, ratio, trials);
							Hypergraph.garbbageCollector.gc();
						}
							break;
							
						case "PIC_prune1":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double ratio = Double.parseDouble(strs[2]);
							int trials = Integer.parseInt(strs[3]);
							
							Clustering_PIC clustering = new Clustering_PIC();
							clustering.runningTime_PIC_prune1(toHigherOrder, ordering, moveStrategy, ratio, trials);
							Hypergraph.garbbageCollector.gc();
						}
							break;
							
						case "PIC_prune2":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double ratio = Double.parseDouble(strs[2]);
							int trials = Integer.parseInt(strs[3]);
							
							Clustering_PIC clustering = new Clustering_PIC();
							clustering.runningTime_PIC_prune2(toHigherOrder, ordering, moveStrategy, ratio, trials);
							Hypergraph.garbbageCollector.gc();
						}
							break;
							
						case "PIC_prune12":
						{
							boolean toHigherOrder = false;
							String ordering = "randomOrder";
							String moveStrategy = "move";
							double ratio = Double.parseDouble(strs[2]);
							int trials = Integer.parseInt(strs[3]);
							
							Clustering_PIC clustering = new Clustering_PIC();
							clustering.runningTime_PIC_prune12(toHigherOrder, ordering, moveStrategy, ratio, trials);
							Hypergraph.garbbageCollector.gc();
						}
							break;
							
						}
					}
				
				break;
				
				case "evaluation":
				{
					
					switch (FUNCTION) {
					case "PIC":
					{
						boolean toHigherOrder = false;
						String ordering = "randomOrder";
						String moveStrategy = "move";
						int maxClusterSize = 0;
						double minRatio = Double.parseDouble(strs[2]);
						double maxRatio = Double.parseDouble(strs[3]);
						int trials = Integer.parseInt(strs[4]);
						
						for (double ratio = minRatio; ratio <= maxRatio; ratio += 0.1) {
							SingleGTruth calculator = new SingleGTruth(false, moveStrategy, maxClusterSize, toHigherOrder, ordering, ratio);
							calculator.DiscoverLoaded = false;
							calculator.GraphLoaded = false;
							calculator.PRF_Prime(trials, maxClusterSize);
							calculator.NMI(trials, maxClusterSize);
							calculator.Purity(trials, maxClusterSize);
							calculator.ARI_Prime(trials, maxClusterSize);
							
							Hypergraph.garbbageCollector.gc();
						}
					}
					break;
					
					case "Baselines":
					{
						if (strs.length == 4) {
							String method = strs[2];
							int maxClusterSize = 0;
							int trials = Integer.parseInt(strs[3]);
							
							SingleGTruth calculator = new SingleGTruth(true, method, maxClusterSize);
							calculator.DiscoverLoaded = false;
							calculator.GraphLoaded = false;
							calculator.PRF_Prime(trials, maxClusterSize);
							calculator.NMI(trials, maxClusterSize);
							calculator.Purity(trials, maxClusterSize);
							calculator.ARI_Prime(trials, maxClusterSize);
						}
					}
					
					break;
					
					}
				}
				break;
				
				case "Scalability":
				{
					
					switch (FUNCTION) {
					case "synCardinalityDistribution":
					{
						double minLambda = Double.parseDouble(strs[2]);
						double maxLambda = Double.parseDouble(strs[3]);
						
						for (double lambda = minLambda; lambda <= maxLambda; lambda += 0.2) {
							lambda = (double) Math.round(lambda * 10d) / 10d;
							ProcessDataSets.synCardinalityDistribution(Hypergraph.getEdgeSize(), lambda);
						}
					}
					break;
					
					case "genRandomHypergraph_onLambda":
					{
						double minLambda = Double.parseDouble(strs[2]);
						double maxLambda = Double.parseDouble(strs[3]);
						
						for (double lambda = minLambda; lambda <= maxLambda; lambda += 0.2) {
							lambda = (double) Math.round(lambda * 10d) / 10d;
							ProcessDataSets.genRandomHypergraph_onLambda(lambda);
						}
					}
					break;
					
					case "genRandomHypergraph_onNodeSize":
					{
						int copies = Integer.parseInt(strs[2]);
						
						ProcessDataSets.genRandomHypergraph_onNodeSize(copies);
					}
					break;
					
					case "statistic_onLambda":
					{
						for (double lambda = 0.1; lambda <= 0.9; lambda += 0.2) {
							lambda = (double) Math.round(lambda * 10d) / 10d;
							Hypergraph.loadsynGraph(lambda, "lambda");
							
							Scalability.statistic_scalability();
						}
					}
					break;
					
					case "statistic_onNodeSize":
					{
						for (double divideRatio = 0.2; divideRatio <= 1.0; divideRatio += 0.2) {
							divideRatio = (double) Math.round(divideRatio * 10d) / 10d;
							Hypergraph.loadsynGraph(divideRatio, "divide");
							
							Scalability.statistic_scalability();
						}
					}
					break;
					
					case "onLambda":
					{
						
						boolean toHigherOrder = Boolean.parseBoolean(strs[2]);
						String ordering = strs[3];
						String moveStrategy = strs[4];
						double ratio = Double.parseDouble(strs[5]);
						int trials = Integer.parseInt(strs[6]);
						
						for (double lambda = 0.1; lambda > 0; lambda -= 0.2) {
//						for (double lambda = minLambda; lambda <= maxLambda; lambda += 0.2) {
							lambda = (double) Math.round(lambda * 10d) / 10d;
							Hypergraph.loadsynGraph(lambda, "lambda");
							
							Scalability.runningTime_onLambda(toHigherOrder, ordering, moveStrategy, ratio, trials, lambda);
						}
					}
					break;
					
					case "onNodeSize":
					{
						
						boolean toHigherOrder = Boolean.parseBoolean(strs[2]);
						String ordering = strs[3];
						String moveStrategy = strs[4];
						double ratio = Double.parseDouble(strs[5]);
						int trials = Integer.parseInt(strs[6]);
						
						for (double divideRatio = 0.2; divideRatio < 1.0; divideRatio += 0.2) {
							divideRatio = (double) Math.round(divideRatio * 10d) / 10d;
							Hypergraph.loadsynGraph(divideRatio, "divide");
							
							Scalability.runningTime_onNodeSize(toHigherOrder, ordering, moveStrategy, ratio, trials, divideRatio);
						}
					}
					break;
					
					}
				}
				break;
				
				case "Statistic_Hypergraph":
					{
						switch (FUNCTION) {
						case "saveCardinality":
							
							Statistic_Hypergraph.saveCardinality();
							
							break;
						}
					}
					
					break;
				}
			}
		}
		
		myObj.close();
	}
}
