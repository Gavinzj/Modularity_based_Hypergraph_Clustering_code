package preprocess.others;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;

import utilities.FilePath_Mon;

// for processing the data sets available in https://www.cs.cornell.edu/~arb/data/
public class RawProcess_ARB {
	private static String dataset = "amazon_reviews-full-full";
	private static String fileOutput_hypergraph = FilePath_Mon.filePathPre + "/amazon_reviews.txt";
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	private static String fileOutput_hyperedge_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/hyperedge_groundTruth.txt";
	
	private static boolean hasNodeLabel = true;
	private static boolean hasHyperedgeLabel = false;
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	private static String fileInput_nverts = fileInputPre + dataset + "-nverts.txt";
	private static String fileInput_simplices = fileInputPre + dataset + "-simplices.txt";
	private static String fileInput_times = fileInputPre + dataset + "-times.txt";
	private static String fileInput_node_label = fileInputPre + dataset + "-node-labels.txt";
	private static String fileInput_simplex_label = fileInputPre + dataset + "-simplex-labels.txt";
	
	
	static HashSet<Integer> nodes;
	static HashSet<String> nodeLabels;
	static HashMap<Integer, String> nodes_Labels;
	
	static int[][] hyperedges_nodes;
	static double[] hyperedge_times;
	static int[] hyperedges_Labels;
	static HashMap<String, Integer> hyperedgeLabels_IDs;
	
	public static void read() throws IOException {
		
		////////////////////////////////////////read nverts.txt
		Path path = Paths.get(fileInput_nverts);
		Scanner scanner = new Scanner(path.toAbsolutePath());
		
		int edgeSize = (int) Files.lines(path).count();
	
		nodes = new HashSet<Integer>();
		hyperedges_nodes = new int[edgeSize][];
		hyperedge_times = new double[edgeSize];
		
		int edgeID = 0;
		int totalCardinality = 0;
		while (scanner.hasNextLine()) {
			// process each line
			int cardinality = Integer.parseInt(scanner.nextLine());
			
			hyperedges_nodes[edgeID++] = new int[cardinality];
			
			totalCardinality += cardinality;
		}
		
		scanner.close();
		
		System.out.println("# of hyperedge " + hyperedges_nodes.length + " edgeID " + edgeID + " total cardinality: " + totalCardinality);
		
		//////////////////////////////////////// read simplices.txt
		path = Paths.get(fileInput_simplices);
		scanner = new Scanner(path.toAbsolutePath());
		
		totalCardinality = 0;
		for (int i = 0; i < hyperedges_nodes.length; i++) {
			int cardinality = hyperedges_nodes[i].length;
			int nodeID;
			for (int j = 0; j < cardinality; j++) {
				nodeID = Integer.parseInt(scanner.nextLine());
				
				nodes.add(nodeID);
				
				hyperedges_nodes[i][j] = nodeID;
				
				totalCardinality++;
			}
		}
		
		scanner.close();
		
		System.out.println("total cardinality: " + totalCardinality + " node size " + nodes.size());
		
		////////////////////////////////////////read times.txt
		path = Paths.get(fileInput_times);
		scanner = new Scanner(path.toAbsolutePath());
		
		edgeID = 0;
		while (scanner.hasNextLine()) {
			// process each line
			double time = Double.parseDouble(scanner.nextLine());
			
			hyperedge_times[edgeID++] = time;
		}
		
		scanner.close();
		
		System.out.println("edgeID: " + edgeID);
		
		
		////////////////////////////////////////read node-labels.txt
		if (hasNodeLabel) {
			path = Paths.get(fileInput_node_label);
			scanner = new Scanner(path.toAbsolutePath());
			
			nodes_Labels = new HashMap<Integer, String>();
			nodeLabels = new HashSet<String>();
			
			String[] strs;
			String label;
			int nodeID = 0;
			while (scanner.hasNextLine()) {
				strs = scanner.nextLine().split(" ");
				
				label = strs[strs.length - 1];
				nodeID = Integer.parseInt(strs[0]);
				
				nodeLabels.add(label);
				nodes_Labels.put(nodeID, label);
			}
			
			scanner.close();
			
			System.out.println("# of node label: " + nodeLabels.size());
		}
		
		////////////////////////////////////////read simplex-label.txt
		if (hasHyperedgeLabel) {
			path = Paths.get(fileInput_simplex_label);
			scanner = new Scanner(path.toAbsolutePath());
			
			hyperedges_Labels = new int[hyperedges_nodes.length];
			hyperedgeLabels_IDs = new HashMap<String, Integer>();
			
			edgeID = 0;
			int labelID = 0;
			while (scanner.hasNextLine()) {
				String[]  strs = scanner.nextLine().split(" ");
				String label = strs[0];
				
				if (!hyperedgeLabels_IDs.containsKey(label)) {
					hyperedgeLabels_IDs.put(label, labelID++);
				}
				
				hyperedges_Labels[edgeID++] = hyperedgeLabels_IDs.get(label);
			}
			
			scanner.close();
			
			System.out.println("# of hyperedge label: " + hyperedgeLabels_IDs.size() + " hyperedge labelID " + labelID);
		}
	}
	
	public static void save() {
		// save hypegraph
		int minHyperedgeSize = 2;
		HashSet<String> existHyperedges = new HashSet<String>(hyperedges_nodes.length);
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph);
			
			for (int edgeID = 0; edgeID < hyperedges_nodes.length; edgeID++) {
				
				if (hyperedges_nodes[edgeID].length < minHyperedgeSize) {
					continue;
				}
				
				Arrays.sort(hyperedges_nodes[edgeID]);
				
				String output = hyperedges_nodes[edgeID][0] + "";
				for (int i = 1; i < hyperedges_nodes[edgeID].length; i ++) {
					output = output + "\t" + hyperedges_nodes[edgeID][i];
				}
				
				if (existHyperedges.contains(output)) {
					continue;
				}
				existHyperedges.add(output);
				
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		
		// save node class ground truth
		if (hasNodeLabel) {
			try {
				FileWriter fwCount = new FileWriter(fileOutput_node_groundTruth);
				
				ArrayList<Integer> nodeIDs = new ArrayList<Integer>();
				
				Iterator<Integer> it = nodes.iterator();
		        while (it.hasNext()) {
		        	nodeIDs.add(it.next());
		        }
		        
		        Collections.sort(nodeIDs);
		        
				for (int i = 0; i < nodeIDs.size(); i++) {
					fwCount.write(nodeIDs.get(i) + "\t" + nodes_Labels.get(nodeIDs.get(i)) + "\n");
				}
				
				fwCount.close();
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			}
		}
		
		
		
		// save hyperedge class ground truth
		if (hasHyperedgeLabel) {
			try {
				FileWriter fwCount = new FileWriter(fileOutput_hyperedge_groundTruth);
		        
				for (int i = 0; i < hyperedges_Labels.length; i++) {
					fwCount.write(i + "\t" + hyperedges_Labels[i] + "\n");
				}
				
				fwCount.close();
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			}
		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
		save();
	}
}
