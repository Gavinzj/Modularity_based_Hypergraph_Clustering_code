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
import java.util.List;
import java.util.Scanner;

import utilities.FilePath_Mon;

// for processing the data sets available in https://www.cs.cornell.edu/~arb/data/
public class RawProcess_dwan {
	private static String dataset = "DAWN";
	private static double sampleRatio = 0.1;
	private static int season = 0;
	private static String fileOutput_hypergraph = FilePath_Mon.filePathPre + "/dawn_s" + season + "_f" + sampleRatio + ".txt";
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_s" + season + 
			"_f" + sampleRatio + ".txt";
	
	private static boolean hasNodeLabel = true;
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "-full/";
	private static String fileInput_simplex_times = fileInputPre + dataset + "-times.txt";
	private static String fileInput_nverts = fileInputPre + dataset + "-nverts.txt";
	private static String fileInput_simplices = fileInputPre + dataset + "-simplices.txt";
	private static String fileInput_node_label = fileInputPre + dataset + "-node-labels.txt";
	private static String fileInput_node_labelHigherOrder = fileInputPre + dataset + "_Final.txt";
	
	static HashSet<Integer> nodes;
	static HashSet<String> nodeHigherOrderLabels;
	static HashSet<String> higherOrderLabels;
	static HashMap<Integer, String> nodes_HigherOrderLabels;
	static HashMap<String, String> label_HigherOrderLabel;
	
	static int[][] hyperedges_nodes;
	static int[] hyperedges_time;
	static HashMap<Integer, Integer> timeSeasons;
	
	public static void read() throws IOException {
		
		////////////////////////////////////////read nverts.txt
		Path path = Paths.get(fileInput_nverts);
		Scanner scanner = new Scanner(path.toAbsolutePath());
		
		int edgeSize = (int) Files.lines(path).count();
	
		nodes = new HashSet<Integer>();
		hyperedges_nodes = new int[edgeSize][];
		
		int edgeID = 0;
		while (scanner.hasNextLine()) {
			int cardinality = Integer.parseInt(scanner.nextLine());
			hyperedges_nodes[edgeID++] = new int[cardinality];
		}
		scanner.close();
		
		System.out.println("# of hyperedge " + edgeID);
		
		////////////////////////////////////////read times.txt
		path = Paths.get(fileInput_simplex_times);
		scanner = new Scanner(path.toAbsolutePath());
	
		hyperedges_time = new int[edgeSize];
		HashSet<Integer> distinctTimes = new HashSet<Integer>();
		
		edgeID = 0;
		while (scanner.hasNextLine()) {
			int time = Integer.parseInt(scanner.nextLine());
			hyperedges_time[edgeID++] = time;
			distinctTimes.add(time);
		}
		scanner.close();
		
		System.out.println("# of hyperedge " + edgeID + " # of distice time " + distinctTimes.size());
		
		List<Integer> times = new ArrayList<Integer>();
		Iterator<Integer> itr_time = distinctTimes.iterator();
		while (itr_time.hasNext()) times.add(itr_time.next());
		
		Collections.sort(times);
		
		timeSeasons = new HashMap<Integer, Integer>();
		int timeID = 0;
		for (int i = 0; i < times.size(); i++) timeSeasons.put(times.get(i), timeID++);
		
		System.out.println(timeSeasons);
		
		//////////////////////////////////////// read simplices.txt
		path = Paths.get(fileInput_simplices);
		scanner = new Scanner(path.toAbsolutePath());
		
		List<List<Integer>> hyperedges = new ArrayList<List<Integer>>(hyperedges_nodes.length);
		for (int i = 0; i < hyperedges_nodes.length; i++) {
			int cardinality = hyperedges_nodes[i].length;
			int nodeID;
			
			List<Integer> hyperedge = new ArrayList<Integer>(cardinality);
			for (int j = 0; j < cardinality; j++) {
				nodeID = Integer.parseInt(scanner.nextLine());
				hyperedge.add(nodeID);
			}
			hyperedges.add(hyperedge);
		}
		scanner.close();
		
		edgeID = 0;
		int totalCardinality = 0;
		for (int i = 0; i < hyperedges.size(); i++) {
			List<Integer> hyperedge = hyperedges.get(i);
			int cardinality = hyperedge.size();
			int time = hyperedges_time[i];
			
			if (cardinality < 2) {
				hyperedges_nodes[i] = new int[0];
				continue;
			}
			if (Math.random() > sampleRatio) {
				hyperedges_nodes[i] = new int[0];
				continue;
			}
			if (timeSeasons.get(time) != season) {
				hyperedges_nodes[i] = new int[0];
				continue;
			}
			
//			System.out.println(cardinality + " " + time);
			
			int nodeID;
			hyperedges_nodes[i] = new int[cardinality];
			for (int j = 0; j < cardinality; j++) {
				nodeID = hyperedge.get(j);
				nodes.add(nodeID);
				hyperedges_nodes[i][j] = nodeID;
				totalCardinality++;
			}
			
			edgeID++;
		}
		
		System.out.println("# of hyperedge " + edgeID + " total cardinality: " + totalCardinality + " node size " + nodes.size());
		
		////////////////////////////////////////read DAWN_Final.xlsx
		if (hasNodeLabel) {
			path = Paths.get(fileInput_node_labelHigherOrder);
			scanner = new Scanner(path.toAbsolutePath());
			
			label_HigherOrderLabel = new HashMap<String, String>();
			higherOrderLabels = new HashSet<String>();
			
			String line;
			String[] strs;
			String label;
			String higherOrderLabel;
			while (scanner.hasNextLine()) {
				line = scanner.nextLine().toLowerCase();
				strs = line.split("\t");
				
//				System.out.println(line);
				
				label = strs[0];
				higherOrderLabel = strs[2];

//				System.out.println("label " + label + " higherOrderLabel " + higherOrderLabel);
				
				higherOrderLabels.add(label);
				label_HigherOrderLabel.put(label, higherOrderLabel);
			}
			scanner.close();
			
			System.out.println("# of higher-order label: " + higherOrderLabels.size());
		}
		
		
		////////////////////////////////////////read node-labels.txt
		if (hasNodeLabel) {
			path = Paths.get(fileInput_node_label);
			scanner = new Scanner(path.toAbsolutePath());
			
			nodes_HigherOrderLabels = new HashMap<Integer, String>();
			nodeHigherOrderLabels = new HashSet<String>();
			
			String line;
			String[] strs;
			String label;
			String higherOrderLabel;
			int nodeID = 0;
			while (scanner.hasNextLine()) {
				line = scanner.nextLine().toLowerCase();
				strs = line.split(" ");
				
				nodeID = Integer.parseInt(strs[0]);
				label = strs[1];
				
				if (!label_HigherOrderLabel.containsKey(label)) {
					higherOrderLabel = "others";
//					System.out.println("no nodeID " + nodeID + " label " + label);
				} else {
					higherOrderLabel = label_HigherOrderLabel.get(label);
				}
				
//				System.out.println("nodeID " + nodeID + " label " + label + " higher-order-label " + higherOrderLabel);
				
				nodeHigherOrderLabels.add(higherOrderLabel);
				nodes_HigherOrderLabels.put(nodeID, higherOrderLabel);
			}
			scanner.close();
			
			System.out.println("# of labeled node: " + nodes_HigherOrderLabels.size());
			System.out.println("# of node label: " + nodeHigherOrderLabels.size());
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
				
				if (existHyperedges.contains(output)) continue;
				existHyperedges.add(output);
				
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of edges " + existHyperedges.size());
		
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
					int nodeID = nodeIDs.get(i);
					fwCount.write(nodeID + "\t" + nodes_HigherOrderLabels.get(nodeID) + "\n");
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
