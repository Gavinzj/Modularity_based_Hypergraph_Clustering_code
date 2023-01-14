package preprocess.others;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import utilities.FilePath_Mon;

// for processing the data sets available in https://www.cs.cornell.edu/~arb/data/
public class RawProcess_Terrorist {
	private static String dataset = "terrorist_attacks-full";
	private static String fileOutput_hypergraph_coLoc = FilePath_Mon.filePathPre + "/terrorist_coLoc.txt";
	private static String fileOutput_hypergraph_coLocOrg = FilePath_Mon.filePathPre + "/terrorist_coLocOrg.txt";
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	private static String fileInput_nodes = fileInputPre + "terrorist_attack_nodes.txt";
	private static String fileInput_coLoc = fileInputPre + "terrorist_attack_loc_edges.txt";
	private static String fileInput_coLocOrg = fileInputPre + "terrorist_attack_loc_org_edges.txt";
	
	static HashMap<String, Integer> nodes_IDs;
	static HashSet<String> labels;
	static HashMap<Integer, String> nodes_Labels;
	
	static HashMap<Integer, HashSet<Integer>> coLoc_attacks;
	static HashMap<Integer, HashSet<Integer>> coLocOrg_attacks;
	
	public static void read() throws IOException {
		
		int nodeCnt = 0;
		nodes_IDs = new HashMap<String, Integer>();
		labels = new HashSet<String>();
		nodes_Labels = new HashMap<Integer, String>();
		
		try (FileReader reader = new FileReader(fileInput_nodes);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, attact, label;
			int node;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
//				if (stop-- <= 0) continue;
				
//				System.out.println(line);
				attact = strs[0];
				label = strs[strs.length - 1];
				
				if (!nodes_IDs.containsKey(attact)) nodes_IDs.put(attact, nodeCnt++);
				node = nodes_IDs.get(attact);
				
				labels.add(label);
				
				if (!nodes_Labels.containsKey(node)) nodes_Labels.put(node, label);
			}
		}
		
		System.out.println("# of nodes " + nodes_IDs.size() + " " + nodeCnt + " " + nodes_Labels.size());
		System.out.println("# of labels " + labels.size());
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		
		coLoc_attacks = new HashMap<Integer, HashSet<Integer>>();
		
		try (FileReader reader = new FileReader(fileInput_coLoc);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, attack1, attack2;
			int node1, node2;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
//				if (stop-- <= 0) continue;
				
//				System.out.println(line);
				attack1 = strs[0];
				attack2 = strs[1];
				
				if (nodes_IDs.containsKey(attack1)) node1 = nodes_IDs.get(attack1);
				else continue;
				
				if (nodes_IDs.containsKey(attack2)) node2 = nodes_IDs.get(attack2);
				else continue;
				
				if (!coLoc_attacks.containsKey(node1)) coLoc_attacks.put(node1, new HashSet<Integer>());
				coLoc_attacks.get(node1).add(node2);
				
				if (!coLoc_attacks.containsKey(node2)) coLoc_attacks.put(node2, new HashSet<Integer>());
				coLoc_attacks.get(node2).add(node1);
				
//				System.out.println(node1 + " " + node2);
			}
		}
		
		System.out.println("# of coloc " + coLoc_attacks.size());
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
				
		
		coLocOrg_attacks = new HashMap<Integer, HashSet<Integer>>();
		
		try (FileReader reader = new FileReader(fileInput_coLocOrg);
		BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, attack1, attack2;
			int node1, node2;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				//if (stop-- <= 0) continue;
				
				//System.out.println(line);
				attack1 = strs[0];
				attack2 = strs[1];
				
				if (nodes_IDs.containsKey(attack1)) node1 = nodes_IDs.get(attack1);
				else continue;
				
				if (nodes_IDs.containsKey(attack2)) node2 = nodes_IDs.get(attack2);
				else continue;
				
				if (!coLocOrg_attacks.containsKey(node1)) coLocOrg_attacks.put(node1, new HashSet<Integer>());
				coLocOrg_attacks.get(node1).add(node2);
				
				if (!coLocOrg_attacks.containsKey(node2)) coLocOrg_attacks.put(node2, new HashSet<Integer>());
				coLocOrg_attacks.get(node2).add(node1);
			}
		}
		
		System.out.println("# of colocorg " + coLocOrg_attacks.size());
	}
	
	public static void save() {
		
		// save hyperedges for coLocOrg
		int minHyperedgeSize = 2;
		HashSet<String> visHyperedges = new HashSet<String>(coLoc_attacks.size());
		List<Integer> nodes;
		
		int counter = 0;
		Iterator<Integer> itr_coloc;
		Iterator<Integer> itr_attacks;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_coLoc);
			
			int coloc;
			String output;
			itr_coloc = coLoc_attacks.keySet().iterator();
			
			while (itr_coloc.hasNext()) {
				coloc = itr_coloc.next();
				
				nodes = new ArrayList<Integer>();
				nodes.add(coloc);
				
				itr_attacks = coLoc_attacks.get(coloc).iterator();
				while(itr_attacks.hasNext()) {
					nodes.add(itr_attacks.next());
				}
				
				if (nodes.size() < minHyperedgeSize) continue;
				
				Collections.sort(nodes);
				
				output = nodes.get(0) + "";
				for (int i = 1; i < nodes.size(); i ++) {
					output = output + "\t" + nodes.get(i);
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
				counter++;
				fwCount.write(output + "\n");
//				System.out.println("coloc " + coloc + " output " + output);
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of coloc hyperedges " + counter);
				
		///////////////////////////////////////////////////////////////////////////////////////////////
		
				
		// save hyperedges for coLocOrg
		minHyperedgeSize = 2;
		visHyperedges = new HashSet<String>(coLocOrg_attacks.size());
		
		counter = 0;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_coLocOrg);
			
			int coloc;
			String output;
			itr_coloc = coLocOrg_attacks.keySet().iterator();
			
			while (itr_coloc.hasNext()) {
				coloc = itr_coloc.next();
				
				nodes = new ArrayList<Integer>();
				nodes.add(coloc);
				
				itr_attacks = coLocOrg_attacks.get(coloc).iterator();
				while(itr_attacks.hasNext()) {
					nodes.add(itr_attacks.next());
				}
				
				if (nodes.size() < minHyperedgeSize) continue;
				
				Collections.sort(nodes);
				
				output = nodes.get(0) + "";
				for (int i = 1; i < nodes.size(); i ++) {
					output = output + "\t" + nodes.get(i);
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
				counter++;
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of colocorg hyperedges " + counter);
				
		////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//save node_groundTruth
		try {
			FileWriter fwCount = new FileWriter(fileOutput_node_groundTruth);
			
			Iterator<Integer> itr_nodes = nodes_Labels.keySet().iterator();
			int node;
			String label;
			while(itr_nodes.hasNext()) {
				node = itr_nodes.next();
				label = nodes_Labels.get(node);
				fwCount.write(node + "\t" + label + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
		save();
	}
}
