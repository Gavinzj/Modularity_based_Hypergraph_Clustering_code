package preprocess.stringDB;

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

public class StringDBParser {
	
	static String fileInput = FilePath_Mon.filePathPre + "/nba_pbp_2001-full/511145.protein.links.full.v11.5.txt";
	static String fileInput_cluster = FilePath_Mon.filePathPre + "/nba_pbp_2001-full/511145.clusters.proteins.v11.5.txt";
	
	static int yearInterval = 1;
	static int minYear, maxYear;
	static ArrayList<Integer> yearsSort = new ArrayList<Integer>();
	static HashMap<String, Integer> venues = new HashMap<String, Integer>();
	static HashMap<String, Integer> authors = new HashMap<String, Integer>();
	static ArrayList<String> publications = new ArrayList<String>();
	
	public static void read() throws IOException {
		
		int proteinCnt = 0;
		HashMap<String, Integer> proteinIDs = new HashMap<String, Integer>();
		
		String[] strs;
		String line;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			
			bufferedReader.readLine();
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split(" ");
				
				String protein1 = strs[0];
				String protein2 = strs[1];
				int cooccurrence = Integer.parseInt(strs[5]);
				int coexpression = Integer.parseInt(strs[7]);
				
				if (cooccurrence > 0 && coexpression > 0) {
					if (!proteinIDs.containsKey(protein1)) proteinIDs.put(protein1, proteinCnt++);
					if (!proteinIDs.containsKey(protein2)) proteinIDs.put(protein2, proteinCnt++);
				}
			}
		}
		
		System.out.println("# of protein " + proteinCnt);
		List<Integer>[] cointeractsWith = new ArrayList[proteinCnt];
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			cointeractsWith[proteinID] = new ArrayList<Integer>();
		}
		
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			
			bufferedReader.readLine();
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split(" ");
				
				String protein1 = strs[0];
				String protein2 = strs[1];
				int cooccurrence = Integer.parseInt(strs[5]);
				int coexpression = Integer.parseInt(strs[7]);
				
				if (cooccurrence > 0 && coexpression > 0) {
					int protein1ID = proteinIDs.get(protein1);
					int protein2ID = proteinIDs.get(protein2);
					cointeractsWith[protein1ID].add(protein2ID);
				}
			}
		}
		
		List<Integer>[] clusters = new ArrayList[proteinCnt];
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			clusters[proteinID] = new ArrayList<Integer>();
		}
		
		try (FileReader reader = new FileReader(fileInput_cluster);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			
			bufferedReader.readLine();
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				
				String clusterStr = strs[1].replace("CL:", "");
				int cluster = Integer.parseInt(clusterStr);
				String protein = strs[2];
				
				if (proteinIDs.containsKey(protein)) {
					int proteinID = proteinIDs.get(protein);
					clusters[proteinID].add(cluster);
				}
			}
		}
		
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			Collections.sort(clusters[proteinID]);
		}
		
		int clusterCnt = 0;
		HashMap<Integer, Integer> distinctClusters = new HashMap<Integer, Integer>();
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			for (int cluster : clusters[proteinID]) {
				if (!distinctClusters.containsKey(cluster)) distinctClusters.put(cluster, clusterCnt++);
			}
		}
		System.out.println("# of clusters " + distinctClusters.size() + " " + clusterCnt);
		
		// for each cluster
		double[] idfs = new double[clusterCnt];
		Iterator<Integer> itr_cluster = distinctClusters.keySet().iterator();
		while (itr_cluster.hasNext()) {
			int thisCluster = itr_cluster.next();
			int thisClusterID = distinctClusters.get(thisCluster);
					
			// for each protein
			double idf = 0;
			for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
				if (clusters[proteinID].contains(thisCluster)) idf++;
			}
			idf = Math.log((double) proteinCnt / idf) / Math.log(2);
			idfs[thisClusterID] = idf;
		}
		
		int[] singleClusterIDs = new int[proteinCnt];
		
		// for each protein
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			List<Integer> thisClusters = clusters[proteinID];
			double clusterSize = thisClusters.size();
			double tf = 1.0 / clusterSize;
			
			double maxTfidf = -1;
			int maxClusterID = -1;
			for (int cluster : thisClusters) {
				int clusterID = distinctClusters.get(cluster);
				double idf = idfs[clusterID];
				double tfidf = tf * idf;
				
				if (maxTfidf < tfidf) {
					maxTfidf = tfidf;
					maxClusterID = clusterID;
				}
			}
			
			singleClusterIDs[proteinID] = maxClusterID;
		}
		
		HashSet<Integer>[] proteinsInClusters = new HashSet[distinctClusters.size()];
		for (int clusterIdx = 0; clusterIdx < distinctClusters.size(); clusterIdx++) {
			proteinsInClusters[clusterIdx] = new HashSet<Integer>();
		}
		
		// for each protein
		for (int proteinID = 0; proteinID < proteinCnt; proteinID++) {
			int clusterID = singleClusterIDs[proteinID];
			proteinsInClusters[clusterID].add(proteinID);
		}
		
		int singleClusterCnt = 0;
		for (int clusterID = 0; clusterID < distinctClusters.size(); clusterID++) {
			if (proteinsInClusters[clusterID].size() <= 0) continue;
			singleClusterCnt++;
//			System.out.println("cluster " + clusterID + " proteins " + proteinsInClusters[clusterID]);
		}
		System.out.println("# of single clusters " + singleClusterCnt);
 	}
	
	public static void constructGraph(int lowYear, int topYear) {
		String fileOutput_publications = FilePath_Mon.filePathPre + "/publications_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_" + lowYear + "_" + topYear + ".txt";
		
		HashSet<Integer> authors = new HashSet<Integer>();
		HashMap<Integer, List<Integer>> nodeLabels = new HashMap<Integer, List<Integer>>();
		ArrayList<String> hyperedges = new ArrayList<String>();
		
		String[] strs;
		int venue, author;
		for (String publication : publications) {
			strs = publication.split("\t");
			venue = Integer.parseInt(strs[0]);
			
//			if (year == topYear) {
			for (int i = 2; i < strs.length; i++) {
				author = Integer.parseInt(strs[i]);
				authors.add(author);
				if (!nodeLabels.containsKey(author)) nodeLabels.put(author, new ArrayList<Integer>());
				nodeLabels.get(author).add(venue);
			}
//			}
			
			String hyperedge = strs[2] + "";
			for (int i = 3; i < strs.length; i++) {
				hyperedge = hyperedge + "\t" + strs[i];
			}
			hyperedges.add(hyperedge);
		}
		
		System.out.println("authors " + authors.size() + " labeled author " + nodeLabels.keySet().size());
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_publications);
			HashSet<String> visHyperedges = new HashSet<String>();
			String hyperedge;
			for (int i = 0; i < hyperedges.size(); i++) {
				hyperedge = hyperedges.get(i);
				if (!visHyperedges.add(hyperedge)) continue;
				
				fw_user.write(hyperedges.get(i) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	////////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_node_groundTruth);
			
			String output = "";
			Iterator<Integer> itr_author = nodeLabels.keySet().iterator();
			List<Integer> itr_label;
			while(itr_author.hasNext()) {
				author = itr_author.next();
				
				output = author + "";
				itr_label = nodeLabels.get(author);
				for (int label : itr_label) {
					output = output + "\t" + label;
				}
				fw_user.write(output + "\n");
			}
			fw_user.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws Exception {
		read();
	}
}