package statistic;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import graph.Edge;
import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.DoubleMergeSort;
import utilities.FilePath_Mon;
import utilities.SingleMergeSort;

public class Statistic_Clustering {
	
	static String fileOutputPre;
	static String fileInput;
	static String dataset;
	static double ratio;
	static double eps = Math.pow(10, -5);
	static double lambda;
	static double gamma;
	static int[] edgeInfos;
	
	// Ground Truth Community
	// [cluster]
	static int[][] clusters_Ground;
	static int[] clusters;
	static int classifiedNums;
	static double avgClusterSize;
	
	public static void loadGroundTruth(int maxClusterSize, double ratio) throws IOException {
		
		Statistic_Clustering.ratio = ratio;
		dataset = Hypergraph.dataset;
		
		if (!Constant.CONNECTED) fileInput = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		else fileInput = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		if (!Constant.CONNECTED) fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "gtclusterInfo_";
		else fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "gtclusterInfo_";

		load(maxClusterSize);
	}
	
	public static void loadDiscover(int trial, int maxClusterSize, String moveStrategy, boolean toHigherOrder, String ordering, double ratio)
			throws IOException, InterruptedException {
		
		Statistic_Clustering.ratio = ratio;
		dataset = Hypergraph.dataset;
		
		String method = "Louvain_" + moveStrategy + "_ordered_" + toHigherOrder + "_order_" + ordering + "_ratio_" + ratio + "_trial_" + trial;
		fileInput = FilePath_Mon.filePathPre + "/clustering/louvain/node_cluster_" + method + ".txt";
		
		if (!Constant.CONNECTED) fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "dvclusterInfo_" + method + "_";
		else fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "dvclusterInfo_" + method + "_";
		
		load(maxClusterSize);
	}
	
	public static void loadBaselineDiscover(int trial, int maxClusterSize, String method, double ratio)
			throws IOException, InterruptedException {
		
		Statistic_Clustering.ratio = ratio;
		dataset = Hypergraph.dataset;
		
		fileInput = FilePath_Mon.filePathPre + "/clustering/baseline/node_cluster_" + method + "_trial_" + trial + ".txt";
		
		if (!Constant.CONNECTED) fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "dvclusterInfo_" + method + "_";
		else fileOutputPre = FilePath_Mon.filePathPre + "/stat/" + "dvclusterInfo_" + method + "_";
		
		load(maxClusterSize);
	}
	
	public static void load(int maxClusterSize) throws IOException {
		
		ArrayList<ArrayList<Integer>> clustersArr = new ArrayList<ArrayList<Integer>>();

		ArrayList<Integer> nodes;
		String[] strs;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line;
			while ((line = bufferedReader.readLine()) != null) {
				
				// process each line
				strs = line.split("\t");
				if (strs.length <= maxClusterSize) continue;
				
				nodes = new ArrayList<Integer>(strs.length);
				for (int i = 0; i < strs.length; i++) {
					nodes.add(Integer.parseInt(strs[i]));
				}
				clustersArr.add(nodes);
			}
		}
		
		avgClusterSize = 0;
		int[] clusterSizes = new int[clustersArr.size()];
		
		clusters_Ground = new int[clustersArr.size()][];
		clusters = new int[Hypergraph.getVertexSize()];
		Arrays.fill(clusters, -1);
		
		int v;
		int clusterSize;
		for (int clusterIdx = 0; clusterIdx < clustersArr.size(); clusterIdx ++) {
			
			nodes = clustersArr.get(clusterIdx);
			clusterSize = nodes.size();
			
			clusterSizes[clusterIdx] = clusterSize;
			avgClusterSize += clusterSize;
			
			clusters_Ground[clusterIdx] = new int[clusterSize];
			
			for (int i = 0; i < clusterSize; i++) {
				v = nodes.get(i);
				clusters_Ground[clusterIdx][i] = v;
				clusters[v] = clusterIdx;
			}
		}
		
		avgClusterSize /= clustersArr.size();
		SingleMergeSort SSort = new SingleMergeSort();
		SSort.sort(clusterSizes);
		
		classifiedNums = 0;
		for (int i = 0; i < clusters.length; i++) {
			if (clusters[i] != -1) classifiedNums++;
		}
		
		System.out.println("cluster number " + clusters_Ground.length + " classified node " + classifiedNums);
		System.out.println("cluster size, avg " + avgClusterSize + " min " + clusterSizes[0] + " max " + clusterSizes[clusterSizes.length - 1] + 
				" mid " + clusterSizes[0 + (clusterSizes.length-0)/2]);
	}

	public static void getInfo() throws IOException, InterruptedException {
//		getNodeDegree();
//		System.out.println();
//		getEdgeCardinality();
//		System.out.println();
		getHModularity();
		System.out.println();
		getBogModularity();
		System.out.println();
//		getDensity();
//		System.out.println();
//		getOverlapness();
//		System.out.println();
//		get2DegOfNodePair();
//		System.out.println();
//		get2Homogeneity();
//		System.out.println();
//		getBoarderEdgeInfo();
//		System.out.println();
	}
	
	public static void getNodeDegree() {
		
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		int clusterNum = clusters_Ground.length;
		
		int[] followArr = new int[clusterNum];
		double[] baseArr = new double[clusterNum];
		for (int clusterID = 0; clusterID < clusterNum; clusterID++) {
			followArr[clusterID] = clusterID;
			baseArr[clusterID] = clusters_Ground[clusterID].length;
		}
		
		DoubleMergeSort Dsort = new DoubleMergeSort();
		Dsort.sort(baseArr, followArr);
		
		String output = "";
		
		for (int i = clusterNum - 1; i > 0; i--) {
			int clusterID = followArr[i];
			int clusterSize = clusters_Ground[clusterID].length;
			int[] nodes = clusters_Ground[clusterID];
			int[] degrees = new int[clusterSize];
			
			// for each node
			for (int j = 0; j < clusterSize; j++) {
				int node = nodes[j];
				
				// for each incident edge
				int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[node];
				if (node == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
				else secondIdx_INC = Hypergraph.INC_head[node + 1];
				
				int degree = secondIdx_INC - firstIdx_INC;
				degrees[j] = degree;
			}
			
			SingleMergeSort Ssort = new SingleMergeSort();
			Ssort.sort(degrees);
			
			int minDegree = degrees[0];
			int maxDegree = degrees[degrees.length - 1];
			int midDegree = degrees[0 + (degrees.length)/2];
			double avgDegree = 0;
			for (int j = 0; j < degrees.length; j++) avgDegree += degrees[j];
			avgDegree /= degrees.length;
			
			output = output + clusterID + "," + clusterSize + "," + minDegree + "," + maxDegree + "," + midDegree + "," + String.format("%.1f",avgDegree) + "\n";
//			System.out.println(clusterID + " size " + clusterSize + " degree min " + minDegree + " max " + maxDegree + " mid " + midDegree + " avg " + avgDegree);
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "nodeDegree.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "nodeDegree.txt");
			}
			
			fwCount.write("clusterID,size,min,max,mid,avg" + "\n");
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void getEdgeCardinality() {
		
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		int clusterNum = clusters_Ground.length;
		
		int[] followArr = new int[clusterNum];
		double[] baseArr = new double[clusterNum];
		for (int clusterID = 0; clusterID < clusterNum; clusterID++) {
			followArr[clusterID] = clusterID;
			baseArr[clusterID] = clusters_Ground[clusterID].length;
		}
		
		DoubleMergeSort Dsort = new DoubleMergeSort();
		Dsort.sort(baseArr, followArr);
		
		String output = "";
		
		for (int i = clusterNum - 1; i > 0; i--) {
			int thisCluster = followArr[i];
			int clusterSize = clusters_Ground[thisCluster].length;
			int[] nodes = clusters_Ground[thisCluster];
			
//			System.out.println("cluster " + thisCluster + " size " + clusterSize);
			
			// for each node
			HashSet<Integer> distinctIncidentEdges = new HashSet<Integer>();
			for (int j = 0; j < clusterSize; j++) {
				int node = nodes[j];
				
				// for each incident edge
				int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[node];
				if (node == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
				else secondIdx_INC = Hypergraph.INC_head[node + 1];
				
				for (int k = firstIdx_INC; k < secondIdx_INC; k++) {
					int edgeID = Hypergraph.INC_eID[k];
					distinctIncidentEdges.add(edgeID);
				}
			}
			
			List<Double> inFractions = new ArrayList<Double>(distinctIncidentEdges.size());
			List<Double> outFractions = new ArrayList<Double>(distinctIncidentEdges.size());
			List<Double> inCardinalitys = new ArrayList<Double>(distinctIncidentEdges.size());
			List<Double> outCardinalitys = new ArrayList<Double>(distinctIncidentEdges.size());
			
			// for each incident edge
			int cnt = 0;
			Iterator<Integer> itr_edge = distinctIncidentEdges.iterator();
			while (itr_edge.hasNext()) {
				int edgeID = itr_edge.next();
				
//				System.out.print("edge " + edgeID + ": ");
				
				// for each node in the edge
				int secondIdx_EINC, firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				double cardinality = secondIdx_EINC - firstIdx_EINC;
				double inFraction = 0;
				double outFraction = 0;
				for (int j = firstIdx_EINC; j < secondIdx_EINC; j++) {
					int node = Hypergraph.EINC_nID[j];
					int cluster = clusters[node];
					
//					System.out.print(node + "(" + cluster + "),");
					
					if (cluster != thisCluster) outFraction++;
					else inFraction++;
				}
//				System.out.println();
				
				inFraction /= cardinality;
				outFraction /= cardinality;
				
				if (inFraction == 1) {
					inFractions.add(inFraction);
					inCardinalitys.add(cardinality);
//					System.out.println("infraction " + inFraction + " cardinality " + cardinality);
				}
				
				if (outFraction > 0) {
					outFractions.add(outFraction);
					outCardinalitys.add(cardinality);
//					System.out.println("outFraction " + outFraction + " cardinality " + cardinality);
				}
				
				cnt++;
			}
			
			// in
			double minCardinality_in, maxCardinality_in, midCardinality_in, avgCardinality_in;
			if (inCardinalitys.size() > 0) {
				Collections.sort(inCardinalitys);
				minCardinality_in = inCardinalitys.get(0);
				maxCardinality_in = inCardinalitys.get(inCardinalitys.size() - 1);
				midCardinality_in = inCardinalitys.get(0 + (inCardinalitys.size())/2);
				avgCardinality_in = 0;
				for (int j = 0; j < inCardinalitys.size(); j++) avgCardinality_in += inCardinalitys.get(j);
				avgCardinality_in /= inCardinalitys.size();
				
			} else {
				minCardinality_in = 0;
				maxCardinality_in = 0;
				midCardinality_in = 0;
				avgCardinality_in = 0;
			}
				
//			System.out.println(thisCluster + " size " + clusterSize);
//			System.out.println("in-fraction size " + inFractions.size() + " cardinality min " + String.format("%.1f",minCardinality_in) + 
//					" max " + String.format("%.1f",maxCardinality_in) + " mid " + String.format("%.1f",midCardinality_in) + " avg " + String.format("%.1f",avgCardinality_in));
			
			// out
			double minFraction_out, maxFraction_out, midFraction_out, avgFraction_out;
			double minCardinality_out, maxCardinality_out, midCardinality_out, avgCardinality_out;
			if (outFractions.size() > 0) {
				Collections.sort(inFractions);
				minFraction_out = outFractions.get(0);
				maxFraction_out = outFractions.get(outFractions.size() - 1);
				midFraction_out = outFractions.get(0 + (outFractions.size())/2);
				avgFraction_out = 0;
				for (int j = 0; j < outFractions.size(); j++) avgFraction_out += outFractions.get(j);
				avgFraction_out /= outFractions.size();
				
				Collections.sort(outCardinalitys);
				minCardinality_out = outCardinalitys.get(0);
				maxCardinality_out = outCardinalitys.get(outCardinalitys.size() - 1);
				midCardinality_out = outCardinalitys.get(0 + (outCardinalitys.size())/2);
				avgCardinality_out = 0;
				for (int j = 0; j < outCardinalitys.size(); j++) avgCardinality_out += outCardinalitys.get(j);
				avgCardinality_out /= outCardinalitys.size();
				
			} else {
				minFraction_out = 0;
				maxFraction_out = 0;
				midFraction_out = 0;
				avgFraction_out = 0;
				
				minCardinality_out = 0;
				maxCardinality_out = 0;
				midCardinality_out = 0;
				avgCardinality_out = 0;
			}
			
//			System.out.println("out-fraction size " + outFractions.size() + " min " + String.format("%.1f",minFraction_out) + 
//					" max " + String.format("%.1f",maxFraction_out) + " mid " + String.format("%.1f",midFraction_out) + " avg " + String.format("%.1f",avgFraction_out));
//			System.out.println("in-fraction size " + outFractions.size() + " cardinality min " + String.format("%.1f",minCardinality_out) + 
//					" max " + String.format("%.1f",maxCardinality_out) + " mid " + String.format("%.1f",midCardinality_out) + " avg " + String.format("%.1f",avgCardinality_out));
			
			output = output + thisCluster + "," + clusterSize + 
					"," + inFractions.size() + "," + String.format("%.1f",minCardinality_in) + "," + String.format("%.1f",maxCardinality_in) + "," + String.format("%.1f",midCardinality_in) + "," + String.format("%.1f",avgCardinality_in) + 
					"," + outFractions.size() + "," + String.format("%.1f",minFraction_out) + "," + String.format("%.1f",maxFraction_out) + "," + String.format("%.1f",midFraction_out) + "," + String.format("%.1f",avgFraction_out) +
					"," + String.format("%.1f",minCardinality_out) + "," + String.format("%.1f",maxCardinality_out) + "," + String.format("%.1f",midCardinality_out) + "," + String.format("%.1f",avgCardinality_out) + "\n";
			
//			System.out.println();
//			System.out.println();
		}
		

		////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "edgeCardinality.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "edgeCardinality.txt");
			}
			
			fwCount.write("clusterID,size,inNum,minC,maxC,midC,avgC,outNum,minF,maxF,midF,avgF,minC,maxC,midC,avgC" + "\n");
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
		
	public static void getHModularity() {
		
		String output = "";
		
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		
		// edge weight
		double totalEdgeWeight = 0;
		double[] edge_weights = new double[m];
		
		double weight;
		int first_idx, second_idx;
		for (int edgeID = 0; edgeID < m; edgeID++) {
			// for each node in the edge
			first_idx = Hypergraph.EINC_head[edgeID];
			second_idx = Hypergraph.getSecondIdx_EINC(edgeID);
			
			// calculate the hyperedge weight
			weight = 0;
			for (int k = first_idx; k < second_idx; k++) {
				weight += Hypergraph.EINC_weight[k];
			}
			totalEdgeWeight += weight;
			edge_weights[edgeID] = weight;
		}
		
//		System.out.println("totalEdgeWeight " + totalEdgeWeight);
		
		// node degree
		double[] node_degrees = new double[n];
    	
    	// for each node
    	int curID = 0;
    	for (curID = 0; curID < n - 1; curID++) {
    		first_idx = Hypergraph.INC_head[curID];
    		second_idx = Hypergraph.INC_head[curID + 1];
			
			// for each incident hyperedge
			weight = 0;
			for (int i = first_idx; i < second_idx; i++) {
				weight += Hypergraph.INC_weight[i];
			}
			node_degrees[curID] = weight;
    	}
    	first_idx = Hypergraph.INC_head[curID];
		second_idx = Hypergraph.INC_weight.length;
		weight = 0;
		for (int i = first_idx; i < second_idx; i++) {
			weight += Hypergraph.INC_weight[i];
		}
		node_degrees[curID] = weight;
		
		// gamma, lambda
        edgeInfos = new int[m];
        int mMinusOne = m - 1;
        for (int edgeID = 0; edgeID < m; edgeID++) {
			if (edgeID == mMinusOne) edgeInfos[edgeID] = (Hypergraph.EINC_nID.length - Hypergraph.EINC_head[edgeID]);
			else edgeInfos[edgeID] = (Hypergraph.EINC_head[edgeID + 1] - Hypergraph.EINC_head[edgeID]);
		}
        
        if (!getGamma(m)) System.out.println("error!");
		
        ////////////////////////////////////////////////////////////////////////////////////
        
		// cluster volumn
		double[] cluster_volumn = new double[n];
		HashSet<Integer> distinctClusterIDs = new HashSet<Integer>(n);
		for (int nodeID = 0; nodeID < n; nodeID++) {
			double degree = node_degrees[nodeID];
			int cluster = clusters[nodeID];
			
			distinctClusterIDs.add(cluster);
			cluster_volumn[cluster] += degree;
		}
		
		// self_loop_weight
		double[] self_loop_weights = new double[n];
		for (int edgeID = 0; edgeID < m; edgeID++) {
			
			// for each node in hyperedge
			first_idx = Hypergraph.EINC_head[edgeID];
			if (edgeID == m - 1) second_idx = Hypergraph.EINC_nID.length;
			else second_idx = Hypergraph.EINC_head[edgeID + 1];
			
			int firstNode = Hypergraph.EINC_nID[first_idx];
			int firstCluster = clusters[firstNode];
			
//			System.out.println("first node " + firstNode + " cluster " + firstCluster);
			
			boolean isInnerEdge = true;
			for (int i = first_idx; i < second_idx; i++) {
				int node = Hypergraph.EINC_nID[i];
				int cluster = clusters[node];
				
				if (cluster != firstCluster) {
					isInnerEdge = false;
					break;
				}
			}
			
			if (isInnerEdge) {
				double edgeWeight = edge_weights[edgeID];
				self_loop_weights[firstCluster] += edgeWeight;
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////
		
		double Q = 0;
		Iterator<Integer> itr_cluster = distinctClusterIDs.iterator();
		while (itr_cluster.hasNext()) {
			int thisCluster = itr_cluster.next();
			
			double vol_c = cluster_volumn[thisCluster];
			double gamma_c = gamma * vol_c;
			Q += (((self_loop_weights[thisCluster] + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c)));
		}
		
		System.out.println("Q " + Q);
		
		output = "H-modularity," + String.format("%.4f", Q) + "\n" + output;
		
		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "Hmodularity.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "Hmodularity.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void getBogModularity() {
		
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		
		// edge weight
		double totalEdgeWeight = 0;
		int[] edge_cardinalities = new int[m];
		double[] edge_weights = new double[m];
		
		double weight;
		int first_idx, second_idx, cardinality;
		for (int edgeID = 0; edgeID < m; edgeID++) {
			// for each node in the edge
			first_idx = Hypergraph.EINC_head[edgeID];
			second_idx = Hypergraph.getSecondIdx_EINC(edgeID);
			
			// cardinality
			cardinality = second_idx - first_idx;
			edge_cardinalities[edgeID] = cardinality;
			
			// edge weight
			weight = 0;
			for (int i = first_idx; i < second_idx; i++) weight += Hypergraph.EINC_weight[i];
			totalEdgeWeight += weight;
			edge_weights[edgeID] = weight;
		}
//		System.out.println("totalEdgeWeight " + totalEdgeWeight);
		
		// node degree
		double[] node_degrees = new double[n];
    	
    	// for each node
    	int curID = 0;
    	for (curID = 0; curID < n - 1; curID++) {
    		// for each incident hyperedge
    		first_idx = Hypergraph.INC_head[curID];
    		second_idx = Hypergraph.INC_head[curID + 1];
			
			weight = 0;
			for (int i = first_idx; i < second_idx; i++) weight += Hypergraph.INC_weight[i];
			node_degrees[curID] = weight;
    	}
    	first_idx = Hypergraph.INC_head[curID];
		second_idx = Hypergraph.INC_weight.length;
		weight = 0;
		for (int i = first_idx; i < second_idx; i++) weight += Hypergraph.INC_weight[i];
		node_degrees[curID] = weight;
		
        ////////////////////////////////////////////////////////////////////////////////////
        
		// cardinality distribution
		int distinctCardinalityNum = 0;
		HashMap<Integer, Integer> distinctCardinalities = new HashMap<Integer, Integer>();
		for (int edgeID = 0; edgeID < m; edgeID++) {
			cardinality = edge_cardinalities[edgeID];
			if (!distinctCardinalities.containsKey(cardinality)) distinctCardinalities.put(cardinality, distinctCardinalityNum++);
		}
		
//		System.out.println("# distinct cardinality " + distinctCardinalityNum);
		
		int[] cardinalityFrequencies = new int[distinctCardinalityNum];
		for (int edgeID = 0; edgeID < m; edgeID++) {
			cardinality = edge_cardinalities[edgeID];
			int idx = distinctCardinalities.get(cardinality);
			cardinalityFrequencies[idx]++;
			
//			System.out.println("cardinality " + cardinality + " frequency " + cardinalityFrequencies[idx]);
		}
		
		// check correctness
		int edgeNum = 0;
		for (int i = 0; i < cardinalityFrequencies.length; i++) edgeNum += cardinalityFrequencies[i];
		if (edgeNum != m) System.out.println("#wrong " + edgeNum + " " + m);

		// cluster volumn
		double totalVolumn = 0;
		double[] cluster_volumn = new double[n];
		HashSet<Integer> distinctClusterIDs = new HashSet<Integer>(n);
		for (int nodeID = 0; nodeID < n; nodeID++) {
			double degree = node_degrees[nodeID];
			int cluster = clusters[nodeID];
			
			distinctClusterIDs.add(cluster);
			cluster_volumn[cluster] += degree;
			totalVolumn += degree;
		}
		
//		System.out.println("total volumn " + totalVolumn);
//		for (int i = 0; i < cluster_volumn.length; i++) {
//			if (cluster_volumn[i] > 0) {
//				System.out.println("clsuter " + i + " vol " + cluster_volumn[i]);
//			}
//		}
		
		// self_loop_weight
		double[] self_loop_weights = new double[n];
		for (int edgeID = 0; edgeID < m; edgeID++) {
			
			// for each node in hyperedge
			first_idx = Hypergraph.EINC_head[edgeID];
			if (edgeID == m - 1) second_idx = Hypergraph.EINC_nID.length;
			else second_idx = Hypergraph.EINC_head[edgeID + 1];
			
			int firstNode = Hypergraph.EINC_nID[first_idx];
			int firstCluster = clusters[firstNode];
			
//					System.out.println("first node " + firstNode + " cluster " + firstCluster);
			
			boolean isInnerEdge = true;
			for (int i = first_idx; i < second_idx; i++) {
				int node = Hypergraph.EINC_nID[i];
				int cluster = clusters[node];
				
				if (cluster != firstCluster) {
					isInnerEdge = false;
					break;
				}
			}
			
			if (isInnerEdge) {
				double edgeWeight = edge_weights[edgeID];
				self_loop_weights[firstCluster] += edgeWeight;
			}
		}
        ////////////////////////////////////////////////////////////////////////////////////
		
		double actual = 0;
		Iterator<Integer> itr_cluster = distinctClusterIDs.iterator();
		while (itr_cluster.hasNext()) {
			int thisCluster = itr_cluster.next();
			double self_loop_weight = self_loop_weights[thisCluster];
			
			actual += self_loop_weight;
		}
		
//		System.out.println("sum actual " + actual);
		
		// expected
		double expected = 0;
		
		// for each cardinality
		Iterator<Integer> itr_cardinality = distinctCardinalities.keySet().iterator();
		while (itr_cardinality.hasNext()) {
			cardinality = itr_cardinality.next();
			int idx = distinctCardinalities.get(cardinality);
			double frequency = cardinalityFrequencies[idx];
			
//			System.out.println("cardinality " + cardinality + " frequency " + frequency);
			
			double sum_pow_d = 0;
			
			// for each cluster
			itr_cluster = distinctClusterIDs.iterator();
			while (itr_cluster.hasNext()) {
				int thisCluster = itr_cluster.next();
				double vol = cluster_volumn[thisCluster];
				double ratio = vol / totalVolumn;
				double ratio_pow_d = Math.pow(ratio, cardinality);
				
//				System.out.println("cluster " + thisCluster + " vol " + vol + " ratio " + ratio + " ratio_pow_d " + ratio_pow_d);
				
				sum_pow_d += ratio_pow_d;
			}
			
//			System.out.println("sum_pow_d " + sum_pow_d);
			
			expected += frequency * sum_pow_d;
			
//			System.out.println("expected " + (frequency * sum_pow_d));
		}
		
//		System.out.println("sum expected " + expected);
		
		double Q = (1 / (double) m) * (actual - expected);
		
		String output = "Bog-modularity," + String.format("%.4f",Q) + "\n";
		
		System.out.println("Bog-modularity " + Q);
		

		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "Bogmodularity.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "Bogmodularity.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static boolean getGamma(int edgeSize) {
		
		// get cardinality distribution
		try {
			FileWriter fwCount = new FileWriter(FilePath_Mon.filePathPre + "/stat/cardinality/cardinality_trial_" + 0 + ".txt");
			for (int edgeID = 0; edgeID < edgeSize; edgeID++) {
				if (edgeInfos[edgeID] >= 2) fwCount.write(edgeInfos[edgeID] + "\n");
				edgeInfos[edgeID] = -1;
			}
			fwCount.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
		
		// ref: https://www.edureka.co/community/358/how-to-execute-a-python-file-with-few-arguments-in-java
		Process process = null;
		try {
			process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/getLambda.py " 
					+ "/data/" + dataset + "/stat/cardinality/cardinality_trial_" + 0 + ".txt" + " " + 2);
		} catch (Exception e) {
			System.out.println("Exception Raised" + e.toString());
		}
		InputStream stdout = process.getInputStream();
		BufferedReader reader = new BufferedReader(new InputStreamReader(stdout, StandardCharsets.UTF_8));
		String line;
		try {
			
			// read the first line
			line = reader.readLine();
			if (line == null) return false;
			
			// read the second line
			line = reader.readLine();
			if (line != null && !line.equals("nan")){
				lambda = Double.parseDouble(line);
				gamma = Math.exp(-1 * lambda);
				
				if (gamma < eps) return false;
				System.out.println(lambda + " " + gamma);
			} else return false;
		} catch (IOException e) {
			System.out.println("Exception in reading output" + e.toString());
		}
				
		return true;
	}
	
	public static void getDensity() {
		// refer to https://arxiv.org/pdf/2101.07480.pdf
		
		String output = "";
		
		int secondIdx_INC, edgeID, firstIdx_EINC, secondIdx_EINC;
		int m = Hypergraph.getEdgeSize();
		int n = Hypergraph.getVertexSize();
		
		double[] density = new double[n];
		double[] inDensity = new double[n];
		Arrays.fill(inDensity, -1);
		double[] outDensity = new double[n];
		Arrays.fill(outDensity, -1);
		
		double avgDensity = 0;
		double avgInDensity = 0;
		double avgOutDensity = 0;
		
		double inDensityCnt = 0;
		double outDensityCnt = 0;
		
		// for each node
		for (int nodeID = 0; nodeID < n; nodeID++) {
			
			HashSet<Integer> incidentEdges = new HashSet<Integer>();
			HashSet<Integer> nodes = new HashSet<Integer>();
			
			HashSet<Integer> inCluster_incidentEdges = new HashSet<Integer>();
			HashSet<Integer> inCluster_nodes = new HashSet<Integer>();
			
			HashSet<Integer> outCluster_incidentEdges = new HashSet<Integer>();
			HashSet<Integer> outCluster_nodes = new HashSet<Integer>();
			
			int clusterID = clusters[nodeID];
			
			// for each incident edge
			if (nodeID == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
			
			for (int j = Hypergraph.INC_head[nodeID]; j < secondIdx_INC; j++) {
				edgeID = Hypergraph.INC_eID[j];
				
				incidentEdges.add(edgeID);
				
//				System.out.print(" for edge " + edgeID + " (");
				
				// for each node in the edge
				firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				double cutoffSize = ratio * (double) (secondIdx_EINC - firstIdx_EINC);
				double fraction = 0;
				
				for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
//					System.out.print(Hypergraph.EINC_nID[k] + ",");
					
					if (clusterID == clusters[Hypergraph.EINC_nID[k]]) fraction++;
					
					nodes.add(Hypergraph.EINC_nID[k]);
				}
				
//				System.out.println(")");
//				System.out.println(incidentEdges.size() + " " + neighbors.size());
				
				if ((fraction + 1) >= cutoffSize) {
//					System.out.println("in");
					
					inCluster_incidentEdges.add(edgeID);
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) inCluster_nodes.add(Hypergraph.EINC_nID[k]);
					
				} else {
//					System.out.println("out");
					
					outCluster_incidentEdges.add(edgeID);
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) outCluster_nodes.add(Hypergraph.EINC_nID[k]);
				}
			}
			
//			System.out.println(inCluster_incidentEdges.size() + " " + inCluster_neighbors.size());
//			System.out.println(outCluster_incidentEdges.size() + " " + outCluster_neighbors.size());
			
			avgDensity += ((double) incidentEdges.size() / (double) nodes.size());
			density[nodeID] = ((double) incidentEdges.size() / (double) nodes.size());
			
			if (inCluster_nodes.size() > 0) {
				avgInDensity += ((double) inCluster_incidentEdges.size() / (double) inCluster_nodes.size());
				inDensityCnt++;
				
				inDensity[nodeID] = ((double) inCluster_incidentEdges.size() / (double) inCluster_nodes.size());
			}
			
			if (outCluster_nodes.size() > 0) {
				avgOutDensity += ((double) outCluster_incidentEdges.size() / (double) outCluster_nodes.size());
				outDensityCnt++;
				
				outDensity[nodeID] = ((double) outCluster_incidentEdges.size() / (double) outCluster_nodes.size());
			}
		}
		
		avgDensity /= n;
		avgInDensity /= inDensityCnt;
		avgOutDensity /= outDensityCnt;
		
		System.out.println("avg node density " + String.format("%.4f",avgDensity) + " in-denstity " + String.format("%.4f",avgInDensity) + 
				" out-density " + String.format("%.4f",avgOutDensity));
		
		output = output + "avgDensity," + String.format("%.4f",avgDensity) + ",avgInDensity," +String.format("%.4f",avgInDensity) 
			+ ",avgOutDensity," +String.format("%.4f",avgOutDensity) + "\n";
		
		////////////////////////////////////////////////////////////////////////////////////////////////

		// for each cluster
		for (int clusterID = 0; clusterID < clusters_Ground.length; clusterID++) {
			
			avgDensity = 0;
			avgInDensity = 0;
			avgOutDensity = 0;
			
			inDensityCnt = 0;
			outDensityCnt = 0;
			
			System.out.print("#cluster " + clusterID + " size " + clusters_Ground[clusterID].length + " (");
//			for (int node : clusters_Ground[clusterID]) System.out.print(node + ",");
			System.out.println(")");
			
			output = output + "cluster," + clusterID + ",size," + clusters_Ground[clusterID].length;
			
			HashSet<Integer> incidentEdges = new HashSet<Integer>();
			HashSet<Integer> neighbors = new HashSet<Integer>();
			
			// for each node in cluster
			for (int i = 0; i < clusters_Ground[clusterID].length; i++) {
				int nodeID = clusters_Ground[clusterID][i];
				
//				System.out.println("for node " + nodeID);
				
				avgDensity += density[nodeID];
				if (inDensity[nodeID] != -1) {
					avgInDensity += inDensity[nodeID];
					inDensityCnt++;
				}
				if (outDensity[nodeID] != -1) {
					avgOutDensity += outDensity[nodeID];
					outDensityCnt++;
				}
				
				// for each incident edge
				if (nodeID == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
				else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
				
				for (int j = Hypergraph.INC_head[nodeID]; j < secondIdx_INC; j++) {
					edgeID = Hypergraph.INC_eID[j];
					
					incidentEdges.add(edgeID);
					
//					System.out.print(" for edge " + edgeID + " (");
					
					// for each node in the edge
					firstIdx_EINC = Hypergraph.EINC_head[edgeID];
					if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
					else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
					
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
//						System.out.print(Hypergraph.EINC_nID[k] + ",");
						
						neighbors.add(Hypergraph.EINC_nID[k]);
					}
					
//					System.out.println();
//					System.out.println(incidentEdges.size() + " " + neighbors.size());
				}
			}
			
			double clusterDensity = (double) incidentEdges.size() / (double) neighbors.size();
			
			avgDensity /= (double) clusters_Ground[clusterID].length;
			
			if (inDensityCnt > 0) {
				avgInDensity /= inDensityCnt;
			} else {
				avgInDensity = -1;
			}
			
			if (outDensityCnt > 0) {
				avgOutDensity /= outDensityCnt;
			} else {
				avgOutDensity = -1;
			}
			
			System.out.println("cluster density " + String.format("%.4f",clusterDensity) + " avg " + String.format("%.4f",avgDensity) + 
					" in-density " + String.format("%.4f",avgInDensity) + " out-density " + String.format("%.4f",avgOutDensity));
			
			output = output + ",clusterDensity," + String.format("%.4f",clusterDensity) + ",avgDensity," + String.format("%.4f",avgDensity) +
				",inDensity," + String.format("%.4f",avgInDensity) + ",outDensity," + String.format("%.4f",avgOutDensity) + "\n";
		}
		

		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "density.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "density.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void getOverlapness() {
		
		String output = "";
		
		int secondIdx_INC, edgeID, firstIdx_EINC, secondIdx_EINC;
		int m = Hypergraph.getEdgeSize();
		int n = Hypergraph.getVertexSize();
		
		double[] overlapness = new double[n];
		double[] inOverlapness = new double[n];
		Arrays.fill(inOverlapness, -1);
		double[] outOverlapness = new double[n];
		Arrays.fill(outOverlapness, -1);
		
		double avgOverlapness = 0;
		double avgInOverlapness = 0;
		double avgOutOverlapness = 0;
		
		double inOverlapnessCnt = 0;
		double outOverlapnessCnt = 0;
		
		// for each node
		for (int nodeID = 0; nodeID < n; nodeID++) {
			
			double totalEdgeSize = 0;
			HashSet<Integer> nodes = new HashSet<Integer>();
			
			double inTotalEdgeSize = 0;
			HashSet<Integer> inCluster_nodes = new HashSet<Integer>();
			
			double outTotalEdgeSize = 0;
			HashSet<Integer> outCluster_nodes = new HashSet<Integer>();
			
			int clusterID = clusters[nodeID];
			
			// for each incident edge
			if (nodeID == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
			
			for (int j = Hypergraph.INC_head[nodeID]; j < secondIdx_INC; j++) {
				edgeID = Hypergraph.INC_eID[j];
				
				// for each node in the edge
				firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				int edgeSize = secondIdx_EINC - firstIdx_EINC;
				
				double cutoffSize = ratio * (double) edgeSize;
				double fraction = 0;
				
				totalEdgeSize += edgeSize;
				for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
					if (clusterID == clusters[Hypergraph.EINC_nID[k]]) fraction++;
					
					nodes.add(Hypergraph.EINC_nID[k]);
				}
				
				if ((fraction + 1) >= cutoffSize) {
					inTotalEdgeSize += edgeSize;
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) inCluster_nodes.add(Hypergraph.EINC_nID[k]);
				} else {
					outTotalEdgeSize += edgeSize;
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) outCluster_nodes.add(Hypergraph.EINC_nID[k]);
				}
			}
			
			avgOverlapness += (totalEdgeSize / (double) nodes.size());
			overlapness[nodeID] = (totalEdgeSize / (double) nodes.size());
			
			if (inCluster_nodes.size() > 0) {
				avgInOverlapness += (inTotalEdgeSize / (double) inCluster_nodes.size());
				inOverlapnessCnt++;
				
				inOverlapness[nodeID] = (inTotalEdgeSize / (double) inCluster_nodes.size());
			}
			
			if (outCluster_nodes.size() > 0) {
				avgOutOverlapness += (outTotalEdgeSize / (double) outCluster_nodes.size());
				outOverlapnessCnt++;
				
				outOverlapness[nodeID] = (outTotalEdgeSize / (double) outCluster_nodes.size());
			}
		}
		
		avgOverlapness /= n;
		avgInOverlapness /= inOverlapnessCnt;
		avgOutOverlapness /= outOverlapnessCnt;
		
		System.out.println("avg node overlapness " + String.format("%.4f",avgOverlapness) + " in-overlapness " + String.format("%.4f",avgInOverlapness) + 
				" out-overlapness " + String.format("%.4f",avgOutOverlapness));
		
		output = output + "avgOverlap," + String.format("%.4f",avgOverlapness) + ",avgInOverlap," +String.format("%.4f",avgInOverlapness) 
		+ ",avgOutOverlap," +String.format("%.4f",avgOutOverlapness) + "\n";
		
		////////////////////////////////////////////////////////////////////////////////////////////////

		// for each cluster
		for (int clusterID = 0; clusterID < clusters_Ground.length; clusterID++) {
			
			avgOverlapness = 0;
			avgInOverlapness = 0;
			avgOutOverlapness = 0;
			
			inOverlapnessCnt = 0;
			outOverlapnessCnt = 0;
			
			System.out.print("#cluster " + clusterID + " size " + clusters_Ground[clusterID].length + " (");
//			for (int node : clusters_Ground[clusterID]) System.out.print(node + ",");
			System.out.println(")");
			
			output = output + "cluster," + clusterID + ",size," + clusters_Ground[clusterID].length;
			
			// for each node in cluster
			for (int i = 0; i < clusters_Ground[clusterID].length; i++) {
				int nodeID = clusters_Ground[clusterID][i];
				
				avgOverlapness += overlapness[nodeID];
				if (inOverlapness[nodeID] != -1) {
					avgInOverlapness += inOverlapness[nodeID];
					inOverlapnessCnt++;
				}
				if (outOverlapness[nodeID] != -1) {
					avgOutOverlapness += outOverlapness[nodeID];
					outOverlapnessCnt++;
				}
			}
			
			avgOverlapness /= (double) clusters_Ground[clusterID].length;
			if (inOverlapnessCnt > 0) {
				avgInOverlapness /= inOverlapnessCnt;
			} else {
				avgInOverlapness = -1;
			}
			if (outOverlapnessCnt > 0) {
				avgOutOverlapness /= outOverlapnessCnt;
			} else {
				avgOutOverlapness = -1;
			}
			
			System.out.println("cluster overlap avg " + String.format("%.4f",avgOverlapness) + 
					" in-overlap " + String.format("%.4f",avgInOverlapness) + " out-overlap " + String.format("%.4f",avgOutOverlapness));
			
			output = output + ",avgOverlap," + String.format("%.4f",avgOverlapness) +
					",inOverlap," + String.format("%.4f",avgInOverlapness) + ",outOverlap," + String.format("%.4f",avgOutOverlapness) + "\n";
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "overlapness.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "overlapness.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void get2DegOfNodePair() {
        
		String output = "";
		
		int secondIdx_INC, edgeID, neighbor, firstIdx_EINC, secondIdx_EINC;
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		
		int[] clusterIdx = new int[n];
        for (int i = 0; i < n; i++) clusterIdx[i] = -1;
        
		double[] degrees = new double[n];
		double[] inDegrees = new double[n];
		Arrays.fill(inDegrees, -1);
		double[] outDegrees = new double[n];
		Arrays.fill(outDegrees, -1);
		
		double avgDegree = 0;
		double avgInDegree = 0;
		double avgOutDegree = 0;
		
		double inDegreeCnt = 0;
		double outDegreeCnt = 0;
		
		// for each node
		for (int nodeID = 0; nodeID < n; nodeID++) {
			
			int clusterID = clusters[nodeID];
			
			// for each incident edge
			if (nodeID == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
			
			HashSet<Integer> distinctNeighbors = new HashSet<Integer>();
			for (int j = Hypergraph.INC_head[nodeID]; j < secondIdx_INC; j++) {
				edgeID = Hypergraph.INC_eID[j];
				
				// for each node in the edge
				firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
					if (nodeID == Hypergraph.EINC_nID[k]) continue;
					distinctNeighbors.add(Hypergraph.EINC_nID[k]);
				}
			}
			
			int distinctNeighborCnt = 0;
			Iterator<Integer> itr = distinctNeighbors.iterator();
			while (itr.hasNext()) {
				clusterIdx[itr.next()] = distinctNeighborCnt++;
			}
			
			int[] degreeWithNeighbor = new int[distinctNeighborCnt];
			for (int j = Hypergraph.INC_head[nodeID]; j < secondIdx_INC; j++) {
				edgeID = Hypergraph.INC_eID[j];
				
				// for each node in the edge
				firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
					neighbor = Hypergraph.EINC_nID[k];
					if (nodeID == neighbor) continue;
					
					degreeWithNeighbor[clusterIdx[neighbor]]++;
				}
			}
			
			double degree = 0;
			double inDegree = 0;
			double outDegree = 0;
			
			double inCnt = 0;
			double outCnt = 0;
			
			// for each neighbor
			itr = distinctNeighbors.iterator();
			while (itr.hasNext()) {
				neighbor = itr.next();
				degree += degreeWithNeighbor[clusterIdx[neighbor]];
				
				if (clusterID == clusters[neighbor]) {
					inDegree += degreeWithNeighbor[clusterIdx[neighbor]];
					inCnt++;
				} else {
					outDegree += degreeWithNeighbor[clusterIdx[neighbor]];
					outCnt++;
				}
				
				clusterIdx[neighbor] = -1;
			}
			
			avgDegree += (degree / distinctNeighborCnt);
			degrees[nodeID] = (degree / distinctNeighborCnt);
			
			if (inCnt > 0) {
				avgInDegree += (inDegree / inCnt);
				inDegreeCnt++;
				
				inDegrees[nodeID] = (inDegree / inCnt);
			}
			
			if (outCnt > 0) {
				avgOutDegree += (outDegree / outCnt);
				outDegreeCnt++;
				
				outDegrees[nodeID] = (outDegree / outCnt);
			}
		}
		
		avgDegree /= n;
		avgInDegree /= inDegreeCnt;
		avgOutDegree /= outDegreeCnt;
		
		System.out.println("2-degree, avg " + String.format("%.4f",avgDegree) + " in-degree " + String.format("%.4f",avgInDegree) + 
				" out-degree " + String.format("%.4f",avgOutDegree));
		
		output = output + "avgDegree," + String.format("%.4f",avgDegree) + ",avgInDegree," +String.format("%.4f",avgInDegree) 
		+ ",avgOutDegree," +String.format("%.4f",avgOutDegree) + "\n";
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// for each cluster
		for (int clusterID = 0; clusterID < clusters_Ground.length; clusterID++) {
			
			avgDegree = 0;
			avgInDegree = 0;
			avgOutDegree = 0;
			
			inDegreeCnt = 0;
			outDegreeCnt = 0;
			
			System.out.print("#cluster " + clusterID + " size " + clusters_Ground[clusterID].length + " (");
//			for (int node : clusters_Ground[clusterID]) System.out.print(node + ",");
			System.out.println(")");
			
			output = output + "cluster," + clusterID + ",size," + clusters_Ground[clusterID].length;
			
			// for each node in cluster
			for (int i = 0; i < clusters_Ground[clusterID].length; i++) {
				int nodeID = clusters_Ground[clusterID][i];
				
				avgDegree += degrees[nodeID];
				if (inDegrees[nodeID] != -1) {
					avgInDegree += inDegrees[nodeID];
					inDegreeCnt++;
				}
				if (outDegrees[nodeID] != -1) {
					avgOutDegree += outDegrees[nodeID];
					outDegreeCnt++;
				}
			}
			
			avgDegree /= (double) clusters_Ground[clusterID].length;
			
			if (inDegreeCnt > 0) avgInDegree /= inDegreeCnt;
			else avgInDegree = -1;
			
			if (outDegreeCnt > 0) avgOutDegree /= outDegreeCnt;
			else avgOutDegree = -1;
			
			System.out.println("cluster 2-degree avg " + String.format("%.4f",avgDegree) + 
					" in-degree " + String.format("%.4f",avgInDegree) + " out-degree " + String.format("%.4f",avgOutDegree));
			
			output = output + ",avg2Degree," + String.format("%.4f",avgDegree) +
					",in2Degree," + String.format("%.4f",avgInDegree) + ",out2Degree," + String.format("%.4f",avgOutDegree) + "\n";
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "2DegreeNode.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "2DegreeNode.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void get2Homogeneity() {
		
		String output = "";
		
		int firstIdx_EINC, secondIdx_EINC;
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		
		int[] clusterIdx = new int[n];
        for (int i = 0; i < n; i++) clusterIdx[i] = -1;
		
        double avgHomogeneity = 0;
        double[] homogeneities = new double[m];
        
		// for each edge
		for (int edgeID = 0; edgeID < m; edgeID++) {
			
			// for each node in the edge
			firstIdx_EINC = Hypergraph.EINC_head[edgeID];
			if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
			
			double degreeSum = 0;
			double pairCnt = 0;
			
			int node1, node2;
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				node1 = Hypergraph.EINC_nID[i];
				
				// for each incident edge
				int secondIdx_INC;
				if (node1 == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
				else secondIdx_INC = Hypergraph.INC_head[node1 + 1];
				
				HashSet<Integer> distinctNeighbors = new HashSet<Integer>();
				for (int j = Hypergraph.INC_head[node1]; j < secondIdx_INC; j++) {
					int incidentEdgeID = Hypergraph.INC_eID[j];
					
					// for each node in the edge
					firstIdx_EINC = Hypergraph.EINC_head[incidentEdgeID];
					if (incidentEdgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
					else secondIdx_EINC = Hypergraph.EINC_head[incidentEdgeID + 1];
					
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
						if (node1 == Hypergraph.EINC_nID[k]) continue;
						distinctNeighbors.add(Hypergraph.EINC_nID[k]);
					}
				}
				
				int distinctNeighborCnt = 0;
				Iterator<Integer> itr = distinctNeighbors.iterator();
				while (itr.hasNext()) {
					clusterIdx[itr.next()] = distinctNeighborCnt++;
				}
				
				int[] degreeWithNeighbor = new int[distinctNeighborCnt];
				for (int j = Hypergraph.INC_head[node1]; j < secondIdx_INC; j++) {
					int incidentEdgeID = Hypergraph.INC_eID[j];
					
					// for each node in the edge
					firstIdx_EINC = Hypergraph.EINC_head[incidentEdgeID];
					if (incidentEdgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
					else secondIdx_EINC = Hypergraph.EINC_head[incidentEdgeID + 1];
					
					for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
						int neighbor = Hypergraph.EINC_nID[k];
						if (node1 == neighbor) continue;
						
						degreeWithNeighbor[clusterIdx[neighbor]]++;
					}
				}
				
				///////////////////////////////////////////////////////////////////////
				
				for (int j = i + 1; j < secondIdx_EINC; j++) {
					node2 = Hypergraph.EINC_nID[j];
					pairCnt++;
					
					if (distinctNeighbors.contains(node2)){
						degreeSum += degreeWithNeighbor[clusterIdx[node2]];
					}
				}
				
				////////////////////////////////////////////////////////////////////////////
				
				itr = distinctNeighbors.iterator();
				while (itr.hasNext()) {
					clusterIdx[itr.next()] = -1;
				}
			}
			
			degreeSum /= pairCnt;
//			System.out.println("degreeSum " + degreeSum);
			
			avgHomogeneity += degreeSum;
					
			homogeneities[edgeID] = degreeSum;
		}
		
		avgHomogeneity /= m;
		
		System.out.println("2-homogeneity avg " + String.format("%.4f",avgHomogeneity));
		
		output = output + "avg2Homogeneity," + String.format("%.4f",avgHomogeneity) + "\n";
		
		/////////////////////////////////////////////////////////////////////////////////
		
		// for each cluster
		for (int clusterID = 0; clusterID < clusters_Ground.length; clusterID++) {
			
			System.out.print("#cluster " + clusterID + " size " + clusters_Ground[clusterID].length + " (");
			System.out.println(")");
			
			output = output + "cluster," + clusterID + ",size," + clusters_Ground[clusterID].length;
			
			double innerHomogeneity = 0;
			double incidentHomogeneity = 0;
			double connectHomogeneity = 0;
			
			double innerCnt = 0;
			double incidentCnt = 0;
			double connectCnt = 0;
			
			for (int edgeID = 0; edgeID < m; edgeID++) {
				
				int second_idx;
				if (edgeID == m - 1) second_idx = Hypergraph.EINC_nID.length;
				else second_idx = Hypergraph.EINC_head[edgeID + 1];
				
				int edgeSize = second_idx - Hypergraph.EINC_head[edgeID];
				double cutoffSize = 1 * (double) edgeSize;
				int nodeCntInEdge = 0;
				
				// for each node in hyperedge
				for (int j = Hypergraph.EINC_head[edgeID]; j < second_idx; j++) {
					if (clusters[Hypergraph.EINC_nID[j]] == clusterID) nodeCntInEdge++;
				}
				
				if (nodeCntInEdge > 0) {
					incidentHomogeneity += homogeneities[edgeID];
					incidentCnt++;
					
					if (nodeCntInEdge >= cutoffSize) {
						innerHomogeneity += homogeneities[edgeID];
						innerCnt++;
					} else {
						connectHomogeneity += homogeneities[edgeID];
						connectCnt++;
					}
				}
			}
			
			if (incidentCnt > 0) incidentHomogeneity /= incidentCnt;
			else incidentHomogeneity = -1;
			
			if (innerCnt > 0) innerHomogeneity /= innerCnt;
			else innerHomogeneity = -1;
			
			if (connectCnt > 0) connectHomogeneity /= connectCnt;
			else connectHomogeneity = -1;
			
			System.out.println("cluster 2-homogeneity incident " + String.format("%.4f",incidentHomogeneity) + 
					" inner " + String.format("%.4f",innerHomogeneity) + " connect " + String.format("%.4f",connectHomogeneity));
			
			output = output + ",incident," + String.format("%.4f",incidentHomogeneity) +
					",inner," + String.format("%.4f",innerHomogeneity) + ",connect," + String.format("%.4f",connectHomogeneity) + "\n";
		}
		

		////////////////////////////////////////////////////////////////////////////////////////////////
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "homogeneity.txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "homogeneity.txt");
			}
			
			fwCount.write(output);
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void getBoarderEdgeInfo() {
		// get the information of boarder edge that how many clusters share a boarder edge
		
		int firstIdx_EINC, secondIdx_EINC;
		int m = Hypergraph.getEdgeSize();
		int n = Hypergraph.getVertexSize();
		
		int boarderEdgeNum = 0;
		List<Double> maxRatios = new ArrayList<Double>();
		
		// for each edge
		for (int edgeID = 0; edgeID < m; edgeID++) {
			
			HashMap<Integer, Integer> nodeFrequency = new HashMap<Integer, Integer>();
			
			// for each node in the edge
			firstIdx_EINC = Hypergraph.EINC_head[edgeID];
			if (edgeID == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
			
			double cardinality = secondIdx_EINC - firstIdx_EINC;
			
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				int clusterID = clusters[Hypergraph.EINC_nID[i]];
				
				int frequency;
				if (nodeFrequency.containsKey(clusterID)) {
					frequency = nodeFrequency.get(clusterID);
					frequency++;
				} else frequency = 1;
				nodeFrequency.put(clusterID, frequency);
			}
			
			if (nodeFrequency.size() > 1) {
//				System.out.println("edge " + edgeID + ": ");
//				for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
//					System.out.println(Hypergraph.EINC_nID[i] + " - " + clusters[Hypergraph.EINC_nID[i]]);
//				}
				
				double maxRatio = -1;
				Iterator<Integer> itr = nodeFrequency.keySet().iterator();
				while (itr.hasNext()) {
					int clusterID = itr.next();
					double frequency = nodeFrequency.get(clusterID);
					
//					System.out.println("cluster " + clusterID + " frequency " + frequency);
					
					if (maxRatio < (frequency / cardinality)) maxRatio = (frequency / cardinality);
				}
				
				boarderEdgeNum++;
				maxRatios.add(maxRatio);
			}
		}
		
		double maxRatio = -1;
		double minRatio = Constant.large;
		double avgRatio = 0;
		for (double ratio : maxRatios) {
			avgRatio += maxRatio;
			if (maxRatio < ratio) maxRatio = ratio;
			if (minRatio > ratio) minRatio = ratio;
		}
		avgRatio /= maxRatios.size();
		
		System.out.println("# of boarder edge " + boarderEdgeNum);
		System.out.println("avg ratio " + avgRatio + " max " + maxRatio + " minRatio " + minRatio);
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
		boolean toHigherOrder = false;
		String ordering = "randomOrder";
		String moveStrategy = "move";
		int maxClusterSize = 0;
		double ratio = 0.5;
		String baselineMethod = "BogCNMRan";
		
		loadGroundTruth(maxClusterSize, ratio);
		getInfo();
		
		for (int trial = 0; trial < 3; trial++) {
			System.out.println("trial " + trial);
			loadDiscover(trial, maxClusterSize, moveStrategy, toHigherOrder, ordering, ratio);
			getInfo();
		}
		
//		for (int trial = 0; trial < 3; trial++) {
//			System.out.println("trial " + trial);
//			loadBaselineDiscover(trial, maxClusterSize, baselineMethod, ratio);
//			getInfo();
//		}
	}
}
