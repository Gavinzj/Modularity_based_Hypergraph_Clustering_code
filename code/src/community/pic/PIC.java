package community.pic;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.DoubleMergeSort;
import utilities.FilePath_Mon;

public class PIC {
	
	int trial;
	double ratio;
	boolean toHigherOrder;
	String ordering;
	
	double increase_total;
	double increasePerPass;
	double Q_prev;
	double Q_curr;
	int n_aggregations = 10;
	
	double lambda;
	double gamma = -1;
	double initGamma;
	double eps = Constant.EPSILON;
	
	int n, m;
	double[] node_degrees;
	
    double totalEdgeWeight;
    int[] INC_eID;
	int[] INC_head;
    int INC_eID_length;
	
	int[] EINC_nID;
	double[] EINC_weight;
	int[] EINC_head;
    int EINC_nID_length;
	double[] edge_weights;
	
	int[] ADJ_nHead;
	int[] ADJ_nID;
	int ADJ_nID_length;
	
	List<Integer>[] nodesInCluster;
	int[] global_cluster;
	int[] clusters;
	double[] cluster_volume;
	
	int[] nodeInOrder;
	int[] nodePriority;
	
	// variables for move()
	double[] neighbor_clusters_weights;
	int[] clusterIdx;
	int[] fractions;
	int[] neighbors;
	int[] nbrClusterIDs;
	
	// variables for rebuildGraph
	List<Integer> distinctIDs;
	List<Integer> newEINC_heads;
	List<Integer> newEINC_nIDs;
	List<Double> newEINC_weight;
	List<Integer>[] newINC_eIDs;
	double[] cluster_weights;
	double[] self_loop_weights;
	int new_n;
	
	// variables for updateGamma
	int[] edgeInfos;
	Process process;
	
	Iterator<Integer> itr_int;
	Iterator<Double> itr_double;
	DoubleMergeSort Dsorter = new DoubleMergeSort();
	
	boolean save;
	String moveStrategy;
	String dataset;
	
	String toPrint;
	int interestedNode = 1203;
	
	public PIC (int trial, boolean toHigherOrder, String ordering, double ratio, boolean save) throws IOException {
		this.trial = trial;
		this.toHigherOrder = toHigherOrder;
		this.ordering = ordering;
		this.ratio = ratio;
		this.save = save;
		this.dataset = Hypergraph.dataset;
	}
	
	public PIC (int trial, boolean toHigherOrder, String ordering, double ratio, double lambda, boolean save) throws IOException {
		this.trial = trial;
		this.toHigherOrder = toHigherOrder;
		this.ordering = ordering;
		this.ratio = ratio;
		this.save = save;
		this.dataset = Hypergraph.dataset;
		this.lambda = lambda;
	}
	
	public void initPart1() {
        
		n = Hypergraph.getVertexSize();
		m = Hypergraph.getEdgeSize();
		
		EINC_head = new int[m];
		for (int i = 0; i < m; i ++) {
			EINC_head[i] = Hypergraph.EINC_head[i];
		}
		
		INC_head = new int[n];
		for (int i = 0; i < n; i ++) {
			INC_head[i] = Hypergraph.INC_head[i];
		}
		
		int cardinality = Hypergraph.EINC_nID.length;
		EINC_nID = new int[cardinality];
		EINC_weight = new double[cardinality];
		INC_eID = new int[cardinality];
		for (int i = 0; i < cardinality; i++) {
			EINC_nID[i] = Hypergraph.EINC_nID[i];
			EINC_weight[i] = Hypergraph.EINC_weight[i];
			INC_eID[i] = Hypergraph.INC_eID[i];
		}
		
		EINC_nID_length = cardinality;
		INC_eID_length = cardinality;
	}
	
	public void initPart2() {
		
		totalEdgeWeight = 0;
		edge_weights = new double[m];
		
		double weight;
		int first_idx, second_idx;
		for (int edgeID = 0; edgeID < m; edgeID++) {
			// for each node in the edge
			first_idx = EINC_head[edgeID];
			second_idx = getSecondIdx_EINC(edgeID);
			
			// calculate the hyperedge weight
			weight = 0;
			for (int k = first_idx; k < second_idx; k++) {
				weight += EINC_weight[k];
			}
			totalEdgeWeight += weight;
			edge_weights[edgeID] = weight;
		}
    	
    	node_degrees = new double[n];
    	
    	// for each node
    	int curID = 0;
    	for (curID = 0; curID < n - 1; curID++) {
    		first_idx = INC_head[curID];
			second_idx = INC_head[curID + 1];
			
			// for each incident hyperedge
			weight = 0;
			for (int i = first_idx; i < second_idx; i++) {
				weight += Hypergraph.INC_weight[i];
			}
			node_degrees[curID] = weight;
    	}
    	first_idx = INC_head[curID];
		second_idx = Hypergraph.INC_weight.length;
		weight = 0;
		for (int i = first_idx; i < second_idx; i++) {
			weight += Hypergraph.INC_weight[i];
		}
		node_degrees[curID] = weight;
		
		if (save) {
            global_cluster = new int[n];
            for (int i = 0; i < n; i++) {
            	global_cluster[i] = i;
            }
        }
    	
		///////////////////////////////////////////////////////////////////////////
		
        // variables for move
		neighbor_clusters_weights = new double[n];
        clusterIdx = new int[n];
        for (int i = 0; i < n; i++) clusterIdx[i] = -1;
        
        // variables for reconstruct
        cluster_weights = new double[n];
        
        // variables for updateGamma
        edgeInfos = new int[m];
        int mMinusOne = m - 1;
        for (int edgeID = 0; edgeID < m; edgeID++) {
			if (edgeID == mMinusOne) edgeInfos[edgeID] = (EINC_nID_length - Hypergraph.EINC_head[edgeID]);
			else edgeInfos[edgeID] = (Hypergraph.EINC_head[edgeID + 1] - Hypergraph.EINC_head[edgeID]);
		}
        
        if (!getGamma(m)) System.out.println("error!");
        
        initGamma = gamma;
        
        ////////////////////////////////////////////////////////////////////////////
        
        nodeInOrder = new int[n];
		clusters = new int[n];
		cluster_volume = new double[n];
		for (int i = 0; i < n; i++) {
			nodeInOrder[i] = i;
			clusters[i] = i;
			cluster_volume[i] += node_degrees[i];
		}
		
		if (toHigherOrder) {
			nodePriority = new int[n];
		}
	}
	
	public String move_AON(int thisNode) {
		
		int thisCluster = clusters[thisNode];
		
//		System.out.println();
//		System.out.println("move node " + thisNode + "(" + thisCluster + ")");
		
		// for each incident edge
		int secondIdx_INC;
		if (thisNode == n - 1) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];
		
		HashMap<Integer, Double> distinctIncidentEdges = new HashMap<Integer, Double>();
		HashMap<Integer, HashSet<Integer>> distinctEdgeClusters = new HashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> distinctNeighborClusters = new HashSet<Integer>();
		for (int i = INC_head[thisNode]; i < secondIdx_INC; i++) {
			int edgeID = INC_eID[i];
			
//			System.out.print("edgeID " + edgeID + ": ");
			
			// for each node in the edge
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			double weight = 0;
			for (int j = EINC_head[edgeID]; j < secondIdx; j++) {
				int neighbor = EINC_nID[j];
				int cluster = clusters[neighbor];
				
//				System.out.print(neighbor + "(" + cluster + "),");
				
				if (!distinctEdgeClusters.containsKey(cluster)) distinctEdgeClusters.put(cluster, new HashSet<Integer>());
				distinctEdgeClusters.get(cluster).add(edgeID);
				distinctNeighborClusters.add(cluster);
				weight += EINC_weight[j];
			}
			distinctIncidentEdges.put(edgeID, weight);
			
//			System.out.println(" weight " + weight);
		}
		
//		System.out.println("neighbor cluster " + distinctNeighborClusters);
//		Iterator<Integer> itr_cluster_1 = distinctEdgeClusters.keySet().iterator();
//		while (itr_cluster_1.hasNext()) {
//			int cluster = itr_cluster_1.next();
//			System.out.println("cluster " + cluster + " edges " + distinctEdgeClusters.get(cluster));
//		}	
		
		if (distinctNeighborClusters.size() == 1) return "";
		
		// try removing other neighbors from original cluster
		double vol = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double eta = (1 / (totalEdgeWeight - (gamma * (vol - thisDegree)))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * vol)));
		
		double self_loop_weight = 0;
		double incident_weight = 0;
		
		// iterate through incident edges
		Iterator<Integer> itr_edges = distinctIncidentEdges.keySet().iterator();
		while (itr_edges.hasNext()) {
			int edgeID = itr_edges.next();
			
			// for each node in the edge
			int firstIdx = EINC_head[edgeID];
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			boolean isInCluster = true;
			boolean isSelfLoop = false;
			if ((secondIdx - firstIdx) == 1 && EINC_nID[firstIdx] == thisNode) {
				isSelfLoop = true;
			} else {
				for (int i = firstIdx; i < secondIdx; i++) {
					int node = EINC_nID[i];
					
					if (clusters[node] != thisCluster) {
						isInCluster = false;
						break;
					}
				}
			}
			
			if (!isInCluster) {
//				System.out.print("invalid edge " + edgeID + ": ");
//				for (int i = firstIdx; i < secondIdx; i++) System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//				System.out.println();
				continue;
			}
			
			double weight = distinctIncidentEdges.get(edgeID);
			if (isSelfLoop) {
				self_loop_weight += weight;
//				System.out.print("self loop edge " + edgeID + ": ");
//				for (int i = firstIdx; i < secondIdx; i++) System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//				System.out.println();
			} else {
				incident_weight += weight;
//				System.out.print("incident edge " + edgeID + ": ");
//				for (int i = firstIdx; i < secondIdx; i++) System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//				System.out.println();
			}
		}
		
		double deltaQ_exit = -1 * ((incident_weight) / totalEdgeWeight) - ((eta * totalEdgeWeight) - 1);
		
//		System.out.println("self_loop_weight " + self_loop_weight + " incident_weight " + incident_weight);
//		System.out.println("deltaQ_exit " + deltaQ_exit);
		
//		// put back
//		vol = cluster_volume[thisCluster] - thisDegree;
//		eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * (vol + thisDegree))));
//		double deltaQ_1 = (incident_weight / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
//		double deltaQ_local_1 = deltaQ_exit + deltaQ_1;
//		System.out.println("deltaQ " + deltaQ_1 + " deltaQ_local " + deltaQ_local_1);
				
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		// try moving to target clsuter
		List<Integer> neighborClustersSorted = new ArrayList<Integer>(distinctNeighborClusters.size());
		Iterator<Integer> itr_cluster = distinctNeighborClusters.iterator();
		while (itr_cluster.hasNext()) neighborClustersSorted.add(itr_cluster.next());
		Collections.sort(neighborClustersSorted);
		
		for (int targetClusterID : neighborClustersSorted){
			if (thisCluster == targetClusterID) continue;
			
//			System.out.println("move to cluster " + targetClusterID);
			
			double toTargetClusterWeight = 0;
			// iterate through incident edges
			Iterator<Integer> itr_edge = distinctEdgeClusters.get(targetClusterID).iterator();
			while (itr_edge.hasNext()) {
				int edgeID = itr_edge.next();
				boolean isInCluster = true;
				
				// for each node in the edge
				int secondIdx = getSecondIdx_EINC(edgeID);
				
				for (int i = EINC_head[edgeID]; i < secondIdx; i++) {
					int nodeID = EINC_nID[i];
					if (nodeID == thisNode) continue;
					
					int cluster = clusters[nodeID];
					if (cluster != targetClusterID) {
						isInCluster = false;
						break;
					}
				}
				
				if (isInCluster) {
					toTargetClusterWeight += distinctIncidentEdges.get(edgeID);
//					System.out.print("valid edge " + edgeID + " weight " + distinctIncidentEdges.get(edgeID) + ": ");
//					for (int i = EINC_head[edgeID]; i < secondIdx; i++) {
//						System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//					}
//					System.out.println();
				} else {
//					System.out.print("invalid edge " + edgeID + " weight " + distinctIncidentEdges.get(edgeID) + ": ");
//					for (int i = EINC_head[edgeID]; i < secondIdx; i++) {
//						System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//					}
//					System.out.println();
				}
			}
			
//			System.out.println("incident_weight " + incident_weight);
        	
			vol = cluster_volume[targetClusterID];
			eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * (vol + thisDegree))));
			double deltaQ = (toTargetClusterWeight / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
			double deltaQ_local = deltaQ + deltaQ_exit;
			if (deltaQ_local > deltaQ_best) {
        		deltaQ_best = deltaQ_local;
        		best_targetCluster = targetClusterID;
            }
			
//			System.out.println(targetClusterID + " deltaQ " + deltaQ + " deltaQ_local " + deltaQ_local);
//			System.out.println("");
		}
		
//		System.out.println("deltaQ_best " + deltaQ_best + " best_targetCluster " + best_targetCluster);
		
		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
        	// decide to move
        	increasePerPass += deltaQ_best;
        	cluster_volume[thisCluster] -= thisDegree;
        	cluster_volume[best_targetCluster] += thisDegree;
            clusters[thisNode] = best_targetCluster;
//            System.out.println("new cluster of " + thisNode + " is " + clusters[thisNode]);
			return deltaQ_best + "," + best_targetCluster;
		}
		
//		System.out.println("=====================\n");
		
		return "";
	}
	
	public String move_split(int thisNode) {
		
		int thisCluster = clusters[thisNode];
		
//		toPrint += "move node " + thisNode + "(" + thisCluster + ")" + "\n";
//		System.out.println("move node " + thisNode + "(" + thisCluster + ")");
		
		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == n - 1) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];
		
		int incidentEdgeNum = secondIdx_INC - firstIdx_INC;
		
		HashMap<Integer, Double> distinctIncidentEdges = new HashMap<Integer, Double>(incidentEdgeNum);
		HashMap<Integer, HashSet<Integer>> distinctEdgeClusters = new HashMap<Integer, HashSet<Integer>>();
		HashMap<Integer, Integer> distinctNeighborClusters = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> idxToCluster = new HashMap<Integer, Integer>();
		int neighborClusterNum = 0;
		for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
			int edgeID = INC_eID[i];
			
//			System.out.print("edgeID " + edgeID + ": ");
			
			// for each node in the edge
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			double weight = 0;
			for (int j = EINC_head[edgeID]; j < secondIdx; j++) {
				int neighbor = EINC_nID[j];
				int cluster = clusters[neighbor];
				
//				System.out.print(neighbor + "(" + cluster + "),");
				
				if (!distinctEdgeClusters.containsKey(cluster)) {
					distinctEdgeClusters.put(cluster, new HashSet<Integer>());
					idxToCluster.put(neighborClusterNum, cluster);
					distinctNeighborClusters.put(cluster, neighborClusterNum++);
				}
				distinctEdgeClusters.get(cluster).add(edgeID);
				
				weight += EINC_weight[j];
			}
			distinctIncidentEdges.put(edgeID, weight);
			
//			System.out.println(" weight " + weight);
		}
		
//		System.out.println("neighbor cluster # " + distinctNeighborClusters.size() + " " + neighborClusterNum);
//		System.out.println("neighbor cluster " + distinctNeighborClusters);
//		Iterator<Integer> itr_cluster_1 = distinctEdgeClusters.keySet().iterator();
//		while (itr_cluster_1.hasNext()) {
//			int cluster = itr_cluster_1.next();
//			System.out.println("cluster " + cluster + " edges " + distinctEdgeClusters.get(cluster));
//		}
		
		if (distinctNeighborClusters.size() == 1) return "";
		
		// try removing other neighbors from original cluster
		double vol = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double eta = (1 / (totalEdgeWeight - (gamma * (vol - thisDegree)))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * vol)));
		
		double incident_weight = 0;
		
		// iterate through incident edges
		Iterator<Integer> itr_edges = distinctIncidentEdges.keySet().iterator();
		while (itr_edges.hasNext()) {
			int edgeID = itr_edges.next();
			double edgeWeight = distinctIncidentEdges.get(edgeID);
			
			// for each node in the edge
			int firstIdx = EINC_head[edgeID];
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			if ((secondIdx - firstIdx) == 1 && EINC_nID[firstIdx] == thisNode) {
				
//				System.out.print("self loop edge " + edgeID + ": ");
//				for (int i = firstIdx; i < secondIdx; i++) System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//				System.out.println();
				
			} else {
//				System.out.print("incident edge " + edgeID + "(" + edgeWeight + "): ");
//				for (int i = firstIdx; i < secondIdx; i++) {
//					int node = EINC_nID[i];
//					int cluster = clusters[node];
//					double nodeWeight = EINC_weight[i];
//					System.out.print(node + "(" + cluster + ";" + nodeWeight + "),");
//				}
//				System.out.println();
				
				double thisWeight = 0;
				double[] fractions = new double[neighborClusterNum];
				for (int i = firstIdx; i < secondIdx; i++) {
					int node = EINC_nID[i];
					double nodeWeight = EINC_weight[i];
					if (node == thisNode) {
						thisWeight = EINC_weight[i];
						continue;
					}
					
					int cluster = clusters[node];
					int clusterIdx = distinctNeighborClusters.get(cluster);
					fractions[clusterIdx] += nodeWeight;
				}
				
//				for (int i = 0; i < fractions.length; i++) {
//					if (fractions[i] > 0) System.out.println("to cluster " + idxToCluster.get(i) + " fraction " + fractions[i] + "/" + edgeWeight + "=" + (fractions[i] / edgeWeight));
//				}
				
				// for the edge incident to nodes in current cluster
				int clusterIdx = distinctNeighborClusters.get(thisCluster);
				double fraction_u = thisWeight;
				double fraction_cNotu = fractions[clusterIdx];
				double fraction_c = fraction_u + fraction_cNotu;
				if (fractions[clusterIdx] > 0) {
					double existingContribution = adjustWeight_LinearLog(edgeWeight, fraction_c);
					double newContribution = adjustWeight_LinearLog(edgeWeight, fraction_u) + adjustWeight_LinearLog(edgeWeight, fraction_cNotu);
					double deltaContribution = existingContribution - newContribution;
					incident_weight += deltaContribution;
					
//					System.out.println("fraction_c " + fraction_c + " fraction_u " + fraction_u + " fraction_cNotu " + fraction_cNotu);
//					System.out.println("existingContribution " + existingContribution + " newContribution " + newContribution + " deltaContribution " + deltaContribution);
				}
			}
		}
		
		double deltaQ_exit = -1 * (incident_weight / totalEdgeWeight) - ((eta * totalEdgeWeight) - 1);
		
//		System.out.println("incident_weight " + incident_weight);
//		System.out.println("deltaQ_exit " + deltaQ_exit);
		
//		// put back
//		vol = cluster_volume[thisCluster] - thisDegree;
//		eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * (vol + thisDegree))));
//		double deltaQ_1 = (incident_weight / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
//		double deltaQ_local_1 = deltaQ_exit + deltaQ_1;
//		System.out.println("deltaQ_local " + deltaQ_local_1);

		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		// try moving to target clsuter
		List<Integer> neighborClustersSorted = new ArrayList<Integer>(distinctNeighborClusters.size());
		Iterator<Integer> itr_cluster = distinctNeighborClusters.keySet().iterator();
		while (itr_cluster.hasNext()) neighborClustersSorted.add(itr_cluster.next());
		Collections.sort(neighborClustersSorted);
		
		for (int targetClusterID : neighborClustersSorted){
			if (thisCluster == targetClusterID) continue;
			
//			System.out.println("move to cluster " + targetClusterID);
			
			double toTargetClusterWeight = 0;
			// iterate through incident edges
			Iterator<Integer> itr_edge = distinctEdgeClusters.get(targetClusterID).iterator();
			while (itr_edge.hasNext()) {
				int edgeID = itr_edge.next();
				double edgeWeight = distinctIncidentEdges.get(edgeID);
				
				// for each node in the edge
				int firstIdx = EINC_head[edgeID];
				int secondIdx = getSecondIdx_EINC(edgeID);
				
//				System.out.print("incident edge " + edgeID + "(" + edgeWeight + "): ");
//				for (int i = firstIdx; i < secondIdx; i++) {
//					int node = EINC_nID[i];
//					int cluster = clusters[node];
//					double nodeWeight = EINC_weight[i];
//					System.out.print(node + "(" + cluster + ";" + nodeWeight + "),");
//				}
//				System.out.println();
				
				double thisWeight = 0;
				double[] fractions = new double[neighborClusterNum];
				for (int i = firstIdx; i < secondIdx; i++) {
					int node = EINC_nID[i];
					double nodeWeight = EINC_weight[i];
					if (node == thisNode) {
						thisWeight = EINC_weight[i];
						continue;
					}
					
					int cluster = clusters[node];
					int clusterIdx = distinctNeighborClusters.get(cluster);
					fractions[clusterIdx] += nodeWeight;
				}
				
//				for (int i = 0; i < fractions.length; i++) {
//					System.out.println("to cluster " + idxToCluster.get(i) + " fraction " + fractions[i] + "/" + edgeWeight + "=" + (fractions[i] / edgeWeight));
//				}
				
				int clusterIdx = distinctNeighborClusters.get(targetClusterID);
				double fraction_u = thisWeight;
				double fraction_cNotu = fractions[clusterIdx];
				double fraction_c = fraction_u + fraction_cNotu;
				if (fractions[clusterIdx] > 0) {
					double existingContribution = adjustWeight_LinearLog(edgeWeight, fraction_u) + adjustWeight_LinearLog(edgeWeight, fraction_cNotu);
					double newContribution = adjustWeight_LinearLog(edgeWeight, fraction_c);
					double deltaContribution = newContribution - existingContribution;
					toTargetClusterWeight += deltaContribution;
					
//					System.out.println("fraction_c " + fraction_c + " fraction_u " + fraction_u + " fraction_cNotu " + fraction_cNotu);
//					System.out.println("existingContribution " + existingContribution + " newContribution " + newContribution + " deltaContribution " + deltaContribution);
				}
			}
			
//			System.out.println("toTargetClusterWeight " + toTargetClusterWeight);
        	
			vol = cluster_volume[targetClusterID];
			eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * thisDegree))) - (1 / (totalEdgeWeight - (gamma * (vol + thisDegree))));
			double deltaQ = (toTargetClusterWeight / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
			double deltaQ_local = deltaQ + deltaQ_exit;
			if (deltaQ_local > deltaQ_best) {
        		deltaQ_best = deltaQ_local;
        		best_targetCluster = targetClusterID;
            }
			
//			System.out.println(targetClusterID + " deltaQ " + deltaQ + " deltaQ_local " + deltaQ_local);
//			System.out.println("");
		}
		
//		System.out.println("deltaQ_best " + deltaQ_best + " best_targetCluster " + best_targetCluster);
		
		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
        	// decide to move
        	increasePerPass += deltaQ_best;
        	cluster_volume[thisCluster] -= thisDegree;
        	cluster_volume[best_targetCluster] += thisDegree;
            clusters[thisNode] = best_targetCluster;
//            System.out.println("new cluster of " + thisNode + " is " + clusters[thisNode]);
            return deltaQ_best + "," + best_targetCluster;
		}
		
//		System.out.println("=====================\n");
		
		return "";
	}
	
	public String move_split_adjust(int thisNode) {
		
		int thisCluster = clusters[thisNode];
		
//		toPrint += "move node " + thisNode + "(" + thisCluster + ")" + "\n";
//		System.out.println("move node " + thisNode + "(" + thisCluster + ")");
		
		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == n - 1) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];
		
		int incidentEdgeNum = secondIdx_INC - firstIdx_INC;
		
		HashMap<Integer, Double> distinctIncidentEdges = new HashMap<Integer, Double>(incidentEdgeNum);
		HashMap<Integer, HashSet<Integer>> distinctEdgeClusters = new HashMap<Integer, HashSet<Integer>>();
		HashMap<Integer, Integer> distinctNeighborClusters = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> idxToCluster = new HashMap<Integer, Integer>();
		int neighborClusterNum = 0;
		for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
			int edgeID = INC_eID[i];
			
//			System.out.print("edgeID " + edgeID + ": ");
			
			// for each node in the edge
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			double weight = 0;
			for (int j = EINC_head[edgeID]; j < secondIdx; j++) {
				int neighbor = EINC_nID[j];
				int cluster = clusters[neighbor];
				
//				System.out.print(neighbor + "(" + cluster + "),");
				
				if (!distinctEdgeClusters.containsKey(cluster)) {
					distinctEdgeClusters.put(cluster, new HashSet<Integer>());
					idxToCluster.put(neighborClusterNum, cluster);
					distinctNeighborClusters.put(cluster, neighborClusterNum++);
				}
				distinctEdgeClusters.get(cluster).add(edgeID);
				
				weight += EINC_weight[j];
			}
			distinctIncidentEdges.put(edgeID, weight);
			
//			System.out.println(" weight " + weight);
		}
		
//		System.out.println("neighbor cluster # " + distinctNeighborClusters.size() + " " + neighborClusterNum);
//		System.out.println("neighbor cluster " + distinctNeighborClusters);
//		Iterator<Integer> itr_cluster_1 = distinctEdgeClusters.keySet().iterator();
//		while (itr_cluster_1.hasNext()) {
//			int cluster = itr_cluster_1.next();
//			System.out.println("cluster " + cluster + " edges " + distinctEdgeClusters.get(cluster));
//		}
		
		if (distinctNeighborClusters.size() == 1) return "";
		
		// try removing other neighbors from original cluster
		double vol_cu, one_minus_gammaExpect_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		double one_minus_gammaExpect_cNou = 1 - (gamma * expect(vol_cNou, totalEdgeWeight, ratio));
		double one_minus_gammaExpect_u = 1 - (gamma * expect(thisDegree, totalEdgeWeight, ratio));
		double one_minus_gammaExpect_c = 1 - (gamma * expect(vol_c, totalEdgeWeight, ratio));
		double residule = (gamma * (ratio - 1)) -1;
		
		double eta = (1.0 / one_minus_gammaExpect_cNou) + (1 / one_minus_gammaExpect_u) - (1 / one_minus_gammaExpect_c);
				
		double incident_weight = 0;
		
		// iterate through incident edges
		Iterator<Integer> itr_edges = distinctIncidentEdges.keySet().iterator();
		while (itr_edges.hasNext()) {
			int edgeID = itr_edges.next();
			double edgeWeight = distinctIncidentEdges.get(edgeID);
			
			// for each node in the edge
			int firstIdx = EINC_head[edgeID];
			int secondIdx = getSecondIdx_EINC(edgeID);
			
			if ((secondIdx - firstIdx) == 1 && EINC_nID[firstIdx] == thisNode) {
				
//				System.out.print("self loop edge " + edgeID + ": ");
//				for (int i = firstIdx; i < secondIdx; i++) System.out.print(EINC_nID[i] + "(" + clusters[EINC_nID[i]] + "),");
//				System.out.println();
				
			} else {
//				System.out.print("incident edge " + edgeID + "(" + edgeWeight + "): ");
//				for (int i = firstIdx; i < secondIdx; i++) {
//					int node = EINC_nID[i];
//					int cluster = clusters[node];
//					double nodeWeight = EINC_weight[i];
//					System.out.print(node + "(" + cluster + ";" + nodeWeight + "),");
//				}
//				System.out.println();
				
				double thisWeight = 0;
				double[] fractions = new double[neighborClusterNum];
				for (int i = firstIdx; i < secondIdx; i++) {
					int node = EINC_nID[i];
					double nodeWeight = EINC_weight[i];
					if (node == thisNode) {
						thisWeight = EINC_weight[i];
						continue;
					}
					
					int cluster = clusters[node];
					int clusterIdx = distinctNeighborClusters.get(cluster);
					fractions[clusterIdx] += nodeWeight;
				}
				
//				for (int i = 0; i < fractions.length; i++) {
//					if (fractions[i] > 0) System.out.println("to cluster " + idxToCluster.get(i) + " fraction " + fractions[i] + "/" + edgeWeight + "=" + (fractions[i] / edgeWeight));
//				}
				
				// for the edge incident to nodes in current cluster
				int clusterIdx = distinctNeighborClusters.get(thisCluster);
				double fraction_u = thisWeight;
				double fraction_cNotu = fractions[clusterIdx];
				double fraction_c = fraction_u + fraction_cNotu;
				if (fraction_cNotu > 0) {
					double existingContribution = adjustWeight_LinearLog(edgeWeight, fraction_c);
					double newContribution = adjustWeight_LinearLog(edgeWeight, fraction_u) + adjustWeight_LinearLog(edgeWeight, fraction_cNotu);
					double deltaContribution = existingContribution - newContribution;
					incident_weight += deltaContribution;
					
//					System.out.println("fraction_c " + fraction_c + " fraction_u " + fraction_u + " fraction_cNotu " + fraction_cNotu);
//					System.out.println("existingContribution " + existingContribution + " newContribution " + newContribution + " deltaContribution " + deltaContribution);
				}
			}
		}
		
		double deltaQ_exit = (incident_weight / totalEdgeWeight) + eta + residule;
		deltaQ_exit = -1 * deltaQ_exit;
		
//		System.out.println("incident_weight " + incident_weight);
//		System.out.println("deltaQ_exit " + deltaQ_exit);
		
//		// put back
//		vol_cu = cluster_volume[thisCluster];
//		vol_c = vol_cu - thisDegree;
//		one_minus_gammaExpect_c = 1 - (gamma * expect(vol_c, totalEdgeWeight, ratio));
//		one_minus_gammaExpect_u = 1 - (gamma * expect(thisDegree, totalEdgeWeight, ratio));
//		one_minus_gammaExpect_cu = 1 - (gamma * expect(vol_cu, totalEdgeWeight, ratio));
//		residule = (gamma * (ratio - 1)) -1;
//		eta = (1.0 / one_minus_gammaExpect_c) + (1 / one_minus_gammaExpect_u) - (1 / one_minus_gammaExpect_cu);
//		double deltaQ_1 = (incident_weight / totalEdgeWeight) + eta + residule;
//		double deltaQ_local_1 = deltaQ_exit + deltaQ_1;
//		System.out.println("deltaQ_local " + deltaQ_local_1);

		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		// try moving to target clsuter
		List<Integer> neighborClustersSorted = new ArrayList<Integer>(distinctNeighborClusters.size());
		Iterator<Integer> itr_cluster = distinctNeighborClusters.keySet().iterator();
		while (itr_cluster.hasNext()) neighborClustersSorted.add(itr_cluster.next());
		Collections.sort(neighborClustersSorted);
		
		for (int targetClusterID : neighborClustersSorted){
			if (thisCluster == targetClusterID) continue;
			
//			System.out.println("move to cluster " + targetClusterID);
			
			double toTargetClusterWeight = 0;
			// iterate through incident edges
			Iterator<Integer> itr_edge = distinctEdgeClusters.get(targetClusterID).iterator();
			while (itr_edge.hasNext()) {
				int edgeID = itr_edge.next();
				double edgeWeight = distinctIncidentEdges.get(edgeID);
				
				// for each node in the edge
				int firstIdx = EINC_head[edgeID];
				int secondIdx = getSecondIdx_EINC(edgeID);
				
//				System.out.print("incident edge " + edgeID + "(" + edgeWeight + "): ");
//				for (int i = firstIdx; i < secondIdx; i++) {
//					int node = EINC_nID[i];
//					int cluster = clusters[node];
//					double nodeWeight = EINC_weight[i];
//					System.out.print(node + "(" + cluster + ";" + nodeWeight + "),");
//				}
//				System.out.println();
				
				double thisWeight = 0;
				double[] fractions = new double[neighborClusterNum];
				for (int i = firstIdx; i < secondIdx; i++) {
					int node = EINC_nID[i];
					double nodeWeight = EINC_weight[i];
					if (node == thisNode) {
						thisWeight = EINC_weight[i];
						continue;
					}
					
					int cluster = clusters[node];
					int clusterIdx = distinctNeighborClusters.get(cluster);
					fractions[clusterIdx] += nodeWeight;
				}
				
//				for (int i = 0; i < fractions.length; i++) {
//					System.out.println("to cluster " + idxToCluster.get(i) + " fraction " + fractions[i] + "/" + edgeWeight + "=" + (fractions[i] / edgeWeight));
//				}
				
				int clusterIdx = distinctNeighborClusters.get(targetClusterID);
				double fraction_u = thisWeight;
				double fraction_cNotu = fractions[clusterIdx];
				double fraction_c = fraction_u + fraction_cNotu;
				if (fraction_cNotu > 0) {
					double existingContribution = adjustWeight_LinearLog(edgeWeight, fraction_u) + adjustWeight_LinearLog(edgeWeight, fraction_cNotu);
					double newContribution = adjustWeight_LinearLog(edgeWeight, fraction_c);
					double deltaContribution = newContribution - existingContribution;
					toTargetClusterWeight += deltaContribution;
					
//					System.out.println("fraction_c " + fraction_c + " fraction_u " + fraction_u + " fraction_cNotu " + fraction_cNotu);
//					System.out.println("existingContribution " + existingContribution + " newContribution " + newContribution + " deltaContribution " + deltaContribution);
				}
			}
			
//			System.out.println("toTargetClusterWeight " + toTargetClusterWeight);
        	
			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			one_minus_gammaExpect_c = 1 - (gamma * expect(vol_c, totalEdgeWeight, ratio));
			one_minus_gammaExpect_u = 1 - (gamma * expect(thisDegree, totalEdgeWeight, ratio));
			one_minus_gammaExpect_cu = 1 - (gamma * expect(vol_cu, totalEdgeWeight, ratio));
			residule = (gamma * (ratio - 1)) -1;
			
			eta = (1.0 / one_minus_gammaExpect_c) + (1 / one_minus_gammaExpect_u) - (1 / one_minus_gammaExpect_cu);
			
			double deltaQ = (toTargetClusterWeight / totalEdgeWeight) + eta + residule;
			double deltaQ_local = deltaQ + deltaQ_exit;
			if (deltaQ_local > deltaQ_best) {
        		deltaQ_best = deltaQ_local;
        		best_targetCluster = targetClusterID;
            }
			
//			System.out.println(targetClusterID + " deltaQ " + deltaQ + " deltaQ_local " + deltaQ_local);
//			System.out.println("");
		}
		
//		System.out.println("deltaQ_best " + deltaQ_best + " best_targetCluster " + best_targetCluster);
		
		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
        	// decide to move
        	increasePerPass += deltaQ_best;
        	cluster_volume[thisCluster] -= thisDegree;
        	cluster_volume[best_targetCluster] += thisDegree;
            clusters[thisNode] = best_targetCluster;
//            System.out.println("new cluster of " + thisNode + " is " + clusters[thisNode]);
            return deltaQ_best + "," + best_targetCluster;
		}
		
//		System.out.println("=====================\n");
		
		return "";
	}
	
	private double expect(double vol_c, double vol_v, double alpha) {
		return (((vol_c / vol_v) - 1) * alpha) + 1;
	}
	
	public double adjustWeight_AON (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			adjustedWeight = edgeWeight;
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Log (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			double invRatio = 1.0 / thisRatio;
			double log = Math.log(invRatio + 1) / Math.log(2);
			adjustedWeight = edgeWeight * (1.0 / log);
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Linear (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			adjustedWeight = edgeWeight * thisRatio;
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_LinearLog (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			double linear = thisRatio;
			double invRatio = 1.0 / thisRatio;
			double log = Math.log(invRatio + 1) / Math.log(2);
			double linear_log = linear * (1.0 / log);
			adjustedWeight = edgeWeight * linear_log;
 		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Sigmoid (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			double invRatio = 1.0 / thisRatio;
			double exp = Math.exp(-1 * (invRatio - 1));
			double sigmod = 2 * ((1.0 / (1+exp)) - 0.5);
			adjustedWeight = edgeWeight * (1 - sigmod);
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Tanh (double thisWeight, double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		if (fraction > 0) {
			double thisRatio = (fraction + thisWeight) / edgeWeight;
			if (thisRatio >= ratio) {
				double invRatio = 1.0 / thisRatio;
				double exp = Math.exp(-1 * 0.5 * (invRatio - 1));
				double sigmod = 2 * (1.0 / (1+exp) - 0.5);
				adjustedWeight = edgeWeight * (1 - sigmod);
			}
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Quadratic (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			double fractionPow = Math.pow(thisRatio, 2);
			adjustedWeight = edgeWeight * fractionPow;
		}
		
		return adjustedWeight;
	}
	
	public double adjustWeight_Cube (double edgeWeight, double fraction) {
		double adjustedWeight = 0;
		
		double thisRatio = fraction / edgeWeight;
		if (thisRatio >= ratio) {
			double fractionPow = Math.pow(thisRatio, 3);
			adjustedWeight = edgeWeight * fractionPow;
		}
		
		return adjustedWeight;
	}
	
	public String move(int i) {
		
		int edgeSize, edgeID, nodeID, clusterID, nodePriority_i, idx, firstIdx_EINC, secondIdx_EINC, secondIdx_INC;
		double deltaQ_local, deltaQ, cutoffSize, potentialCutOffSize, fraction;
		
		double deg_i = node_degrees[i];
		int cur_ClusterID = clusters[i];
		
		// get self loop weight
		double self_loop_weight = 0;
		
		////////////////////////////////////////////////////////////////////////////////////////
		
		// get neighbor
		int distinctNeighborCnt = 0;
		
		// get neighbor clusters
		int distinctClusterCnt = 0;
		
		if (i == n - 1) secondIdx_INC = ADJ_nID_length;
		else secondIdx_INC = ADJ_nHead[i + 1];
		
		if (toHigherOrder) {
			
			for (int j = ADJ_nHead[i]; j < secondIdx_INC; j++) {
				
				nodeID = ADJ_nID[j];
				neighbors[distinctNeighborCnt++] = nodeID;
				
				clusterID = clusters[nodeID];
				if (clusterIdx[clusterID] == -1) {
					nbrClusterIDs[distinctClusterCnt] = clusterID;
					clusterIdx[clusterID] = distinctClusterCnt++;
				}
			}
			
			if (distinctClusterCnt == 1 && nbrClusterIDs[0] == cur_ClusterID) {
				for (int j = 0; j < distinctNeighborCnt; j++) neighbors[j] = -1;
				
				nbrClusterIDs[0] = -1;
				clusterIdx[cur_ClusterID] = -1;
				neighbor_clusters_weights[cur_ClusterID] = 0;
				return "";
			}
			
		} else {
			
			for (int j = ADJ_nHead[i]; j < secondIdx_INC; j++) {
				clusterID = clusters[ADJ_nID[j]];
				
				if (clusterIdx[clusterID] == -1) {
					nbrClusterIDs[distinctClusterCnt] = clusterID;
					clusterIdx[clusterID] = distinctClusterCnt++;
				}
			}
			
			if (distinctClusterCnt == 1 && nbrClusterIDs[0] == cur_ClusterID) {
				nbrClusterIDs[0] = -1;
				clusterIdx[cur_ClusterID] = -1;
				neighbor_clusters_weights[cur_ClusterID] = 0;
				return "";
			}
		}
		
		// for each incident edge
		if (i == n - 1) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[i + 1];
		
		for (int j = INC_head[i]; j < secondIdx_INC; j++) {
			edgeID = INC_eID[j];
			
			firstIdx_EINC = EINC_head[edgeID];
			if (edgeID == m - 1) secondIdx_EINC = EINC_nID_length;
			else secondIdx_EINC = EINC_head[edgeID + 1];
			
			edgeSize = secondIdx_EINC - firstIdx_EINC;
			
			// check if self loop
			if (edgeSize == 1 && EINC_nID[firstIdx_EINC] == i) {
				
				self_loop_weight += edge_weights[edgeID];
				
			} else {
				
				cutoffSize = ratio * (double) edgeSize;
				potentialCutOffSize = edgeSize;
				
				for (int k = firstIdx_EINC; k < secondIdx_EINC; k++) {
					nodeID = EINC_nID[k];
					if (nodeID == i) continue;
					
					idx = clusterIdx[clusters[nodeID]];
					
					if (fractions[idx] < 0) fractions[idx] = 1;
					else fractions[idx]++;
				}
				
				for (int k = 0; k < distinctClusterCnt; k++) {
					fraction = fractions[k];
					fractions[k] = -1;
					
					if (fraction < 0) continue;
					
					if ((fraction + 1) >= cutoffSize) {
						neighbor_clusters_weights[nbrClusterIDs[k]] += edge_weights[edgeID];
					}
					
//					potentialCutOffSize -= fraction;
//					if (potentialCutOffSize < cutoffSize) {
//						for (; k < distinctClusterCnt; k++) {
//							fractions[k] = -1;
//						}
//						break;
//					}
				}
			}
		}
		
		// try removing i from cluster[i]
		double vol = cluster_volume[cur_ClusterID];
		double eta = (1 / (totalEdgeWeight - (gamma * (vol - deg_i)))) + (1 / (totalEdgeWeight - (gamma * deg_i))) - (1 / (totalEdgeWeight - (gamma * vol)));
		double deltaQ_exit = -1 * (neighbor_clusters_weights[cur_ClusterID] / totalEdgeWeight) - ((eta * totalEdgeWeight) - 1);
//		System.out.println(deltaQ_exit);
		
//		// put back
//		vol = cluster_volume[cur_ClusterID] - deg_i;
//		eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * deg_i))) - (1 / (totalEdgeWeight - (gamma * (vol + deg_i))));
//		deltaQ = (neighbor_clusters_weights[cur_ClusterID] / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
//		System.out.println((deltaQ + deltaQ_exit));
		
		int best_ClusterID = cur_ClusterID;
		double deltaQ_best = 0;
		
		// for each neighbor
		if (toHigherOrder) {
			
			clusterIdx[cur_ClusterID] = -1;
			
			nodePriority_i = nodePriority[i];
			
			for (int j = 0; j < distinctNeighborCnt; j++) {
				nodeID = neighbors[j];
				neighbors[j] = -1;
				if (nodePriority[nodeID] < nodePriority_i) continue;
				
				clusterID = clusters[nodeID];
	        	if (clusterIdx[clusterID] == -1) continue;
	        	
	        	vol = cluster_volume[clusterID];
	        	eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * deg_i))) - (1 / (totalEdgeWeight - (gamma * (vol + deg_i))));
	        	deltaQ = (neighbor_clusters_weights[clusterID] / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
	        	
	        	deltaQ_local = deltaQ + deltaQ_exit;
	        	if (deltaQ_local > deltaQ_best) {
	        		deltaQ_best = deltaQ_local;
	        		best_ClusterID = clusterID;
	            }
	        	
	        	clusterIdx[clusterID] = -1;
	        	neighbor_clusters_weights[clusterID] = 0;
			}
			
			// for each cluster
			for (int j = 0; j < distinctClusterCnt; j++) {
				clusterID = nbrClusterIDs[j];
				
				nbrClusterIDs[j] = -1;
				clusterIdx[clusterID] = -1;
				neighbor_clusters_weights[clusterID] = 0;
			}
			
		} else {
			
			// for each cluster
			for (int j = 0; j < distinctClusterCnt; j++) {
				clusterID = nbrClusterIDs[j];
				nbrClusterIDs[j] = -1;
				
				if (clusterID == cur_ClusterID) {
					continue;
				}
	        	
	        	vol = cluster_volume[clusterID];
	        	eta = (1 / (totalEdgeWeight - (gamma * vol))) + (1 / (totalEdgeWeight - (gamma * deg_i))) - (1 / (totalEdgeWeight - (gamma * (vol + deg_i))));
	        	deltaQ = (neighbor_clusters_weights[clusterID] / totalEdgeWeight) + ((eta * totalEdgeWeight) - 1);
	        	
	        	deltaQ_local = deltaQ + deltaQ_exit;
	        	if (deltaQ_local > deltaQ_best) {
	        		deltaQ_best = deltaQ_local;
	        		best_ClusterID = clusterID;
	            }
	        	
//	        	System.out.println("toTargetClusterWeight " + neighbor_clusters_weights[clusterID]);
//	        	System.out.println("cluster " + clusterID + " deltaQ " + deltaQ + " deltaQ_local " + deltaQ_local);
	        	
	        	clusterIdx[clusterID] = -1;
	        	neighbor_clusters_weights[clusterID] = 0;
			}
			
			clusterIdx[cur_ClusterID] = -1;
			neighbor_clusters_weights[cur_ClusterID] = 0;
		}
        
//		System.out.println("deltaQ_best " + deltaQ_best + " best_ClusterID " + best_ClusterID);
		
        if (deltaQ_best > 0 && best_ClusterID != cur_ClusterID) {
        	// decide to move
        	increasePerPass += deltaQ_best;
        	cluster_volume[cur_ClusterID] -= deg_i;
        	cluster_volume[best_ClusterID] += deg_i;
            clusters[i] = best_ClusterID;
            
            return deltaQ_best + "," + best_ClusterID;
        }
        
        return "";
	}

	public void rebuildGraph(boolean firstReConstruct) {
		
		distinctIDs = new ArrayList<Integer>();
		
		new_n = 0; // the number of nodes in new hypergraph
		int clusterID;
		for (int i = 0; i < n; i++) {
			clusterID = clusters[i];
			
			if (clusterIdx[clusterID] == -1) {
				clusterIdx[clusterID] = new_n++;
				distinctIDs.add(clusterID);
			}
		}
		
		if (firstReConstruct) {
			nodesInCluster = new ArrayList[new_n];
			newINC_eIDs = new ArrayList[new_n];
			
			// update self_loop_weights
			self_loop_weights = new double[new_n];
		}
		
		for (int i = 0; i < new_n; i++) {
			nodesInCluster[i] = new ArrayList<Integer>();
			newINC_eIDs[i] = new ArrayList<Integer>();
			
			// update node_degrees
			node_degrees[i] = 0;
		}
		
		int newID;
		for (int i = 0; i < n; i++) {
			newID = clusterIdx[clusters[i]];
			clusters[i] = newID;
			nodesInCluster[newID].add(i);
		}
		
//		interestedNode = clusters[interestedNode];
		
		for (int c : distinctIDs) clusterIdx[c] = -1;
		
		// update global_cluster
		if (save) {
			int prevID;
			for (int i = 0; i < global_cluster.length; i++) {
				prevID = global_cluster[i];
				newID = clusters[prevID];
				global_cluster[i] = newID;
			}
		}
		
		int edgeCnt = 0;
		int edge_head = 0;
		int cluster_v, second_idx_u, second_idx_v;
		double weight;
		
		// store the nodes and their weights in each hyperedge
		newEINC_heads = new ArrayList<Integer>();
		newEINC_nIDs = new ArrayList<Integer>();
		newEINC_weight = new ArrayList<Double>();
		
		// for each cluster
		for (int cluster_u = 0; cluster_u < new_n; cluster_u++) {
			
			// for each node in the cluster
			int edgeID;
			for (int node : nodesInCluster[cluster_u]) {
				
				// for each incident hyperedge
				if (node == n - 1) second_idx_u = INC_eID_length;
				else second_idx_u = INC_head[node + 1];
				
				for (int edgeIdx = INC_head[node]; edgeIdx < second_idx_u; edgeIdx++) {
					
					edgeID = INC_eID[edgeIdx];
					if (edgeInfos[edgeID] != -1) continue;
					edgeInfos[edgeID] = 0;
					
					// store the new nodes' ID
					distinctIDs = new ArrayList<Integer>();
					
					// for each node old node in this hyperedge
					if (edgeID == m - 1) second_idx_v = EINC_nID_length;
					else second_idx_v = EINC_head[edgeID + 1];
					
					for (int i = EINC_head[edgeID]; i < second_idx_v; i++) {
						// match old node to new node
						cluster_v = clusters[EINC_nID[i]];
						
						if (clusterIdx[cluster_v] == -1) {
							distinctIDs.add(cluster_v);
							clusterIdx[cluster_v] = 0;
						}
						
						cluster_weights[cluster_v] += EINC_weight[i];
					}
					
					if (distinctIDs.size() > 1) {
						
						newEINC_heads.add(edge_head);
						
						for (int i = 0; i < distinctIDs.size(); i++) {
							clusterID = distinctIDs.get(i);
							clusterIdx[clusterID] = -1;
							
							weight = cluster_weights[clusterID];
							cluster_weights[clusterID] = 0;
							
							node_degrees[clusterID] += weight;
							
							// new node "clusterID" has incident hyperedge with ID "head"
							newINC_eIDs[clusterID].add(edgeCnt);
							
							newEINC_nIDs.add(clusterID);
							newEINC_weight.add(weight);
							edge_head++;
						}
						
						// new node set in hyperedge with ID "edgeCnt"
						edgeCnt++;
						
					} else {
						clusterID = distinctIDs.get(0);
						clusterIdx[clusterID] = -1;
						
						weight = cluster_weights[clusterID];
						cluster_weights[clusterID] = 0;
						
						self_loop_weights[clusterID] += weight;
					}
				}
			}
		}
		
		for (int j = 0; j < m; j++) {
			edgeInfos[j] = -1;
		}
		
		double self_loop_weight, vol_c, gamma_c;
		Q_curr = 0;
		
		// iterate self-loop edge
		for (clusterID = 0; clusterID < new_n; clusterID++) {
			if (self_loop_weights[clusterID] > 0) {
				
				self_loop_weight = self_loop_weights[clusterID];
				self_loop_weights[clusterID] = 0;
				
				if (save) {
					// calculate modularity
					vol_c = node_degrees[clusterID];
					gamma_c = gamma * vol_c;
					Q_curr += (((self_loop_weight + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c)));
//					System.out.println("cluster " + clusterID + " vol_c " + vol_c + " gamma_c " + gamma_c + " Q " + (((self_loop_weight + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c))));
				}
				
				node_degrees[clusterID] += self_loop_weight;
				
				// new node "clusterID" has incident hyperedge with ID "head"
				newINC_eIDs[clusterID].add(edgeCnt);
				edgeCnt++;
				
				// self-loop hyperedge
				newEINC_heads.add(edge_head);
				newEINC_nIDs.add(clusterID);
				newEINC_weight.add(self_loop_weight);
				edge_head++;
			}
		}
		
		// update n
		n = new_n;
		
		// update m
		m = newEINC_heads.size();
		
		int totalCardinality = edge_head;
		
		// update EINC_nID_length
		EINC_nID_length = totalCardinality;
		
		// update edge_weights
		for (int i = 0; i < m; i++) {
			edge_weights[i] = 0;
		}
		
		// update EINC_head
		int idx = 0;
		itr_int = newEINC_heads.iterator();
		while (itr_int.hasNext()) {
			EINC_head[idx++] = itr_int.next();
		}
		
		// update EINC_nID, EINC_weight
		idx = 0;
		itr_int = newEINC_nIDs.iterator();
		itr_double = newEINC_weight.iterator();
		while (itr_int.hasNext()) {
			EINC_nID[idx] = itr_int.next();
			EINC_weight[idx] = itr_double.next();
			idx++;
		}
		
		// for each edge
		if (m > 1) {
			int mMinusOne = m - 1;
			int curEdgeID = 0;
			int nextHead = EINC_head[curEdgeID + 1];
			idx = 0;
			itr_double = newEINC_weight.iterator();
			while (itr_double.hasNext()) {
				weight = 0;
				while (idx++ < nextHead) {
					weight += itr_double.next();
				}
				idx--;
				
				edge_weights[curEdgeID] = weight;
				curEdgeID++;
				
				if (curEdgeID < mMinusOne) nextHead = EINC_head[curEdgeID + 1];
				else {
					nextHead = newEINC_weight.size();
					break;
				}
			}
			weight = 0;
			while (idx++ < nextHead) weight += itr_double.next();
			edge_weights[curEdgeID] = weight;
			
		} else {
			edge_weights[0] = newEINC_weight.iterator().next();
		}
		
		// update INC_head, INC_eID, INC_weight
        INC_eID_length = totalCardinality;
        
        int node_head = 0;
        for (int nodeID = 0; nodeID < n; nodeID ++) {
			INC_head[nodeID] = node_head;
			for (int edgeID : newINC_eIDs[nodeID]) {
				INC_eID[node_head] = edgeID;
				node_head++;
			}
		}
        
//        // check correctness
// 		double totalNode_weights = 0;
// 		for (int i = 0; i < n; i++) totalNode_weights += node_degrees[i];
// 		double totalEdge_weights = 0;
// 		for (int i = 0; i < EINC_nID_length; i++) totalEdge_weights += EINC_weight[i];
// 		double totalWeights = 0;
// 		for (int i = 0; i < m; i++) totalWeights += edge_weights[i];
// 		System.out.println("# new nodes " + n + " # of new hyperedges " + m + " Q_curr " + Q_curr + 
// 				" total weight" + totalNode_weights + "  " + totalEdge_weights + " " + totalWeights);
	}
	
	public String louvain(String moveStrategy) throws Exception {
		
		Random random = new Random();
		
		int count_aggregations = 0;
		this.moveStrategy = moveStrategy;
		
		System.out.println(dataset + " " + this.lambda + " " + gamma);
		System.out.println(trial + " " + toHigherOrder + " " + ordering + " " + moveStrategy + " " + ratio + " " + save);
		
//		printNodeInfo();
		
		// variables for ordering
		int j, temp_node;
		
		// variables for rebuildGraph
		boolean firstReConstruct = true;
		
		// variables for updateGamma
		double timeUpdateGamma;
		
		// variables for move
		boolean increase;
		
		Q_prev = -1;
		
		double startTime, endTime;
		double extractClusterTime = 0;
		Hypergraph.garbbageCollector.gc();
		long startMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
		do {
			System.out.println("itr " + count_aggregations);
			count_aggregations++;
			
			////////////////////////////////////////////////////////////////////////
			
			startTime = System.currentTimeMillis();
			
			switch (ordering) {
				
			case "randomOrder":
			{
				// random the order
				for (int i = 0; i < n; i++) {
					j = random.nextInt(n);
					temp_node = nodeInOrder[i];
					nodeInOrder[i] = nodeInOrder[j];
					nodeInOrder[j] = temp_node;
				}
			}
				break;
			}
			
			switch (moveStrategy) {
			case "move":
			{
				int count_move = 0;
				
				increase_total = 0;
				increase = true;
				while (increase) {
					count_move++;
					
					increase = false;
					increasePerPass = 0;
					
					for (int idx = 0; idx < nodeInOrder.length; idx++) {
//						toPrint = "";
						move_split_adjust(nodeInOrder[idx]);
//						if (nodeInOrder[idx] == interestedNode) {
//							System.out.println(toPrint);
//						}
					}
					
//					System.out.println("pass " + count_move + " increasePerPass " + increasePerPass);
					increase_total += increasePerPass;
					if (increasePerPass > eps) increase = true;
					if (count_move > n_aggregations) break;
				}
			}
			break;
			}
			
//			System.out.println("increase_total " + increase_total);
			if (increase_total <= eps) break;
			
			rebuildGraph(firstReConstruct);
			firstReConstruct = false;
			
			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);
			
			/////////////////////////////////////////////////////////////////////
			
			timeUpdateGamma = updateGamma();
			if (timeUpdateGamma == -999) break;
			extractClusterTime += (timeUpdateGamma * Constant.RUNNING_TIME_UNIT);
			
			/////////////////////////////////////////////////////////////////////
			
			if (count_aggregations > n_aggregations) break;
			
			nodeInOrder = new int[n];
			for (int i = 0; i < n; i++) {
				nodeInOrder[i] = i;
				clusters[i] = i;
				cluster_volume[i] = node_degrees[i];
			}
			
//			printNodeInfo();
			
		} while (true);
		
		long endMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
		long memoryUse = endMem - startMem;
				
//		System.out.println("running time " + String.format("%.8f", (extractClusterTime / Constant.RUNNING_TIME_UNIT)));
			
		if (save) {
			saveClusters();
			saveModularity();
		}
		
		return extractClusterTime + "," + memoryUse;
	}
	
	public double updateGamma() {
		
		// get cardinality distribution
		try {
			FileWriter fwCount = new FileWriter(FilePath_Mon.filePathPre + "/stat/cardinality/cardinality_trial_" + trial + ".txt");
			int cardinality;
			int mMinusOne = m - 1;
			for (int edgeID = 0; edgeID < m; edgeID++) {
				
				if (edgeID == mMinusOne) cardinality = (EINC_nID_length - EINC_head[edgeID]);
				else cardinality = (EINC_head[edgeID + 1] - EINC_head[edgeID]);
				
				if (cardinality >= 2) fwCount.write(cardinality + "\n");
			}
			fwCount.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -999;
		}
		
		double extractClusterTime = -999;
		// ref: https://www.edureka.co/community/358/how-to-execute-a-python-file-with-few-arguments-in-java
		try {
			process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/getLambda.py " 
					+ "/data/" + dataset + "/stat/cardinality/cardinality_trial_" + trial + ".txt" + " " + 2);
		} catch (Exception e) {
			System.out.println("Exception Raised" + e.toString());
		}
		InputStream stdout = process.getInputStream();
		BufferedReader reader = new BufferedReader(new InputStreamReader(stdout, StandardCharsets.UTF_8));
		String line;
		try {
			
			// read the first line
			line = reader.readLine();
			if (line != null){
				extractClusterTime = Double.parseDouble(line.split(",")[1]);
			} else {
				return -999;
			}
			
			// read the second line
			line = reader.readLine();
			if (line != null && !line.equals("nan")){
				this.lambda = Double.parseDouble(line);
				this.gamma = Math.exp(-1 * lambda);
				
				if (gamma < eps) return -999;
				
//				System.out.println(lambda + " " + gamma);
			} else {
				return -999;
			}
		} catch (IOException e) {
			System.out.println("Exception in reading output" + e.toString());
		}
		
		return extractClusterTime;
	}
	
	public boolean getGamma(int edgeSize) {
		
		// get cardinality distribution
		try {
			FileWriter fwCount = new FileWriter(FilePath_Mon.filePathPre + "/stat/cardinality/cardinality_trial_" + trial + ".txt");
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
		try {
			process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/getLambda.py " 
					+ "/data/" + dataset + "/stat/cardinality/cardinality_trial_" + trial + ".txt" + " " + 2);
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
				this.lambda = Double.parseDouble(line);
				this.gamma = Math.exp(-1 * lambda);
				
				if (gamma < eps) return false;
//				System.out.println(lambda + " " + gamma);
			} else return false;
		} catch (IOException e) {
			System.out.println("Exception in reading output" + e.toString());
		}
				
		return true;
	}
	
	public double getHModularity() {
		
		int clusterNum = n;
		
		List<Integer>[] nodesInCluster = new ArrayList[n];
		for (int i = 0; i < clusterNum; i++) {
			nodesInCluster[i] = new ArrayList<Integer>();
		}
		
		for (int i = 0; i < global_cluster.length; i++) {
			nodesInCluster[global_cluster[i]].add(i);
		}
		
		// for each cluster
		double[] cluster_vol = new double[clusterNum];
		for (int cluster = 0; cluster < clusterNum; cluster++) {
			
			// for each node
			int secondIdx_INC, firstIdx_INC;
			for (int node : nodesInCluster[cluster]) {
				// for each incident edge
				firstIdx_INC = Hypergraph.INC_head[node];
				if (node == Hypergraph.getVertexSize() - 1) secondIdx_INC = Hypergraph.INC_eID.length;
				else secondIdx_INC = Hypergraph.INC_head[node + 1];
				
				double weight = 0;
				for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
					weight += Hypergraph.INC_weight[i];
				}
				cluster_vol[cluster] += weight;
			}
		}
		
//		double totalClusterVol = 0;
//		// for each cluster
//		for (int cluster = 0; cluster < clusterNum; cluster++) {
//			totalClusterVol += cluster_vol[cluster];
//		}
//		System.out.println("totalClusterVol " + totalClusterVol);
		
		// for each edge
		double[] self_loop_weights = new double[n];
		HashSet<Integer> distinctClusters;
		int secondIdx_EINC, firstIdx_EINC;
		for (int edgeID = 0; edgeID < Hypergraph.getEdgeSize(); edgeID++) {
			
			distinctClusters = new HashSet<Integer>();
			
			// for each node in the edge
			firstIdx_EINC = Hypergraph.EINC_head[edgeID];
			if (edgeID == Hypergraph.getEdgeSize() - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
			
//			System.out.print("edge " + edgeID + ": ");
			
			double weight = 0;
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				int cluster = global_cluster[Hypergraph.EINC_nID[i]];
				distinctClusters.add(cluster);
				
				weight += Hypergraph.EINC_weight[i];
				
//				System.out.print(Hypergraph.EINC_nID[i] + "(" + global_cluster[Hypergraph.EINC_nID[i]] + ")");
			}
//			System.out.println("weight " + weight);
			 
			// check is inner edge
			if (distinctClusters.size() <= 1) {
				int cluster = distinctClusters.iterator().next();
				
//				System.out.println("add to cluster " + cluster);
				
				self_loop_weights[cluster] += weight;
			}
		}
		
//		double totalSelf_loop_weight = 0;
//		// for each cluster
//		for (int cluster = 0; cluster < clusterNum; cluster++) {
//			totalSelf_loop_weight += self_loop_weights[cluster];
//		}
//		System.out.println("totalSelf_loop_weight" + totalSelf_loop_weight);
		
		double Q = 0;
		
		for (int cluster = 0; cluster < clusterNum; cluster++) {
			
			// calculate modularity
			double vol_c = node_degrees[cluster];
			double gamma_c = initGamma * vol_c;
			double self_loop_weight = self_loop_weights[cluster];
			
			Q += (((self_loop_weight + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c)));
//			System.out.println("cluster " + cluster + " vol_c " + vol_c + " gamma_c " + gamma_c + " Q " + (((self_loop_weight + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c))));
		}
		
		System.out.println("getModularity " + Q);
		
		return Q;
	}
	
	public void saveClusters() {
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/node_cluster_Louvain_" + this.moveStrategy + "_ordered_"
					+ "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + ratio + "_trial_" + (this.trial) + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/node_cluster_Louvain_" + this.moveStrategy + "_ordered_"
					+ "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + ratio + "_connect_trial_" + (this.trial) + ".txt";
		}
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			
			nodesInCluster = new ArrayList[n];
			for (int i = 0; i < n; i++) {
				nodesInCluster[i] = new ArrayList<Integer>();
			}
			
			for (int i = 0; i < global_cluster.length; i++) {
				nodesInCluster[global_cluster[i]].add(i);
			}
			
			for (int i = 0; i < nodesInCluster.length; i++) {
				if (nodesInCluster[i].size() > 0) {
					fwCount.write(nodesInCluster[i].stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
				}
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}

	public void saveModularity() {
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/modularity_Louvain_" + this.moveStrategy + "_ordered_"
					+ "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + this.ratio + "_trial_" + this.trial + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/modularity_Louvain_" + this.moveStrategy + "_ordered_"
					+ "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + this.ratio + "_connect_trial_" + this.trial + ".txt";
		}
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			
			double modularity = 0;
			
			fwCount.write("modularity\t" + modularity + "\n");
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public int getSecondIdx_EINC(int edgeID) {
		if (edgeID == m - 1) return EINC_nID_length;
		else return EINC_head[edgeID + 1];
	}
	
	public void printNodeInfo() {
		// cluster
		List<Integer>[] nodesInCluster = new ArrayList[n];
		for (int i = 0; i < n; i++) {
			nodesInCluster[i] = new ArrayList<Integer>();
		}
		
		for (int i = 0; i < global_cluster.length; i++) {
			nodesInCluster[global_cluster[i]].add(i);
		}
		
		int[] clusterSize = new int[n];
		for (int i = 0; i < n; i++) {
			clusterSize[i] = nodesInCluster[i].size();
		}
		
		// node
		// for each node
		for (int nodeID = 1200; nodeID < 1220; nodeID++) {
			int cluster = global_cluster[nodeID];
			
			// for each incident edge
			int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[nodeID];
			if (nodeID == Hypergraph.getVertexSize() - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
			
			int degree = secondIdx_INC - firstIdx_INC;
			
			double avgCardinality = 0;
			for (int k = firstIdx_INC; k < secondIdx_INC; k++) {
				int edgeID = Hypergraph.INC_eID[k];
				
				// for each node in the edge
				int secondIdx_EINC, firstIdx_EINC = Hypergraph.EINC_head[edgeID];
				if (edgeID == Hypergraph.getEdgeSize() - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
				
				double cardinality = secondIdx_EINC - firstIdx_EINC;
				avgCardinality += cardinality;
			}
			avgCardinality /= degree;
			
			System.out.println("node " + nodeID + "(" + cluster + ";" + clusterSize[cluster] + ") degree " + degree + " avgincidentCardinality " + avgCardinality);
		}
		
		System.out.println("===================\n");
	}
	
	public static void main(String arg[]) throws Exception {
		Hypergraph.loadGraph();
		
		boolean toHigherOrder = false;
		String ordering = "randomOrder";
		String moveStrategy = "move";
		double ratio = 0.3;
		boolean save = true;
		
//		LouvainHelper LouvainHelper = new LouvainHelper();
//		System.out.println("ratio " + LouvainHelper.findRatioParallel(toHigherOrder, ordering, moveStrategy, 14));
		
		double avgRunningTime = 0;
		int trials = 1;
		for (int trial = 0; trial < trials; trial++) {
			System.out.println("trial " + trial);
			PIC l = new PIC(trial, toHigherOrder, ordering, ratio, save);
			l.initPart1();
			l.initPart2();
			avgRunningTime += Double.parseDouble(l.louvain(moveStrategy).split(",")[0]);
		}
		
		avgRunningTime /= trials;
		
		System.out.println("average running time " + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)));
	}
}
