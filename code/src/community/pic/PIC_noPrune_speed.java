package community.pic;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;

// calculate the running time of PIC algorithm with no optimization technique
// Note: on the hypergraphs with very large sizes, use
// 1) call constructADJ_largeGraph() in initPart1() and pic() to construct the adjacent list
// 2) call move_split_largeGraph() in pic() to move a node 
public class PIC_noPrune_speed {

	int trial;
	double ratio;
	String ordering;

	double incPerIteration;
	int n_aggregations = 10;

	double lambda;
	double gamma = -1;
	double initGamma;
	double eps = Constant.EPSILON;

	int n, m;
	int endN, endM;
	double[] node_degrees;

	double totalEdgeWeight;
	double[] INC_weight;
	int[] INC_eID;
	int[] INC_head;
	int INC_eID_length;

	int[] EINC_nID;
	double[] EINC_weight;
	int[] EINC_head;
	int EINC_nID_length;
	double[] edge_weights;

	int curTimeStamp;
	int[] lastChangeTimes;

	int[] ADJ_nID;
	int[] ADJ_head;
	int ADJ_nID_length;
	
	byte[] largeSplit;
	int[][] largeADJ_nID;
	int[][] largeADJ_head;
	int[] largeADJ_nID_length;
	int[] largeEndNode;

	// variables for move
	int iteration;
	boolean[] visited;
	double[] fractions;
	int[] incident_clusteridxs;
	double[] incident_weights;
	int[] possibleNeighborClusters;
	int[] impossibleNeighborClusters;

	List<Integer>[] nodesInCluster;
	int[] global_cluster;
	int[] clusters;
	double[] cluster_volume;

	int[] nodeInOrder;
	int[] nodePriority;

	// variables for rebuildGraph
	int[] newEINC_head;
	int[] newEINC_nID;
	double[] newEINC_weight;
	List<Integer>[] newINC_eIDs;
	List<Double>[] newINC_weights;
	int[] clusterIdxs;
	double[] cluster_weights;
	double[] self_loop_weights;
	int new_n;

	boolean usePruning = false;
	boolean save;
	String dataset;

	public PIC_noPrune_speed(int trial, String ordering, double ratio,
			boolean save) throws IOException {
		this.trial = trial;
		this.ordering = ordering;
		this.ratio = ratio;
		this.save = save;
		this.dataset = Hypergraph.dataset;
	}

	public void initPart1() {
		
		// initialize the program by 
		// 1) copying the hypergraph in a form of incident matrix and,
		// 2) constructing the adjacent list
		
		n = Hypergraph.getVertexSize();
		endN = n - 1;
		m = Hypergraph.getEdgeSize();
		endM = m - 1;

		EINC_head = new int[m];
		for (int i = 0; i < m; i++) {
			EINC_head[i] = Hypergraph.EINC_head[i];
		}

		INC_head = new int[n];
		for (int i = 0; i < n; i++) {
			INC_head[i] = Hypergraph.INC_head[i];
		}

		int totalCardinality = Hypergraph.EINC_nID.length;
		EINC_nID = new int[totalCardinality];
		EINC_weight = new double[totalCardinality];
		INC_eID = new int[totalCardinality];
		INC_weight = new double[totalCardinality];
		for (int i = 0; i < totalCardinality; i++) {
			EINC_nID[i] = Hypergraph.EINC_nID[i];
			EINC_weight[i] = Hypergraph.EINC_weight[i];
			INC_eID[i] = Hypergraph.INC_eID[i];
			INC_weight[i] = Hypergraph.INC_weight[i];
		}

		EINC_nID_length = totalCardinality;
		INC_eID_length = totalCardinality;

		// variables for move
		visited = new boolean[n];
		
		// construct the adjacent list 
		constructADJ(true);		// uncomment this for normal-sized hypergraph
//		constructADJ_largeGraph(true);	// uncomment this for very large hypergraph
	}

	public double initPart2() {
		
		// initialize the program by 
		// 1) calculating the gamma of the initial hypergraph
		// 2) calculating the node degrees and edge weights
		// 3) initializing the variables for later computations
		
		// calculate gamma
		double timeUpdateGamma = calculateGamma();
		double runningTime = (timeUpdateGamma * Constant.RUNNING_TIME_UNIT);

		initGamma = gamma;
		
		double time = System.currentTimeMillis();
		
		edge_weights = new double[m];
		totalEdgeWeight = Constant.INITIAL_NODEONEDGE_WEIGHT * m;
		
		for (int edge = 0; edge < m; edge++) {
			edge_weights[edge] = Constant.INITIAL_NODEONEDGE_WEIGHT;
		}

		node_degrees = new double[n];

		// for each node
		double weight;
		int second_idx;
		int node = 0;
		for (node = 0; node < n - 1; node++) {
			// for each incident hyperedge
			second_idx = INC_head[node + 1];

			weight = 0;
			for (int i = INC_head[node]; i < second_idx; i++) weight += INC_weight[i];
			node_degrees[node] = weight;
		}
		second_idx = INC_weight.length;
		weight = 0;
		for (int i = INC_head[node]; i < second_idx; i++) weight += INC_weight[i];
		node_degrees[node] = weight;

		if (save) {
			global_cluster = new int[n];
			for (int i = 0; i < n; i++) {
				global_cluster[i] = i;
			}
		}

		///////////////////////////////////////////////////////////////////////////

		// variables for reconstruct
		newEINC_head = new int[m];
		newEINC_nID = new int[EINC_nID.length];
		newEINC_weight = new double[EINC_weight.length];
		cluster_weights = new double[n];
		clusterIdxs = new int[n];
		for (int i = 0; i < n; i++) clusterIdxs[i] = -1;

		////////////////////////////////////////////////////////////////////////////
		
		nodeInOrder = new int[n];
		clusters = new int[n];
		cluster_volume = new double[n];
		for (int i = 0; i < n; i++) {
			nodeInOrder[i] = i;
			clusters[i] = i;
			cluster_volume[i] += node_degrees[i];
		}
		
		runningTime += (System.currentTimeMillis() - time);
		
		return runningTime;
	}

	public void move_split(int thisNode) {

		// current cluster of the node
		int thisCluster = clusters[thisNode];

		// adjacent list of the node
		int secondIdx_ADJ;
		if (thisNode == endN) secondIdx_ADJ = ADJ_nID_length;
		else secondIdx_ADJ = ADJ_head[thisNode + 1];

		int cluster;
		int possibleNeighborClusterNum = 0;
		int distinctNeighborClusterNum = 0;
		
		// find the neighbor clusters
		// add cluster of thisNode
		possibleNeighborClusters[possibleNeighborClusterNum++] = thisCluster;
		clusterIdxs[thisCluster] = distinctNeighborClusterNum++;
		// add cluster of adjacent neighbors
		for (int i = ADJ_head[thisNode]; i < secondIdx_ADJ; i++) {
			cluster = clusters[ADJ_nID[i]];
			if (clusterIdxs[cluster] == -1) {
				possibleNeighborClusters[possibleNeighborClusterNum++] = cluster;
				clusterIdxs[cluster] = distinctNeighborClusterNum++;
			}
		}
		
		// if all neighbor clusters are the same as the current cluster
		if (distinctNeighborClusterNum == 1) {
			clusterIdxs[thisCluster] = -1;
			return;
		}

		int incident_clusterNum = 0;

		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == endN) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];

		// calculate the delta support
		int edge, node, clusterIdx;
		double thisWeight, edgeWeight;
		double ratio_u, ratio_cNotu, ratio_c;
		double contribute_u, contribute_cNotu, contribute_c;
		for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
			edge = INC_eID[i];
			thisWeight = INC_weight[i];
			edgeWeight = edge_weights[edge];

			// for each node in the edge
			int secondIdx, firstIdx = EINC_head[edge];
			if (edge == endM) secondIdx = EINC_nID_length;
			else secondIdx = EINC_head[edge + 1];

			incident_clusterNum = 0;
			for (int j = firstIdx; j < secondIdx; j++) {
				node = EINC_nID[j];
				if (node == thisNode) continue;

				clusterIdx = clusterIdxs[clusters[node]];
				fractions[clusterIdx] += EINC_weight[j];
				if (!visited[clusterIdx]) {
					incident_clusteridxs[incident_clusterNum++] = clusterIdx;
					visited[clusterIdx] = true;
				}
			}

			ratio_u = thisWeight / edgeWeight;
			if (ratio_u >= ratio) contribute_u = adjustWeight_LinearLog(edgeWeight, ratio_u);
			else contribute_u = 0;

			for (int j = 0; j < incident_clusterNum; j++) {

				clusterIdx = incident_clusteridxs[j];
				visited[clusterIdx] = false;

				ratio_cNotu = fractions[clusterIdx] / edgeWeight;
				fractions[clusterIdx] = 0;

				if (ratio_cNotu >= ratio) {
					contribute_cNotu = adjustWeight_LinearLog(edgeWeight, ratio_cNotu);
				} else
					contribute_cNotu = 0;

				ratio_c = ratio_u + ratio_cNotu;
				if (ratio_c >= ratio) {
					contribute_c = adjustWeight_LinearLog(edgeWeight, ratio_c);
				} else {
					contribute_c = 0;
				}

				incident_weights[clusterIdx] += (contribute_c - contribute_u - contribute_cNotu);
			}
		}

		// try removing the node from current cluster
		double vol_cu, delta_eta_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		
		double eta_cu;
		double CC = 1.0 - gamma;
		double oneMinusGamma = 1.0 - gamma;
		double eta_cNou = ratio * (1.0 - (vol_cNou / totalEdgeWeight));
		double delta_eta_cNou =  (1.0 - eta_cNou) * (1.0 - eta_cNou) * CC / (oneMinusGamma + (gamma * eta_cNou));
		double eta_u = ratio * (1.0 - (thisDegree / totalEdgeWeight));
		double delta_eta_u = (1.0 - eta_u) * (1.0 - eta_u) * CC / (oneMinusGamma + (gamma * eta_u));		
		double eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
		double delta_eta_c = (1.0 - eta_c) * (1.0 - eta_c) * CC / (oneMinusGamma + (gamma * eta_c));
		

		clusterIdx = clusterIdxs[thisCluster];
		clusterIdxs[thisCluster] = -1;
		double deltaQ_exit = (incident_weights[clusterIdx] / totalEdgeWeight) + delta_eta_cNou + delta_eta_u - delta_eta_c;
		incident_weights[clusterIdx] = 0;
		deltaQ_exit = -1 * deltaQ_exit;

		// try merging the node to each of its distinct neighbor clusters
		int targetClusterID;
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;

		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			targetClusterID = possibleNeighborClusters[i];
			if (thisCluster == targetClusterID) continue;

			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			
			eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
			delta_eta_c = (1.0 - eta_c) * (1.0 - eta_c) * CC / (oneMinusGamma + (gamma * eta_c));
			eta_u = ratio * (1.0 - (thisDegree / totalEdgeWeight));
			delta_eta_u = (1.0 - eta_u) * (1.0 - eta_u) * CC / (oneMinusGamma + (gamma * eta_u));
			eta_cu = ratio * (1.0 - (vol_cu / totalEdgeWeight));
			delta_eta_cu = (1.0 - eta_cu) * (1.0 - eta_cu) * CC / (oneMinusGamma + (gamma * eta_cu));

			clusterIdx = clusterIdxs[targetClusterID];
			clusterIdxs[targetClusterID] = -1;
			double deltaQ = (incident_weights[clusterIdx] / totalEdgeWeight) + delta_eta_c + delta_eta_u - delta_eta_cu;
			incident_weights[clusterIdx] = 0;
			double deltaQ_local = deltaQ + deltaQ_exit;
			
			// update the best modularity gain and cluster
			if (deltaQ_local > deltaQ_best) {
				deltaQ_best = deltaQ_local;
				best_targetCluster = targetClusterID;
			}
		}

		// if there is positive modularity gain, change the clustering
		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
			// decide to move
			incPerIteration += deltaQ_best;
			cluster_volume[thisCluster] -= thisDegree;
			cluster_volume[best_targetCluster] += thisDegree;
			clusters[thisNode] = best_targetCluster;

			return;
		}

		return;
	}

	public void move_split_largeGraph(int thisNode) {

		// current cluster of the node
		int thisCluster = clusters[thisNode];

		int split = largeSplit[thisNode];
		int endNode = largeEndNode[split];
		
		// adjacent list of the node
		int secondIdx_ADJ;
		if (thisNode == endNode) secondIdx_ADJ = largeADJ_nID_length[split];
		else secondIdx_ADJ = largeADJ_head[split][thisNode + 1];

		int cluster;
		int possibleNeighborClusterNum = 0;
		int distinctNeighborClusterNum = 0;
		
		// find the neighbor clusters
		// add cluster of thisNode
		possibleNeighborClusters[possibleNeighborClusterNum++] = thisCluster;
		clusterIdxs[thisCluster] = distinctNeighborClusterNum++;
		visited[thisCluster] = true;
		// add cluster of adjacent neighbors
		for (int i = largeADJ_head[split][thisNode]; i < secondIdx_ADJ; i++) {
			cluster = clusters[largeADJ_nID[split][i]];
			if (!visited[cluster]) {
				possibleNeighborClusters[possibleNeighborClusterNum++] = cluster;
				clusterIdxs[cluster] = distinctNeighborClusterNum++;
				visited[cluster] = true;
			}
		}
		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			visited[possibleNeighborClusters[i]] = false;
		}
		
		// if all neighbor clusters are the same as the current cluster
		if (distinctNeighborClusterNum == 1) {
			clusterIdxs[thisCluster] = -1;
			return;
		}

		int incident_clusterNum = 0;

		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == endN) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];

		// calculate the delta support
		int edge, node, clusterIdx;
		double thisWeight, edgeWeight;
		double ratio_u, ratio_cNotu, ratio_c;
		double contribute_u, contribute_cNotu, contribute_c;
		for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
			edge = INC_eID[i];
			thisWeight = INC_weight[i];
			edgeWeight = edge_weights[edge];

			// for each node in the edge
			int secondIdx, firstIdx = EINC_head[edge];
			if (edge == endM) secondIdx = EINC_nID_length;
			else secondIdx = EINC_head[edge + 1];

			incident_clusterNum = 0;
			for (int j = firstIdx; j < secondIdx; j++) {
				node = EINC_nID[j];
				if (node == thisNode) continue;

				clusterIdx = clusterIdxs[clusters[node]];
				fractions[clusterIdx] += EINC_weight[j];
				if (!visited[clusterIdx]) {
					incident_clusteridxs[incident_clusterNum++] = clusterIdx;
					visited[clusterIdx] = true;
				}
			}

			ratio_u = thisWeight / edgeWeight;
			if (ratio_u >= ratio) contribute_u = adjustWeight_LinearLog(edgeWeight, ratio_u);
			else contribute_u = 0;

			for (int j = 0; j < incident_clusterNum; j++) {

				clusterIdx = incident_clusteridxs[j];
				visited[clusterIdx] = false;

				ratio_cNotu = fractions[clusterIdx] / edgeWeight;
				fractions[clusterIdx] = 0;

				if (ratio_cNotu >= ratio) {
					contribute_cNotu = adjustWeight_LinearLog(edgeWeight, ratio_cNotu);
				} else
					contribute_cNotu = 0;

				ratio_c = ratio_u + ratio_cNotu;
				if (ratio_c >= ratio) {
					contribute_c = adjustWeight_LinearLog(edgeWeight, ratio_c);
				} else {
					contribute_c = 0;
				}

				incident_weights[clusterIdx] += (contribute_c - contribute_u - contribute_cNotu);
			}
		}

		// try removing the node from current cluster
		double vol_cu, delta_eta_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		
		double eta_cu;
		double CC = 1.0 - gamma;
		double oneMinusGamma = 1.0 - gamma;
		double eta_cNou = ratio * (1.0 - (vol_cNou / totalEdgeWeight));
		double delta_eta_cNou =  (1.0 - eta_cNou) * (1.0 - eta_cNou) * CC / (oneMinusGamma + (gamma * eta_cNou));
		double eta_u = ratio * (1.0 - (thisDegree / totalEdgeWeight));
		double delta_eta_u = (1.0 - eta_u) * (1.0 - eta_u) * CC / (oneMinusGamma + (gamma * eta_u));		
		double eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
		double delta_eta_c = (1.0 - eta_c) * (1.0 - eta_c) * CC / (oneMinusGamma + (gamma * eta_c));

		clusterIdx = clusterIdxs[thisCluster];
		clusterIdxs[thisCluster] = -1;
//		double deltaQ_exit = (incident_weights[clusterIdx] / totalEdgeWeight) + eta;
		double deltaQ_exit = (incident_weights[clusterIdx] / totalEdgeWeight) + delta_eta_cNou + delta_eta_u - delta_eta_c;
		incident_weights[clusterIdx] = 0;
		deltaQ_exit = -1 * deltaQ_exit;

		// try merging the node to each of its distinct neighbor clusters
		int targetClusterID;
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;

		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			targetClusterID = possibleNeighborClusters[i];
			if (thisCluster == targetClusterID) continue;

			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			
			eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
			delta_eta_c = (1.0 - eta_c) * (1.0 - eta_c) * CC / (oneMinusGamma + (gamma * eta_c));
			eta_u = ratio * (1.0 - (thisDegree / totalEdgeWeight));
			delta_eta_u = (1.0 - eta_u) * (1.0 - eta_u) * CC / (oneMinusGamma + (gamma * eta_u));
			eta_cu = ratio * (1.0 - (vol_cu / totalEdgeWeight));
			delta_eta_cu = (1.0 - eta_cu) * (1.0 - eta_cu) * CC / (oneMinusGamma + (gamma * eta_cu));

			clusterIdx = clusterIdxs[targetClusterID];
			clusterIdxs[targetClusterID] = -1;
//			double deltaQ = (incident_weights[clusterIdx] / totalEdgeWeight) + eta;
			double deltaQ = (incident_weights[clusterIdx] / totalEdgeWeight) + delta_eta_c + delta_eta_u - delta_eta_cu;
			incident_weights[clusterIdx] = 0;
			double deltaQ_local = deltaQ + deltaQ_exit;
			
			if (deltaQ_local > deltaQ_best) {
				deltaQ_best = deltaQ_local;
				best_targetCluster = targetClusterID;
			}
		}

		// if there is positive modularity gain, change the clustering
		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
			// decide to move
			incPerIteration += deltaQ_best;
			cluster_volume[thisCluster] -= thisDegree;
			cluster_volume[best_targetCluster] += thisDegree;
			clusters[thisNode] = best_targetCluster;

			return;
		}

		return;
	}

	private double expect(double vol_c) {
		double eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
		return  ((1.0 - eta_c) * (1.0 - eta_c)) * (1 - gamma) / (1.0 - gamma + (gamma * eta_c));
	}

	public double adjustWeight_LinearLog(double edgeWeight, double thisRatio) {
		double log = Math.log((1.0 / thisRatio) + 1) / Math.log(2);
		double linear_log = thisRatio * (1.0 / log);

		return edgeWeight * linear_log;
	}

	public void rebuildGraph(boolean firstReConstruct) {

		int arrayIndex = 0;

		// re-label the clusters, calculate the number of nodes in new hypergraph
		new_n = 0;
		int cluster;
		for (int node = 0; node < n; node++) {
			cluster = clusters[node];

			if (clusterIdxs[cluster] == -1) {
				clusterIdxs[cluster] = new_n++;
				Hypergraph.array[arrayIndex++] = cluster;
			}
		}

		// initialize the variables for graph compression
		if (firstReConstruct) {
			nodesInCluster = new ArrayList[new_n];
			newINC_eIDs = new ArrayList[new_n];
			newINC_weights = new ArrayList[new_n];
			self_loop_weights = new double[new_n];
		}

		for (int i = 0; i < new_n; i++) {
			nodesInCluster[i] = new ArrayList<Integer>();
			newINC_eIDs[i] = new ArrayList<Integer>();
			newINC_weights[i] = new ArrayList<Double>();
			node_degrees[i] = 0;
		}

		// group the nodes in each cluster
		for (int node = 0; node < n; node++) {
			cluster = clusterIdxs[clusters[node]];
			clusters[node] = cluster;
			nodesInCluster[cluster].add(node);
		}
		
		for (int i = 0; i < arrayIndex; i++) {
			clusterIdxs[Hypergraph.array[i]] = -1;
		}

		// update global_cluster
		if (save) {
			int prevID, newID;
			for (int i = 0; i < global_cluster.length; i++) {
				prevID = global_cluster[i];
				newID = clusters[prevID];
				global_cluster[i] = newID;
			}
		}

		// get the hyperedges in the compressed graph
		int newEINC_headArrIndex = 0;
		int newEINC_nIDArrIndex = 0;
		int newEINC_weightArrIndex = 0;

		// for each hyperedge
		int edgeCnt = 0;
		int edge_head = 0;
		int second_idx;
		double weight;
		for (int edge = 0; edge < m; edge++) {
			// for each old node in this hyperedge
			if (edge == endM) second_idx = EINC_nID_length;
			else second_idx = EINC_head[edge + 1];

			// store the new nodes' ID
			arrayIndex = 0;
			for (int i = EINC_head[edge]; i < second_idx; i++) {
				// match old node to new node
				cluster = clusters[EINC_nID[i]];

				if (clusterIdxs[cluster] == -1) {
					Hypergraph.array[arrayIndex++] = cluster;
					clusterIdxs[cluster] = 0;
				}
				cluster_weights[cluster] += EINC_weight[i];
			}

			// if the hyperedge contains more than one distinct new nodes
			if (arrayIndex > 1) {
				newEINC_head[newEINC_headArrIndex++] = edge_head;

				int c;
				for (int i = 0; i < arrayIndex; i++) {
					c = Hypergraph.array[i];
					clusterIdxs[c] = -1;

					weight = cluster_weights[c];
					cluster_weights[c] = 0;

					node_degrees[c] += weight;

					// new node "c" has incident hyperedge with ID "edgeCnt"
					newINC_eIDs[c].add(edgeCnt);
					newINC_weights[c].add(weight);

					newEINC_nID[newEINC_nIDArrIndex++] = c;
					newEINC_weight[newEINC_weightArrIndex++] = weight;
					edge_head++;
				}

				edgeCnt++;

			} else {
				// if the hyperedge contains only one distinct new node
				cluster = Hypergraph.array[0];
				clusterIdxs[cluster] = -1;

				weight = cluster_weights[cluster];
				cluster_weights[cluster] = 0;

				self_loop_weights[cluster] += weight;
			}
		}

		// iterate self-loop edge
		double self_loop_weight;
		for (cluster = 0; cluster < new_n; cluster++) {
			if (self_loop_weights[cluster] > 0) {

				self_loop_weight = self_loop_weights[cluster];
				self_loop_weights[cluster] = 0;

				node_degrees[cluster] += self_loop_weight;

				// new node "cluster" has incident hyperedge with ID "edgeCnt"
				newINC_eIDs[cluster].add(edgeCnt);
				newINC_weights[cluster].add(self_loop_weight);
				edgeCnt++;

				// self-loop hyperedge
				newEINC_head[newEINC_headArrIndex++] = edge_head;
				newEINC_nID[newEINC_nIDArrIndex++] = cluster;
				newEINC_weight[newEINC_weightArrIndex++] = self_loop_weight;
				edge_head++;
			}
		}

		// update n
		n = new_n;
		endN = n - 1;

		// update m
		m = newEINC_headArrIndex;
		endM = m - 1;
		
		int totalCardinality = edge_head;

		// update EINC_nID_length
		EINC_nID_length = totalCardinality;

		// update EINC_head
		for (int i = 0; i < newEINC_headArrIndex; i++) {
			EINC_head[i] = newEINC_head[i];
		}

		// update EINC_nID
		for (int i = 0; i < newEINC_nIDArrIndex; i++) {
			EINC_nID[i] = newEINC_nID[i];
		}
		
		// update EINC_weight
		for (int i = 0; i < newEINC_weightArrIndex; i++) {
			EINC_weight[i] = newEINC_weight[i];
		}

		// update edge_weights
		for (int edge = 0; edge < m; edge++) {
			// for each node in this hyperedge
			if (edge == endM) second_idx = EINC_nID_length;
			else second_idx = EINC_head[edge + 1];

			weight = 0;
			for (int i = EINC_head[edge]; i < second_idx; i++) weight += EINC_weight[i];
			edge_weights[edge] = weight;
		}

		// update INC_head, INC_eID, INC_weight
		INC_eID_length = totalCardinality;

		int node_head = 0;
		int edge;
		for (int node = 0; node < n; node++) {
			INC_head[node] = node_head;

			for (int i = 0; i < newINC_eIDs[node].size(); i++) {
				edge = newINC_eIDs[node].get(i);
				weight = newINC_weights[node].get(i);

				INC_eID[node_head] = edge;
				INC_weight[node_head] = weight;
				node_head++;
			}
		}
	}

	public void constructADJ(boolean firstConstruct) {

		// construct the adjacent list 
		
		// initialize the variables for constructing the adjacent list
		int lastHead = 0;
		int head = 0;
		int arrayIndex = 0;
		if (firstConstruct) {
			ADJ_head = new int[n];
		}

		// for each node
		int maxDegree = -1;
		int secondIdx_INC, secondIdx_EINC, edgeID, neighbor, degree;
		for (int thisNode = 0; thisNode < n; thisNode++) {

			ADJ_head[thisNode] = head;
			lastHead = head;

			// for each incident edge
			if (thisNode == endN) secondIdx_INC = INC_eID_length;
			else secondIdx_INC = INC_head[thisNode + 1];
			
			for (int i = INC_head[thisNode]; i < secondIdx_INC; i++) {
				edgeID = INC_eID[i];

				// for each node in the edge
				if (edgeID == endM) secondIdx_EINC = EINC_nID_length;
				else secondIdx_EINC = EINC_head[edgeID + 1];

				for (int j = EINC_head[edgeID]; j < secondIdx_EINC; j++) {
					neighbor = EINC_nID[j];
					if (neighbor == thisNode) continue;
					
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						Hypergraph.array[arrayIndex++] = neighbor;
						head++;
					}
				}
			}
			
			degree = head - lastHead;
			if (maxDegree < degree) maxDegree = degree;
			
			for (int i = lastHead; i < head; i++) {
				visited[Hypergraph.array[i]] = false;
			}
		}

		if (firstConstruct || arrayIndex > ADJ_nID.length) ADJ_nID = new int[arrayIndex];
		ADJ_nID_length = arrayIndex;
		
		// store the adjacent list
		for (int i = 0; i < arrayIndex; i++) {
			ADJ_nID[i] = Hypergraph.array[i];
		}

		// initialize the variables for the later node move using the information of 
		// the maximum number of adjacent neighbor
		if (firstConstruct || (maxDegree + 1) > fractions.length) {
			fractions = new double[maxDegree + 1];
			incident_clusteridxs = new int[maxDegree + 1];
			incident_weights = new double[maxDegree + 1];
			possibleNeighborClusters = new int[maxDegree + 1];
			impossibleNeighborClusters = new int[maxDegree + 1];
		}
	}
	
	public void constructADJ_largeGraph(boolean firstConstruct) {
		
		// construct the adjacency matrix 
		// for a very large hypergraph, the adjacent list needs to be splitted into
		// several smaller one
		
		int maxSplitNum = 32;
		System.out.println("#constructADJ_largeGraph");
		
		// initialize the variables for constructing the adjacent list
		byte splitCnt = 0;
		int lastHead = 0;
		int maxArraySize = Constant.maxArraySize;
		int arrayIndex = 0;
		if (firstConstruct) {
			largeSplit = new byte[n];
			largeADJ_head = new int[maxSplitNum][];
			largeADJ_nID = new int[maxSplitNum][];
			largeADJ_nID_length = new int[maxSplitNum];
			largeEndNode = new int[maxSplitNum];
			
			largeADJ_head[splitCnt] = new int[n];
		}
		
		// for each node
		boolean full = false;
		int maxDegree = -1;
		int secondIdx_INC, secondIdx_EINC, edgeID, neighbor, degree;
		for (int thisNode = 0; thisNode < n; thisNode++) {
			
			largeADJ_head[splitCnt][thisNode] = arrayIndex;
			lastHead = arrayIndex;
			
			// for each incident edge
			if (thisNode == endN) secondIdx_INC = INC_eID_length;
			else secondIdx_INC = INC_head[thisNode + 1];
			
			for (int i = INC_head[thisNode]; i < secondIdx_INC; i++) {
				edgeID = INC_eID[i];

				// for each node in the edge
				if (edgeID == endM) secondIdx_EINC = EINC_nID_length;
				else secondIdx_EINC = EINC_head[edgeID + 1];

				for (int j = EINC_head[edgeID]; j < secondIdx_EINC; j++) {
					neighbor = EINC_nID[j];
					if (neighbor == thisNode) continue;
					
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						Hypergraph.array[arrayIndex++] = neighbor;
						
						if (arrayIndex >= maxArraySize) {
							full = true;
							break;
						}
					}
				}
				
				if (full) break;
			}
			
			if (full) {
				System.out.println("full, lastHead " + lastHead);
				
				for (int i = lastHead; i < arrayIndex; i++) {
					visited[Hypergraph.array[i]] = false;
				}
				
				largeEndNode[splitCnt] = thisNode - 1;
				if (firstConstruct || lastHead > largeADJ_nID[splitCnt].length) {
					largeADJ_nID[splitCnt] = new int[lastHead];
				}
				largeADJ_nID_length[splitCnt] = lastHead;
				
				for (int i = 0; i < lastHead; i++) {
					largeADJ_nID[splitCnt][i] = Hypergraph.array[i];
				}
				
				System.out.println("split " + splitCnt + " end node " + largeEndNode[splitCnt] + " nID_length " + largeADJ_nID_length[splitCnt]);
				
				////////////////////////////////////////////////////////
				
				splitCnt++;
				thisNode--;
				full = false;
				
				////////////////////////////////////////////////////////
				
				lastHead = 0;
				arrayIndex = 0;
				
				largeADJ_head[splitCnt] = new int[n];
				
			} else {
				largeSplit[thisNode] = splitCnt;
				
				degree = arrayIndex - lastHead;
				if (maxDegree < degree) maxDegree = degree;
				
//				System.out.println("node " + thisNode + " split " + largeSplit[thisNode] + " node size " + nodeSize + " degree " + degree);
				
				for (int i = lastHead; i < arrayIndex; i++) {
					visited[Hypergraph.array[i]] = false;
				}
			}
		}
		for (int i = lastHead; i < arrayIndex; i++) {
			visited[Hypergraph.array[i]] = false;
		}
		largeEndNode[splitCnt] = n - 1;
		if (firstConstruct || arrayIndex > largeADJ_nID[splitCnt].length) {
			largeADJ_nID[splitCnt] = new int[arrayIndex];
		}
		
		// store the adjacent list
		largeADJ_nID_length[splitCnt] = arrayIndex;
		for (int i = 0; i < arrayIndex; i++) {
			largeADJ_nID[splitCnt][i] = Hypergraph.array[i];
		}
		
		System.out.println("split " + splitCnt + " end node " + largeEndNode[splitCnt] + " nID_length " + largeADJ_nID_length[splitCnt]);
		
		long totVol = 0;
		for (int split = 0; split <= splitCnt; split++) {
			totVol += largeADJ_nID_length[split];
//			System.out.println("split " + split + " largeADJ_nID_length[split] " + largeADJ_nID_length[split]);
		}
		System.out.println("totVol " + totVol + " maxDegree " + maxDegree);
		
		// initialize the variables for the later node move using the information of 
		// the maximum number of adjacent neighbor
		if (firstConstruct || (maxDegree + 1) > fractions.length) {
			fractions = new double[maxDegree + 1];
			incident_clusteridxs = new int[maxDegree + 1];
			incident_weights = new double[maxDegree + 1];
			possibleNeighborClusters = new int[maxDegree + 1];
			impossibleNeighborClusters = new int[maxDegree + 1];
		}
	}

	public String pic() throws Exception {

		// PIC follows the Louvain-style framework that iteratively maximizes the PI modularity
		
		Random random = new Random();

		int count_aggregations = 0;	// count the # of rounds

		System.out.println("PIC_noPrune_speed " + dataset + " " + lambda + " " + gamma);
		System.out.println(trial + " " + ordering + " " + ratio + " " + save);

		// variables for shuffling the order of nodes
		int randIdx, randNode;

		// variables for rebuildGraph
		boolean firstReConstruct = true;

		// variables for move
		double increase_total = 0;	// calculate the total modularity gain during iterations (a round)
		boolean hasIncrease;	// true if there is positive modularity gain in an iteration

		double startTime, endTime;
		double extractClusterTime = 0;
		do {
			System.out.println("round " + count_aggregations);
			count_aggregations++;

			////////////////////////////////////////////////////////////////////////
			// shuffle the node order

			startTime = System.currentTimeMillis();

			switch (ordering) {
			case "randomOrder": {
				// shuffle the node order
				for (int i = 0; i < n; i++) {
					randIdx = random.nextInt(n);
					randNode = nodeInOrder[i];
					nodeInOrder[i] = nodeInOrder[randIdx];
					nodeInOrder[randIdx] = randNode;
				}
			}
				break;
			}

			//////////////////////////////////////////////////////////////////////
			// iteratively mobilize each node from its own cluster to its neighbors' clusters
			
			iteration = 0;	// reset the number of iterations
			curTimeStamp = 0;	// reset the number of node moves
			increase_total = 0;	// reset the total modularity gain
			hasIncrease = true;
			
			while (hasIncrease) {
				iteration++;
				hasIncrease = false;
				incPerIteration = 0;	// calculate the total modularity gain in an iteration
				
				// each iteration scans all the nodes in the graph
				for (int idx = 0; idx < nodeInOrder.length; idx++) {
					curTimeStamp++;

					// move
					move_split(nodeInOrder[idx]);		// uncomment this for normal-sized hypergraph
//					move_split_largeGraph(nodeInOrder[idx]);	// uncomment this for very large hypergraph
					
				}
				
				increase_total += incPerIteration;
				// whether the total modularity gained in an iteration is negligible
				if (incPerIteration > eps) hasIncrease = true;
				// whether reaching the maximum number of iterations
				if (iteration > n_aggregations) break;
			}
				
			//////////////////////////////////////////////////////////////////////
			// compress all the nodes in each cluster into a supernode; update the graph and construct new adjacent list

			if (increase_total <= eps) break;

			// compressed graph, compress all the nodes in each cluster into a supernode
			rebuildGraph(firstReConstruct);
			
			// construct the adjacent list 
			constructADJ(false);		// uncomment this for normal-sized hypergraph
//			constructADJ_largeGraph(false);	// uncomment this for very large hypergraph
			
			firstReConstruct = false;

			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);

			/////////////////////////////////////////////////////////////////////

			// calculate gamma
			double timeUpdateGamma = calculateGamma();
			extractClusterTime += (timeUpdateGamma * Constant.RUNNING_TIME_UNIT);

			/////////////////////////////////////////////////////////////////////

			// whether reaching the maximum number of rounds
			if (count_aggregations > n_aggregations) break;

			startTime = System.currentTimeMillis();
			// reset node moving orders, clusters and cluster volumes
			nodeInOrder = new int[n];
			for (int i = 0; i < n; i++) {
				nodeInOrder[i] = i;
				clusters[i] = i;
				cluster_volume[i] = node_degrees[i];
			}
			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);

		} while (true);

		if (save) {
			// save the clustering result to a txt file
			saveClusters();
		}

		return extractClusterTime + "";
	}

	public double calculateGamma() {
	
		// calculate the gamma
		
		int cardinality;
		int nonTrivialEdgeNum = 0;
		double time = System.currentTimeMillis();
		for (int edgeID = 0; edgeID < m; edgeID++) {
			if (edgeID == endM) cardinality = (EINC_nID_length - EINC_head[edgeID]);
			else cardinality = (EINC_head[edgeID + 1] - EINC_head[edgeID]);

			if (cardinality >= 2) nonTrivialEdgeNum++;
		}
		int trivialEdgeNum = m - nonTrivialEdgeNum;
		
		double gamma = (double) ((EINC_nID_length-trivialEdgeNum) - 2 * nonTrivialEdgeNum) / (double) ((EINC_nID_length-trivialEdgeNum) - nonTrivialEdgeNum);
		
		double runningTime = (System.currentTimeMillis() - time) / Constant.RUNNING_TIME_UNIT;
//		System.out.println("calculateGamma " + gamma);
		
		this.gamma = gamma;
		return runningTime;
	}
	
	public void saveClusters() {

		// save the clustering result to a txt file
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_cluster_pic_move"
					+ "_ordered_" + "false_order_" + this.ordering + "_ratio_" + ratio + "_trial_"
					+ (this.trial) + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_cluster_pic_move"
					+ "_ordered_" + "false_order_" + this.ordering + "_ratio_" + ratio
					+ "_connect_trial_" + (this.trial) + ".txt";
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
					fwCount.write(
							nodesInCluster[i].stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
				}
			}

			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
}
