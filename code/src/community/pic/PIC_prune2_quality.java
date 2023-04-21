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

// hypergraph clustering using PIC with optimization technique 2
public class PIC_prune2_quality {

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

	// node - hyperedge matrix
	double totalEdgeWeight;
	double[] INC_weight;
	int[] INC_eID;
	int[] INC_head;
	int INC_eID_length;

	// hyperedge - node matrix
	int[] EINC_nID;
	double[] EINC_weight;
	int[] EINC_head;
	int EINC_nID_length;
	double[] edge_weights;

	double[] maxFractionsLazy;
	double[] maxFractionChanges;

	// adjacent list
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

	// variables for clusterings
	List<Integer>[] nodesInCluster;
	int[] global_cluster;
	int[] clusters;
	double[] cluster_volume;

	int[] nodeInOrder;
	int[] nodePriority;

	// variables for rebuildGraph (graph compression)
	int[] nodeArray;
	int[] newEINC_head;
	int[] newEINC_nID;
	double[] newEINC_weight;
	List<Integer>[] newINC_eIDs;
	List<Double>[] newINC_weights;
	int[] clusterIdxs;
	double[] cluster_weights;
	double[] self_loop_weights;
	int new_n;

	boolean usePruning = true;
	boolean save;
	String dataset;

	// PIC algorithm with optimization technique 2
	public PIC_prune2_quality(int trial, String ordering, double ratio,
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
		
		// number of node and hyperedges
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
		constructADJ(true);
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
		
		// calculate hyperedge weights
		totalEdgeWeight = 0;
		edge_weights = new double[m];

		double weight;
		int second_idx;
		for (int edge = 0; edge < m; edge++) {
			// for each node in the edge
			if (edge == endM) second_idx = EINC_nID_length;
			else second_idx = EINC_head[edge + 1];

			weight = 0;
			for (int k = EINC_head[edge]; k < second_idx; k++) weight += EINC_weight[k];
			totalEdgeWeight += weight;
			edge_weights[edge] = weight;
		}

		// calculate node weights
		node_degrees = new double[n];

		// for each node
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
		nodeArray = new int[n];
		newEINC_head = new int[m];
		newEINC_nID = new int[EINC_nID.length];
		newEINC_weight = new double[EINC_weight.length];
		cluster_weights = new double[n];
		clusterIdxs = new int[n];
		for (int i = 0; i < n; i++) clusterIdxs[i] = -1;

		////////////////////////////////////////////////////////////////////////////
		
		// initialize the node moving order, clustering and cluster volumes
		nodeInOrder = new int[n];
		clusters = new int[n];
		cluster_volume = new double[n];
		for (int i = 0; i < n; i++) {
			nodeInOrder[i] = i;
			clusters[i] = i;
			cluster_volume[i] += node_degrees[i];
		}

		///////////////////////////////////////////////////////////////////////////////

		// use optimization technique
		if (usePruning) {
			ES_construct();
		}
		
		runningTime += (System.currentTimeMillis() - time);
		
		return runningTime;
	}

	public String move_split_prune2(int thisNode) {

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
			return "";
		}
		
		int thisClusterIdx = clusterIdxs[thisCluster];
		int incident_clusterNum = 0;

		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == endN) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];

		// calculate the delta support (using optimization technique 2)
		int edge, node, clusterIdx, secondIdx;
		double thisWeight, edgeWeight;
		double maxFraction, fraction_cNotu;
		double ratio_u, ratio_cNotu, ratio_c;
		double contribute_u, contribute_cNotu, contribute_c;
		for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
			edge = INC_eID[i];
			thisWeight = INC_weight[i];
			edgeWeight = edge_weights[edge];

			// maxRatioLazy_c = (thisWeight + maxFractionsLazy[edge] + maxFractionChanges[edge]) / edgeWeight;
			if ((thisWeight + maxFractionsLazy[edge] + maxFractionChanges[edge]) / edgeWeight < ratio) continue;

			// for each node in the edge
			if (edge == endM) secondIdx = EINC_nID_length;
			else secondIdx = EINC_head[edge + 1];

			incident_clusterNum = 0;
			for (int j = EINC_head[edge]; j < secondIdx; j++) {
				node = EINC_nID[j];
				if (node == thisNode) continue;

				clusterIdx = clusterIdxs[clusters[node]];
				fractions[clusterIdx] += EINC_weight[j];
				if (!visited[clusterIdx]) {
					incident_clusteridxs[incident_clusterNum++] = clusterIdx;
					visited[clusterIdx] = true;
				}
			}

			if (maxFractionChanges[edge] > 0) {
				//needUpdate
				
				maxFraction = thisWeight;
				for (int j = 0; j < incident_clusterNum; j++) {
					clusterIdx = incident_clusteridxs[j];

					if (clusterIdx == thisClusterIdx) fraction_cNotu = fractions[clusterIdx] + thisWeight;
					else fraction_cNotu = fractions[clusterIdx];

					if (maxFraction < fraction_cNotu) maxFraction = fraction_cNotu;
				}
				maxFractionsLazy[edge] = maxFraction;
				maxFractionChanges[edge] = 0;

				// maxRatio_c = (thisWeight + maxFractionsLazy[edge]) / edgeWeight;
				if ((thisWeight + maxFractionsLazy[edge]) / edgeWeight < ratio) {
					for (int j = 0; j < incident_clusterNum; j++) {
						clusterIdx = incident_clusteridxs[j];
						visited[clusterIdx] = false;
						fractions[clusterIdx] = 0;
					}
					continue;
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

				ratio_c = ratio_u + ratio_cNotu;
				if (ratio_c >= ratio) {
					contribute_c = adjustWeight_LinearLog(edgeWeight, ratio_c);
					if (ratio_cNotu >= ratio) contribute_cNotu = adjustWeight_LinearLog(edgeWeight, ratio_cNotu);
					else contribute_cNotu = 0;
				} else {
					contribute_c = 0;
					contribute_cNotu = 0;
				}
				
				incident_weights[clusterIdx] += (contribute_c - contribute_u - contribute_cNotu);
			}
		}
		
		// try removing the node from current cluster
		double vol_cu, eta_cu, delta_eta_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		
		double eta_cNou = ratio * (1.0 - (vol_cNou / totalEdgeWeight));
		double delta_eta_cNou = ((1.0 - eta_cNou) * (1.0 - eta_cNou)) * (1 - gamma) / (1.0 - gamma + (gamma * eta_cNou));
		double eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
		double delta_eta_c = ((1.0 - eta_c) * (1.0 - eta_c)) * (1 - gamma) / (1.0 - gamma + (gamma * eta_c));
		
		clusterIdx = clusterIdxs[thisCluster];
		clusterIdxs[thisCluster] = -1;
		double incident_weight = incident_weights[clusterIdx];
		incident_weights[clusterIdx] = 0;
		double deltaQ_exit_short = (incident_weight / totalEdgeWeight) + delta_eta_cNou - delta_eta_c;
		deltaQ_exit_short = -1 * deltaQ_exit_short;

		// try merging the node to each of its distinct neighbor clusters
		int targetClusterID;
		double deltaQ_short, deltaQ_local_short;
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		// for each neighbor cluster
		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			targetClusterID = possibleNeighborClusters[i];
			if (thisCluster == targetClusterID) continue;

			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			
			eta_c = ratio * (1.0 - (vol_c / totalEdgeWeight));
			delta_eta_c = ((1.0 - eta_c) * (1.0 - eta_c)) * (1 - gamma) / (1.0 - gamma + (gamma * eta_c));
			eta_cu = ratio * (1.0 - (vol_cu / totalEdgeWeight));
			delta_eta_cu = ((1.0 - eta_cu) * (1.0 - eta_cu)) * (1 - gamma) / (1.0 - gamma + (gamma * eta_cu));		

			clusterIdx = clusterIdxs[targetClusterID];
			clusterIdxs[targetClusterID] = -1;
			deltaQ_short = (incident_weights[clusterIdx] / totalEdgeWeight) + delta_eta_c - delta_eta_cu;
			incident_weights[clusterIdx] = 0;
			deltaQ_local_short = deltaQ_short + deltaQ_exit_short;
			
			// update the best modularity gain and cluster
			if (deltaQ_local_short > deltaQ_best) {
				deltaQ_best = deltaQ_local_short;
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

			for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
				edge = INC_eID[i];
				maxFractionChanges[edge] += INC_weight[i];
			}
			
            return deltaQ_best + "," + best_targetCluster;
		}
		
		return "";
	}

	public double adjustWeight_AON(double edgeWeight, double thisRatio) {
		if (thisRatio >= 1.0) return edgeWeight;
		return 0;
	}

	public double adjustWeight_LinearLog(double edgeWeight, double thisRatio) {
		return edgeWeight * thisRatio * Math.log(2) / Math.log((1.0 / thisRatio) + 1);
	}

	public double adjustWeight_Quadratic(double edgeWeight, double thisRatio) {
		double fractionPow = Math.pow(thisRatio, 2);
		return edgeWeight * fractionPow;
	}

	public double adjustWeight_Exp(double edgeWeight, double thisRatio) {
		double fractionExp = (Math.exp(thisRatio)-1) / (Math.exp(1)-1);
		return edgeWeight * fractionExp;
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
				nodeArray[arrayIndex++] = cluster;
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
			clusterIdxs[nodeArray[i]] = -1;
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
					nodeArray[arrayIndex++] = cluster;
					clusterIdxs[cluster] = 0;
				}
				cluster_weights[cluster] += EINC_weight[i];
			}

			// if the hyperedge contains more than one distinct new nodes
			if (arrayIndex > 1) {
				newEINC_head[newEINC_headArrIndex++] = edge_head;

				int c;
				for (int i = 0; i < arrayIndex; i++) {
					c = nodeArray[i];
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
				cluster = nodeArray[0];
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
		List<Integer> adj_nIDs = new ArrayList<Integer>();
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
			if (thisNode == n - 1) secondIdx_INC = INC_eID_length;
			else secondIdx_INC = INC_head[thisNode + 1];
			
			for (int i = INC_head[thisNode]; i < secondIdx_INC; i++) {
				edgeID = INC_eID[i];

				// for each node in the edge
				if (edgeID == m - 1) secondIdx_EINC = EINC_nID_length;
				else secondIdx_EINC = EINC_head[edgeID + 1];

				for (int j = EINC_head[edgeID]; j < secondIdx_EINC; j++) {
					neighbor = EINC_nID[j];
					if (neighbor == thisNode) continue;
					
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						adj_nIDs.add(neighbor);
						head++;
					}
				}
			}
			
			degree = head - lastHead;
			if (maxDegree < degree) maxDegree = degree;
			
			for (int i = lastHead; i < head; i++) {
				visited[adj_nIDs.get(i)] = false;
			}
		}

		
		if (firstConstruct || adj_nIDs.size() > ADJ_nID.length) {
			ADJ_nID = new int[adj_nIDs.size()];
		}
		ADJ_nID_length = adj_nIDs.size();

		// store the adjacent list
		int idx = 0;
		for (int adj_nID : adj_nIDs) {
			ADJ_nID[idx++] = adj_nID;
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
	
	public String pic() throws Exception {

		// PIC follows the Louvain-style framework that iteratively maximizes the PI modularity
		
		Random random = new Random();

		int count_aggregations = 0;	// count the # of rounds

		System.out.println("PIC_prune2_quality " + dataset + " " + lambda + " " + gamma);
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
			increase_total = 0;	// reset the total modularity gain
			hasIncrease = true;
			
			while (hasIncrease) {
				iteration++;
				hasIncrease = false;
				incPerIteration = 0;	// calculate the total modularity gain in an iteration
				
				// each iteration scans all the nodes in the graph
				for (int idx = 0; idx < nodeInOrder.length; idx++) {

					// move
					move_split_prune2(nodeInOrder[idx]);
					
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
			constructADJ(false);
			
			// refresh the variables for optimization techniques
			if (usePruning) {
				ES_refresh();
			}
			
			firstReConstruct = false;

			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);

			/////////////////////////////////////////////////////////////////////

			// calculate gamma
			double updateGammaTime = calculateGamma();
			extractClusterTime += (updateGammaTime * Constant.RUNNING_TIME_UNIT);

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

		System.out.println(String.format("%.4f", extractClusterTime));

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

	private void ES_construct() {

		maxFractionsLazy = new double[m];
		maxFractionChanges = new double[m];

		for (int edge = 0; edge < m; edge++) {
			maxFractionsLazy[edge] = EINC_weight[EINC_head[edge]];
		}
	}

	private void ES_refresh() {

		int firstIdx_EINC, secondIdx_EINC;

		double maxFraction, fraction;
		for (int edge = 0; edge < m; edge++) {

			// for each node in the edge
			firstIdx_EINC = EINC_head[edge];
			if (edge == endM) secondIdx_EINC = EINC_nID_length;
			else secondIdx_EINC = EINC_head[edge + 1];

			maxFraction = -1;
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				fraction = EINC_weight[i];
				if (maxFraction < fraction) {
					maxFraction = fraction;
				}
			}

			maxFractionsLazy[edge] = maxFraction;
			maxFractionChanges[edge] = 0;
		}
	}
}
