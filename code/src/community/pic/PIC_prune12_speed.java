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

public class PIC_prune12_speed {

	int trial;
	double ratio;
	boolean toHigherOrder;
	String ordering;

	double increasePerPass;
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

	double[] maxFractionsLazy;
	double[] maxFractionChanges;
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
	int round;
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

	// variables for updateGamma
	Process process;

	// count-min sketch
	boolean useSketch = true;

	boolean save;
	String moveStrategy;
	String dataset;

	public PIC_prune12_speed(int trial, boolean toHigherOrder, String ordering, double ratio,
			boolean save) throws IOException {
		this.trial = trial;
		this.toHigherOrder = toHigherOrder;
		this.ordering = ordering;
		this.ratio = ratio;
		this.save = save;
		this.dataset = Hypergraph.dataset;
	}

	public PIC_prune12_speed(int trial, boolean toHigherOrder, String ordering, double ratio,
			double lambda, boolean save) throws IOException {
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

		int cardinality = Hypergraph.EINC_nID.length;
		EINC_nID = new int[cardinality];
		EINC_weight = new double[cardinality];
		INC_eID = new int[cardinality];
		INC_weight = new double[cardinality];
		for (int i = 0; i < cardinality; i++) {
			EINC_nID[i] = Hypergraph.EINC_nID[i];
			EINC_weight[i] = Hypergraph.EINC_weight[i];
			INC_eID[i] = Hypergraph.INC_eID[i];
			INC_weight[i] = Hypergraph.INC_weight[i];
		}

		EINC_nID_length = cardinality;
		INC_eID_length = cardinality;

		// variables for move
		visited = new boolean[n];
		constructADJ(true);
//		constructADJ_largeGraph(true);
	}

	public double initPart2() {
		
//		if (!getGamma(m)) System.out.println("error!");
				
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

		if (toHigherOrder) {
			nodePriority = new int[n];
		}

		///////////////////////////////////////////////////////////////////////////////

		if (useSketch) {
			NS_construct();
			ES_construct();
		}
		
		runningTime += (System.currentTimeMillis() - time);
		
		return runningTime;
	}

	public void move_split_prune12(int thisNode) {

		int thisCluster = clusters[thisNode];

//		System.out.println("move node " + thisNode + "(" + thisCluster + ")");

		int secondIdx_ADJ;
		if (thisNode == endN) secondIdx_ADJ = ADJ_nID_length;
		else secondIdx_ADJ = ADJ_head[thisNode + 1];

		int cluster;
		boolean hasImpossible = false;
		int distinctNeighborClusterNum = 0;
		int possibleNeighborClusterNum = 0;
		int impossibleNeighborClusterNum = 0;
		
		if (round > 1) {
			if ((curTimeStamp - lastChangeTimes[thisCluster]) > n) {
				// add cluster of thisNode
				clusterIdxs[thisCluster] = distinctNeighborClusterNum++;
				// add cluster of adjacent neighbors
				for (int i = ADJ_head[thisNode]; i < secondIdx_ADJ; i++) {
					cluster = clusters[ADJ_nID[i]];
					if (clusterIdxs[cluster] == -1) {
						
						if ((curTimeStamp - lastChangeTimes[cluster]) > n) {
							impossibleNeighborClusters[impossibleNeighborClusterNum++] = cluster;
						} else {
							possibleNeighborClusters[possibleNeighborClusterNum++] = cluster;
						}
						
						clusterIdxs[cluster] = distinctNeighborClusterNum++;
					}
				}
				
				if (impossibleNeighborClusterNum > 0) hasImpossible = true;
				
				if (possibleNeighborClusterNum > 0) {
					// canSkip = false;
				} else {
					clusterIdxs[thisCluster] = -1;
					for (int i = 0; i < impossibleNeighborClusterNum; i++) {
						cluster = impossibleNeighborClusters[i];
						clusterIdxs[cluster] = -1;
					}
					
					return;
				}
			} else {
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
			}
			
		} else {
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
		}
		
		if (distinctNeighborClusterNum == 1) {
			clusterIdxs[thisCluster] = -1;
			
			return;
		}
		
		int thisClusterIdx = clusterIdxs[thisCluster];
		int incident_clusterNum = 0;

		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == endN) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];

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
		
		// try removing from original cluster
		double vol_cu, one_minus_gammaExpect_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		double one_minus_gammaExpect_cNou = 1 - (gamma * ((((vol_cNou / totalEdgeWeight) - 1) * ratio) + 1));
		double one_minus_gammaExpect_c = 1 - (gamma * ((((vol_c / totalEdgeWeight) - 1) * ratio) + 1));
		double eta_short = (1.0 / one_minus_gammaExpect_cNou) - (1 / one_minus_gammaExpect_c);
		
		clusterIdx = clusterIdxs[thisCluster];
		clusterIdxs[thisCluster] = -1;
		double incident_weight = incident_weights[clusterIdx];
		incident_weights[clusterIdx] = 0;
		double deltaQ_exit_short = (incident_weight / totalEdgeWeight) + eta_short;
		deltaQ_exit_short = -1 * deltaQ_exit_short;

		int targetClusterID;
		double deltaQ_short, deltaQ_local_short;
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			targetClusterID = possibleNeighborClusters[i];
			if (thisCluster == targetClusterID) continue;

//			System.out.println("move to cluster " + targetClusterID);

			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			one_minus_gammaExpect_c = 1 - (gamma * ((((vol_c / totalEdgeWeight) - 1) * ratio) + 1)); 
			one_minus_gammaExpect_cu = 1 - (gamma * ((((vol_cu / totalEdgeWeight) - 1) * ratio) + 1)); 
			eta_short = (1.0 / one_minus_gammaExpect_c) - (1 / one_minus_gammaExpect_cu);

			clusterIdx = clusterIdxs[targetClusterID];
			clusterIdxs[targetClusterID] = -1;
			deltaQ_short = (incident_weights[clusterIdx] / totalEdgeWeight) + eta_short;
			incident_weights[clusterIdx] = 0;
			deltaQ_local_short = deltaQ_short + deltaQ_exit_short;
			if (deltaQ_local_short > deltaQ_best) {
				deltaQ_best = deltaQ_local_short;
				best_targetCluster = targetClusterID;
			}
		}
		
		if (hasImpossible) {
			for (int i = 0; i < impossibleNeighborClusterNum; i++) {
				cluster = impossibleNeighborClusters[i];
				clusterIdx = clusterIdxs[cluster];
				clusterIdxs[cluster] = -1;
				incident_weights[clusterIdx] = 0;
			}
		}

		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
			// decide to move
			increasePerPass += deltaQ_best;
			cluster_volume[thisCluster] -= thisDegree;
			cluster_volume[best_targetCluster] += thisDegree;
			clusters[thisNode] = best_targetCluster;
			
			lastChangeTimes[thisCluster] = curTimeStamp;
			lastChangeTimes[best_targetCluster] = curTimeStamp;
//            System.out.println("new cluster of " + thisNode + " is " + clusters[thisNode]);

			for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
				edge = INC_eID[i];
				maxFractionChanges[edge] += INC_weight[i];
			}

//            System.out.println("=====================\n");
			
            return;
//            return "";
		}

//		System.out.println("=====================\n");
		
		return;
	}

	public void move_split_prune12_largeGraph(int thisNode) {

		int thisCluster = clusters[thisNode];

//		System.out.println("move node " + thisNode + "(" + thisCluster + ")");

		int split = largeSplit[thisNode];
		int endNode = largeEndNode[split];
		
		int secondIdx_ADJ;
		if (thisNode == endNode) secondIdx_ADJ = largeADJ_nID_length[split];
		else secondIdx_ADJ = largeADJ_head[split][thisNode + 1];

		int cluster;
		boolean hasImpossible = false;
		int distinctNeighborClusterNum = 0;
		int possibleNeighborClusterNum = 0;
		int impossibleNeighborClusterNum = 0;
		
		if (round > 1) {
			if ((curTimeStamp - lastChangeTimes[thisCluster]) > n) {
				// add cluster of thisNode
				clusterIdxs[thisCluster] = distinctNeighborClusterNum++;
				visited[thisCluster] = true;
				// add cluster of adjacent neighbors
				for (int i = largeADJ_head[split][thisNode]; i < secondIdx_ADJ; i++) {
					cluster = clusters[largeADJ_nID[split][i]];
					if (!visited[cluster]) {
						
						if ((curTimeStamp - lastChangeTimes[cluster]) > n) {
							impossibleNeighborClusters[impossibleNeighborClusterNum++] = cluster;
						} else {
							possibleNeighborClusters[possibleNeighborClusterNum++] = cluster;
						}
						
						clusterIdxs[cluster] = distinctNeighborClusterNum++;
						visited[cluster] = true;
					}
				}
				
				if (impossibleNeighborClusterNum > 0) hasImpossible = true;
				
				if (possibleNeighborClusterNum > 0) {
					// canSkip = false;
					visited[thisCluster] = false;
					for (int i = 0; i < impossibleNeighborClusterNum; i++) {
						visited[impossibleNeighborClusters[i]] = false;
					}
					for (int i = 0; i < possibleNeighborClusterNum; i++) {
						visited[possibleNeighborClusters[i]] = false;
					}
				} else {
					visited[thisCluster] = false;
					clusterIdxs[thisCluster] = -1;
					for (int i = 0; i < impossibleNeighborClusterNum; i++) {
						cluster = impossibleNeighborClusters[i];
						visited[cluster] = false;
						clusterIdxs[cluster] = -1;
					}
					
					return;
				}
			} else {
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
			}
			
		} else {
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
		}
		
		if (distinctNeighborClusterNum == 1) {
			clusterIdxs[thisCluster] = -1;
			
			return;
		}
		
		int thisClusterIdx = clusterIdxs[thisCluster];
		int incident_clusterNum = 0;

		// for each incident edge
		int secondIdx_INC, firstIdx_INC = INC_head[thisNode];
		if (thisNode == endN) secondIdx_INC = INC_eID_length;
		else secondIdx_INC = INC_head[thisNode + 1];

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
		
		// try removing from original cluster
		double vol_cu, one_minus_gammaExpect_cu;
		double vol_c = cluster_volume[thisCluster];
		double thisDegree = node_degrees[thisNode];
		double vol_cNou = vol_c - thisDegree;
		double one_minus_gammaExpect_cNou = 1 - (gamma * ((((vol_cNou / totalEdgeWeight) - 1) * ratio) + 1));
		double one_minus_gammaExpect_c = 1 - (gamma * ((((vol_c / totalEdgeWeight) - 1) * ratio) + 1));
		double eta_short = (1.0 / one_minus_gammaExpect_cNou) - (1 / one_minus_gammaExpect_c);
		
		clusterIdx = clusterIdxs[thisCluster];
		clusterIdxs[thisCluster] = -1;
		double incident_weight = incident_weights[clusterIdx];
		incident_weights[clusterIdx] = 0;
		double deltaQ_exit_short = (incident_weight / totalEdgeWeight) + eta_short;
		deltaQ_exit_short = -1 * deltaQ_exit_short;

		int targetClusterID;
		double deltaQ_short, deltaQ_local_short;
		int best_targetCluster = thisCluster;
		double deltaQ_best = 0;
		
		for (int i = 0; i < possibleNeighborClusterNum; i++) {
			targetClusterID = possibleNeighborClusters[i];
			if (thisCluster == targetClusterID) continue;

//			System.out.println("move to cluster " + targetClusterID);

			vol_c = cluster_volume[targetClusterID];
			vol_cu = thisDegree + vol_c;
			one_minus_gammaExpect_c = 1 - (gamma * ((((vol_c / totalEdgeWeight) - 1) * ratio) + 1)); 
			one_minus_gammaExpect_cu = 1 - (gamma * ((((vol_cu / totalEdgeWeight) - 1) * ratio) + 1)); 
			eta_short = (1.0 / one_minus_gammaExpect_c) - (1 / one_minus_gammaExpect_cu);

			clusterIdx = clusterIdxs[targetClusterID];
			clusterIdxs[targetClusterID] = -1;
			deltaQ_short = (incident_weights[clusterIdx] / totalEdgeWeight) + eta_short;
			incident_weights[clusterIdx] = 0;
			deltaQ_local_short = deltaQ_short + deltaQ_exit_short;
			if (deltaQ_local_short > deltaQ_best) {
				deltaQ_best = deltaQ_local_short;
				best_targetCluster = targetClusterID;
			}
		}
		
		if (hasImpossible) {
			for (int i = 0; i < impossibleNeighborClusterNum; i++) {
				cluster = impossibleNeighborClusters[i];
				clusterIdx = clusterIdxs[cluster];
				clusterIdxs[cluster] = -1;
				incident_weights[clusterIdx] = 0;
			}
		}

		if (deltaQ_best > 0 && best_targetCluster != thisCluster) {
			// decide to move
			increasePerPass += deltaQ_best;
			cluster_volume[thisCluster] -= thisDegree;
			cluster_volume[best_targetCluster] += thisDegree;
			clusters[thisNode] = best_targetCluster;
			
			lastChangeTimes[thisCluster] = curTimeStamp;
			lastChangeTimes[best_targetCluster] = curTimeStamp;
//            System.out.println("new cluster of " + thisNode + " is " + clusters[thisNode]);

			for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
				edge = INC_eID[i];
				maxFractionChanges[edge] += INC_weight[i];
			}

//            System.out.println("=====================\n");
			
            return;
		}

//		System.out.println("=====================\n");
		
		return;
	}

	public double adjustWeight_AON(double edgeWeight, double thisRatio) {
		if (thisRatio >= 1.0) return edgeWeight;
		return 0;
	}

	public double adjustWeight_LinearLog(double edgeWeight, double thisRatio) {
		double log = Math.log((1.0 / thisRatio) + 1) / Math.log(2);
		double linear_log = thisRatio * (1.0 / log);

		return edgeWeight * linear_log;
	}

	public double adjustWeight_Sigmoid(double edgeWeight, double thisRatio) {
		double exp = Math.exp(1 - (1.0 / thisRatio));
		double sigmod = 2 * ((1.0 / (1 + exp)) - 0.5);

		return edgeWeight * (1 - sigmod);
	}

	public double adjustWeight_Quadratic(double edgeWeight, double thisRatio) {
		double fractionPow = Math.pow(thisRatio, 2);
		return edgeWeight * fractionPow;
	}

	public void rebuildGraph(boolean firstReConstruct) {

		int arrayIndex = 0;

		new_n = 0; // the number of nodes in new hypergraph
		int cluster;
		for (int node = 0; node < n; node++) {
			cluster = clusters[node];

			if (clusterIdxs[cluster] == -1) {
				clusterIdxs[cluster] = new_n++;
				Hypergraph.array[arrayIndex++] = cluster;
			}
		}

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

		// store the nodes and their weights in each hyperedge
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

			if (arrayIndex > 1) {
				newEINC_head[newEINC_headArrIndex++] = edge_head;

				int c;
				for (int i = 0; i < arrayIndex; i++) {
					c = Hypergraph.array[i];
					clusterIdxs[c] = -1;

					weight = cluster_weights[c];
					cluster_weights[c] = 0;

					node_degrees[c] += weight;

					// new node "cluster" has incident hyperedge with ID "head"
					newINC_eIDs[c].add(edgeCnt);
					newINC_weights[c].add(weight);

					newEINC_nID[newEINC_nIDArrIndex++] = c;
					newEINC_weight[newEINC_weightArrIndex++] = weight;
					edge_head++;
				}

				// new node set in hyperedge with ID "edgeCnt"
				edgeCnt++;

			} else {
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

				// new node "clusterID" has incident hyperedge with ID "head"
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

//        // check correctness
//		double totalNodeDegree = 0;
// 		for (int i = 0; i < n; i++) totalNodeDegree += node_degrees[i];
//		double totalNode_weights = 0;
// 		for (int i = 0; i < INC_eID_length; i++) totalNode_weights += INC_weight[i];
//		double totalEdge_weights = 0;
// 		for (int i = 0; i < EINC_nID_length; i++) totalEdge_weights += EINC_weight[i];
//		double totalWeights = 0;
// 		for (int i = 0; i < m; i++) totalWeights += edge_weights[i];
//		System.out.println("# new nodes " + n + " # of new hyperedges " + m + " Q_curr " + Q_curr + " total weight "
//				+ totalNodeDegree + " " + totalNode_weights + " " + totalEdge_weights + " " + totalWeights);
	}

	public void constructADJ(boolean firstConstruct) {

		int lastHead = 0;
		int head = 0;
		int arrayIndex = 0;
		if (firstConstruct) {
			ADJ_head = new int[n];
		}

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
		
		for (int i = 0; i < arrayIndex; i++) {
			ADJ_nID[i] = Hypergraph.array[i];
		}
		
//		System.out.println("totVol " + ADJ_nID_length + " maxDegree " + maxDegree);

		if (firstConstruct || (maxDegree + 1) > fractions.length) {
			fractions = new double[maxDegree + 1];
			incident_clusteridxs = new int[maxDegree + 1];
			incident_weights = new double[maxDegree + 1];
			possibleNeighborClusters = new int[maxDegree + 1];
			impossibleNeighborClusters = new int[maxDegree + 1];
		}
	}
	
	public void constructADJ_largeGraph(boolean firstConstruct) {
		
		int maxSplitNum = 32;
		System.out.println("#constructADJ_largeGraph");
		
		byte splitCnt = 0;
		int lastHead = 0;
		int maxArraySize = Hypergraph.maxArraySize;
		int arrayIndex = 0;
		if (firstConstruct) {
			largeSplit = new byte[n];
			largeADJ_head = new int[maxSplitNum][];
			largeADJ_nID = new int[maxSplitNum][];
			largeADJ_nID_length = new int[maxSplitNum];
			largeEndNode = new int[maxSplitNum];
			
			largeADJ_head[splitCnt] = new int[n];
		}
		
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
		
		if (firstConstruct || (maxDegree + 1) > fractions.length) {
			fractions = new double[maxDegree + 1];
			incident_clusteridxs = new int[maxDegree + 1];
			incident_weights = new double[maxDegree + 1];
			possibleNeighborClusters = new int[maxDegree + 1];
			impossibleNeighborClusters = new int[maxDegree + 1];
		}
	}

	public String louvain(String moveStrategy) throws Exception {

		Random random = new Random();

		int count_aggregations = 0;
		this.moveStrategy = moveStrategy;

		System.out.println(dataset + " " + lambda + " " + gamma);
		System.out.println(trial + " " + toHigherOrder + " " + ordering + " " + moveStrategy + " " + ratio + " " + save);

		// variables for ordering
		int randIdx, randNode;

		// variables for rebuildGraph
		boolean firstReConstruct = true;

		// variables for updateGamma
		double timeUpdateGamma;

		// variables for move
		double increase_total = 0;
		boolean hasIncrease;

//		double randomTime = 0;
//		double moveTime = 0;
//		double rebuildTime = 0;
//		double matchingTime = 0;
//		double otherTime = 0;
//		double timer = 0;

		double startTime, endTime;
		double extractClusterTime = 0;
		do {
			System.out.println("itr " + count_aggregations);
			count_aggregations++;

			////////////////////////////////////////////////////////////////////////

			startTime = System.currentTimeMillis();

//			timer = System.currentTimeMillis();
			switch (ordering) {
			case "randomOrder": {
				// random the order
				for (int i = 0; i < n; i++) {
					randIdx = random.nextInt(n);
					randNode = nodeInOrder[i];
					nodeInOrder[i] = nodeInOrder[randIdx];
					nodeInOrder[randIdx] = randNode;
				}
			}
				break;
			}
//			randomTime += (System.currentTimeMillis() - timer) / Constant.RUNNING_TIME_UNIT;

//			timer = System.currentTimeMillis();
			switch (moveStrategy) {
			case "move": {
				round = 0;
				curTimeStamp = 0;
				
				increase_total = 0;
				hasIncrease = true;
				while (hasIncrease) {
					round++;
					hasIncrease = false;
					increasePerPass = 0;

//					double step1Time = 0;
//					double step2Time = 0;
//					double step21Time = 0;
//					double step22Time = 0;
//					double step23Time = 0;
//					double step3Time = 0;

//					double total = 0;
//					double sum1 = 0;
//					double sum2 = 0;
//					double sum3 = 0;
					
//					System.out.println("round " + round);
					
					for (int idx = 0; idx < nodeInOrder.length; idx++) {
						curTimeStamp++;
						
//						String[] strs;
//						String str;

						move_split_prune12(nodeInOrder[idx]);
//						move_split_prune12_largeGraph(nodeInOrder[idx]);
						
//						str = move_split_opt(nodeInOrder[idx]);
//						if (!str.equals("") && !str.equals("NaN,NaN,NaN,NaN")) {
//							strs = str.split(",");
//							total += Double.parseDouble(strs[0]);
//							sum1 += Double.parseDouble(strs[1]);
//							sum2 += Double.parseDouble(strs[2]);
//							sum3 += Double.parseDouble(strs[3]);
//						}

//						str = move_split_opt(nodeInOrder[idx]);
//						strs = str.split(",");
//						step1Time += Double.parseDouble(strs[0]);
//						step2Time += Double.parseDouble(strs[1]);
//						step21Time += Double.parseDouble(strs[2]);
//						step22Time += Double.parseDouble(strs[3]);
//						step23Time += Double.parseDouble(strs[4]);
//						step3Time += Double.parseDouble(strs[5]);

//						String str1 = move_split(nodeInOrder[idx]);
//						String str2 = move_split_prune123(nodeInOrder[idx]);
//						if (str1.equals("") && str2.equals("")) continue;
//						
//						if (str1.equals("") && !str2.equals("")) {
//							System.out.println("wrong3 str1 " + str1 + " str2 " + str2);
//							continue;
//						}
//						
//						if (!str1.equals("") && str2.equals("")) {
//							System.out.println("wrong4 str1 " + str1 + " str2 " + str2);
//							continue;
//						}
////						System.out.println(str1 + " " + str2);
//						
//						String[] str1s = str1.split(",");
//						String[] str2s = str2.split(",");
//						
//						double deltaQ1 = Double.parseDouble(str1s[0]);
//						double deltaQ2 = Double.parseDouble(str2s[0]);
//						if (Math.abs(deltaQ1 - deltaQ2) >= eps) {
//							System.out.println("wrong1 str1 " + str1 + " str2 " + str2);
//						}
//						
//						int c1 = Integer.parseInt(str1s[1]);
//						int c2 = Integer.parseInt(str2s[1]);
//						if (c1 != c2) {
//							System.out.println("wrong2 str1 " + str1 + " str2 " + str2);
//						}
					}

//					System.out.println(String.format("%.4f", step1Time) + " " + String.format("%.4f", step2Time) + " "
//							+ String.format("%.4f", step21Time) + " " + String.format("%.4f", step22Time) + " "
//							+ String.format("%.4f", step23Time) + " " + String.format("%.4f", step3Time));

//					System.out.println("total " + total + 
//							" sum1 " + sum1);

//					System.out.println("pass " + count_move + " increasePerPass " + increasePerPass);
					increase_total += increasePerPass;
					if (increasePerPass > eps) hasIncrease = true;
					if (round > n_aggregations) break;
				}
			}
				break;
			}
//			moveTime += (System.currentTimeMillis() - timer) / Constant.RUNNING_TIME_UNIT;

//			System.out.println("increase_total " + increase_total);
			if (increase_total <= eps) break;

//			timer = System.currentTimeMillis();
			rebuildGraph(firstReConstruct);
			constructADJ(false);
//			constructADJ_largeGraph(false);
			if (useSketch) {
				NS_refresh();
				ES_refresh();
			}
			firstReConstruct = false;
//			rebuildTime += (System.currentTimeMillis() - timer) / Constant.RUNNING_TIME_UNIT;

			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);

			/////////////////////////////////////////////////////////////////////

			timeUpdateGamma = calculateGamma();
			if (timeUpdateGamma == -999) break;
			extractClusterTime += (timeUpdateGamma * Constant.RUNNING_TIME_UNIT);

//			matchingTime += timeUpdateGamma;

			/////////////////////////////////////////////////////////////////////

			if (count_aggregations > n_aggregations) break;

//			timer = System.currentTimeMillis();

			startTime = System.currentTimeMillis();
			nodeInOrder = new int[n];
			for (int i = 0; i < n; i++) {
				nodeInOrder[i] = i;
				clusters[i] = i;
				cluster_volume[i] = node_degrees[i];
			}
			endTime = System.currentTimeMillis();
			extractClusterTime += (endTime - startTime);

//			otherTime += (System.currentTimeMillis() - timer) / Constant.RUNNING_TIME_UNIT;

		} while (true);

//		System.out.println("running time " + String.format("%.8f", (extractClusterTime / Constant.RUNNING_TIME_UNIT)));

//		System.out.println("randomTime " + String.format("%.4f", randomTime) + " moveTime " + String.format("%.4f", moveTime)
//						+ " rebuildTime " + String.format("%.4f", rebuildTime) + " matchingTime "
//						+ String.format("%.4f", matchingTime) + " otherTime " + String.format("%.4f", otherTime));

		if (save) {
			saveClusters();
			saveModularity();
		}

		return extractClusterTime + "";
	}

	public double calculateGamma() {
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
		double extractClusterTime = (System.currentTimeMillis() - time) / Constant.RUNNING_TIME_UNIT;
//		System.out.println("calculateGamma " + gamma);
		
		this.gamma = gamma;
		return extractClusterTime;
	}
	
	public void saveClusters() {

		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_cluster_pic_" + this.moveStrategy
					+ "_ordered_" + "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + ratio + "_trial_"
					+ (this.trial) + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_cluster_pic_" + this.moveStrategy
					+ "_ordered_" + "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + ratio
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

	public void saveModularity() {
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/modularity_Louvain_" + this.moveStrategy
					+ "_ordered_" + "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + this.ratio
					+ "_trial_" + this.trial + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/louvain/modularity_Louvain_" + this.moveStrategy
					+ "_ordered_" + "" + this.toHigherOrder + "_order_" + this.ordering + "_ratio_" + this.ratio
					+ "_connect_trial_" + this.trial + ".txt";
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

	private void NS_construct() {
		lastChangeTimes = new int[n];
	}

	private void NS_refresh() {
		for (int i = 0; i < n; i ++) {
			lastChangeTimes[i] = 0;
		}
	}

	private void ES_construct() {

//		int firstIdx_EINC, secondIdx_EINC;

//		// cut off degree
//		int[] baseArr = new int[m];
//		int[] followArr = new int[m];
//		
//		double avgCardinality = 0;
//		for (int edge = 0; edge < m; edge++) {
//			
//			// for each node in the edge
//			firstIdx_EINC = EINC_head[edge];
//			if (edge == m - 1) secondIdx_EINC = EINC_nID_length;
//			else secondIdx_EINC = EINC_head[edge + 1];
//			
//			avgCardinality += (secondIdx_EINC - firstIdx_EINC);
//			
//			baseArr[edge] = secondIdx_EINC - firstIdx_EINC;
//			followArr[edge] = edge;
//		}
//		Dsorter.sort(baseArr, followArr);
//		
//		System.out.println("cardinality min " + baseArr[0] + " max " + baseArr[baseArr.length - 1] + 
//				" avg " + String.format("%.4f", (avgCardinality / (double) m)));

		//////////////////////////////////////////////////////////

//		maxDominNodeNum = (int) Math.floor(1.0 / ratio);
		maxFractionsLazy = new double[m];
		maxFractionChanges = new double[m];
//		dominNodes = new int[m];

//		System.out.println("maxDominNodeNum " + maxDominNodeNum);

		for (int edge = 0; edge < m; edge++) {
			maxFractionsLazy[edge] = EINC_weight[EINC_head[edge]];
//			dominNodes[edge] = -1;
		}
	}

	private void ES_refresh() {

		int firstIdx_EINC, secondIdx_EINC;

//		// cut off degree
//		int[] baseArr = new int[m];
//		int[] followArr = new int[m];
//		int firstIdx_EINC, secondIdx_EINC;
//		double avgCardinality = 0;
//		for (int edge = 0; edge < m; edge++) {
//			
//			// for each node in the edge
//			firstIdx_EINC = EINC_head[edge];
//			if (edge == m - 1) secondIdx_EINC = EINC_nID_length;
//			else secondIdx_EINC = EINC_head[edge + 1];
//			
//			avgCardinality += (secondIdx_EINC - firstIdx_EINC);
//			
//			baseArr[edge] = secondIdx_EINC - firstIdx_EINC;
//			followArr[edge] = edge;
//		}
//		Dsorter.sort(baseArr, followArr);
//		
//		System.out.println("cardinality min " + baseArr[0] + " max " + baseArr[baseArr.length - 1] + 
//				" avg " + String.format("%.4f", (avgCardinality / (double) m)));

		//////////////////////////////////////////////////////////
		double maxFraction;
		int maxIdx;
		for (int edge = 0; edge < m; edge++) {

			// for each node in the edge
			firstIdx_EINC = EINC_head[edge];
			if (edge == endM) secondIdx_EINC = EINC_nID_length;
			else secondIdx_EINC = EINC_head[edge + 1];

			maxFraction = -1;
			maxIdx = -1;
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				if (maxFraction < EINC_weight[i]) {
					maxFraction = EINC_weight[i];
					maxIdx = i;
				}
			}

			maxFractionsLazy[edge] = maxFraction;
			maxFractionChanges[edge] = 0;
//			if (maxRatio >= ratio) {
//				dominNodes[edge] = EINC_nID[maxIdx];
//			} else {
//				dominNodes[edge] = -1;
//			}
		}
	}

	public static void main(String arg[]) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		Hypergraph.loadGraph();
		
		int trials = 3;
		boolean toHigherOrder = false;
		String ordering = "randomOrder";
		String moveStrategy = "move";
		double ratio = 0.6;
		boolean save = false;

//		LouvainHelper LouvainHelper = new LouvainHelper();
//		System.out.println("ratio " + LouvainHelper.findRatioParallel(toHigherOrder, ordering, moveStrategy, 14));

		double avgRunningTime = 0;
		double avgAllMemoryUse = 0;
		String[] strs;
		for (int trial = 0; trial < trials; trial++) {
			System.out.println("trial " + trial);
			PIC_prune12_speed l = new PIC_prune12_speed(trial, toHigherOrder, ordering,
					ratio, save);
			l.initPart1();

			Hypergraph.garbbageCollector.gc();
			long startMem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			l.initPart2();

			strs = l.louvain(moveStrategy).split(",");
			long endMem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			long memoryUse = endMem - startMem;

			avgRunningTime += Double.parseDouble(strs[0]);
			avgAllMemoryUse += memoryUse;
		}
		avgRunningTime /= trials;
		avgAllMemoryUse /= trials;

		System.out.println("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)));
		System.out.println("AllMemoryUse," + String.format("%.8f", (avgAllMemoryUse / Constant.MEMORY_UNIT)));
	}
}
