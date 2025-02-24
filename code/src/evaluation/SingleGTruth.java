package evaluation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.SingleMergeSort;

public class SingleGTruth {
	
	String fileInput_groundTruth;
	String fileInput_discover_pre;
	String method;
	boolean isBaseline;
	
	int minClusterSize = 0;	// regard the cluster with size less than or equals to minClusterSize as outlier, skip
	
	// variables for hypergraph
	double gamma;
	double ratio;
	double vol_H;
	double[] node_weights;
	
	ArrayList<ArrayList<Integer>> allComponents;
	int[] component;
	double[] component_weight;
	
	// Ground Truth Community
	// [cluster]
	int[][] clusters_Ground;
	int[] vertex2Cluster_Ground;
	int classifiedNums_Ground;
	double minClusterSize_Ground;
	double maxClusterSize_Ground;
	double avgClusterSize_Ground;
	double midClusterSize_Ground;
	
	// Discovered Community
	// [trial][cluster]
	List<Integer>[][] clusters_Discover;	// cluster to nodes
	// [trial]
	int[][] vertex2Cluster_Discover;
	int[] classifiedNums_Discover;
	double[] minClusterSizes_Discover;
	double[] maxClusterSizes_Discover;
	double[] avgClusterSizes_Discover;
	double[] midClusterSizes_Discover;
	
	// calculate the metric scores in parallel
	LoadDiscover LoadDiscover;
	CountDownLatch latch;
	
	public boolean DiscoverLoaded;
	public boolean ComponentLoaded;
	
	/////////////////////////////////// for PRF_Prime ///////////////////////////////////
	// [trial]
	double[][] ConfusionMatrix;	//[TP, FN, FP, TN]
	double[] precisions_Prime;	// precision
	double[] recalls_Prime;	// recall
	double[] f1s_Prime;	// f1-measure
	double[] aris_Prime;	// adjusted rand index
	double[] jcc_Prime;	// jaccard index
	double[] balRi_Prime;	// balanced accuracy
	PRF_Prime PRF_Prime;
	
	double[] avgTPFNs;
	double[] avgTPFPs;
	
	/////////////////////////////////// for NMI ///////////////////////////////////
	// [trials]
	double[] NMIs;	// nmi
	NMI NMI;
	
	/////////////////////////////////// for Purity ///////////////////////////////////
	// [trials]
	double[] Purities;	// purity
	Purity Purity;
	
	/////////////////////////////////// for ARIs_Prime ///////////////////////////////////
	// [trials]
	double[] ARIs_Prime;	// adjusted rand index
	ARI_Prime ARI_Prime;
	double[] avgDiscoverCoverCnts;
	double[] avgGroundCoverCnts;
	double[] avgNoCoverCnts;

	// for PIC
	public SingleGTruth(boolean isBaseline, int minClusterSize, String ordering, double ratio) throws IOException, InterruptedException {
		
		// to set the file paths linking to the ground truth clustering and discovered clusterings, then
		// load the ground truth clustering
		
		// file paths to ground truth clustering
		// Note that here, the value of global variable FilePath_Mon.filePathPre has been modified
		// based on the user input using the command processor in exp/Panel.java
		if (!Constant.CONNECTED) fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		else fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		this.isBaseline = isBaseline;
		
		// file path to discovered clustering
		fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/pic/";
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		this.ratio = ratio;
		
		if (!Constant.CONNECTED) method = "pic_move_ordered_false_order_" + ordering + "_ratio_" + ratio;
		else method = "pic_move_ordered_false_order_" + ordering + "_ratio_" + ratio + "_connect";
		
		// load the ground truth clustering
		loadGroundTruth();
		
		DiscoverLoaded = false;
		ComponentLoaded = false;
	}
	
	public void loadGroundTruth() throws IOException {
		
		// a program to load the ground truth clustering
		
		ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();

		ArrayList<Integer> nodes;
		String[] strs;
		try (FileReader reader = new FileReader(fileInput_groundTruth);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line;
			while ((line = bufferedReader.readLine()) != null) {
				
				// process each line
				strs = line.split("\t");
				if (strs.length <= minClusterSize) continue;
				
				nodes = new ArrayList<Integer>(strs.length);
				for (int i = 0; i < strs.length; i++) {
					nodes.add(Integer.parseInt(strs[i]));
				}
				clusters.add(nodes);
			}
		}
		
		avgClusterSize_Ground = 0;
		int[] clusterSizes = new int[clusters.size()];
		
		clusters_Ground = new int[clusters.size()][];
		vertex2Cluster_Ground = new int[Hypergraph.getVertexSize()];
		Arrays.fill(vertex2Cluster_Ground, -1);
		
		int v;
		int clusterSize;
		for (int clusterIdx = 0; clusterIdx < clusters.size(); clusterIdx ++) {
			
			nodes = clusters.get(clusterIdx);
			clusterSize = nodes.size();
			
			clusterSizes[clusterIdx] = clusterSize;
			avgClusterSize_Ground += clusterSize;
			
			clusters_Ground[clusterIdx] = new int[clusterSize];
			
			for (int i = 0; i < clusterSize; i++) {
				v = nodes.get(i);
				clusters_Ground[clusterIdx][i] = v;
				vertex2Cluster_Ground[v] = clusterIdx;
			}
		}
		
		avgClusterSize_Ground /= clusters.size();
		SingleMergeSort SSort = new SingleMergeSort();
		SSort.sort(clusterSizes);
		minClusterSize_Ground = clusterSizes[0];
		maxClusterSize_Ground = clusterSizes[clusterSizes.length - 1];
		midClusterSize_Ground = clusterSizes[0 + (clusterSizes.length-0)/2];
		
		classifiedNums_Ground = 0;
		for (int i = 0; i < vertex2Cluster_Ground.length; i++) {
			if (vertex2Cluster_Ground[i] != -1) classifiedNums_Ground++;
		}
		
		System.out.println("ground truth cluster # " + clusters_Ground.length + " min " + minClusterSize_Ground + 
				" max " + maxClusterSize_Ground + " mid " + midClusterSize_Ground + " avg " + avgClusterSize_Ground);
	} 

	private class LoadDiscover extends Thread{
		private int trial;
		
		LoadDiscover(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}
		
		public void run() {
			
			ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();
			
			try {
				int classifiedNum= 0;
				
				String line;
				String[] strs;
				int v;
				ArrayList<Integer> nodes;
				try (FileReader reader = new FileReader(fileInput_discover_pre + "node_cluster_" + method + "_trial_" + trial + ".txt");
						BufferedReader bufferedReader = new BufferedReader((reader))) {
					while ((line = bufferedReader.readLine()) != null) {
						
						// process each line
						strs = line.split("\t");
						
						if (strs.length <= minClusterSize) continue;
						
						nodes = new ArrayList<Integer>(strs.length);
						
						for (int i = 0; i < strs.length; i++) {
							v = Integer.parseInt(strs[i]);
							nodes.add(v);
							classifiedNum++;
						}
						
						clusters.add(nodes);
					}
				}

				int size;
				double avgSize = 0;
				int[] sizes = new int[clusters.size()];
				
				clusters_Discover[trial] = new ArrayList[clusters.size()];
				Arrays.fill(vertex2Cluster_Discover[trial], -1);
				
				for (int clusterID = 0; clusterID < clusters.size(); clusterID++) {
					
					size = clusters.get(clusterID).size();
					avgSize += size;
					sizes[clusterID] = size;
					
					clusters_Discover[trial][clusterID] = clusters.get(clusterID);
					
					for (int node : clusters.get(clusterID)) vertex2Cluster_Discover[trial][node] = clusterID;
				}
						 
				avgSize /= clusters.size();
				SingleMergeSort SSort = new SingleMergeSort();
				SSort.sort(sizes);
				
				classifiedNums_Discover[trial] = classifiedNum;
				minClusterSizes_Discover[trial] = sizes[0];
				maxClusterSizes_Discover[trial] = sizes[sizes.length - 1];
				midClusterSizes_Discover[trial] = sizes[(sizes.length)/2];
				avgClusterSizes_Discover[trial] = avgSize;

				System.out.println("trial " + trial + " number of discovered cluster " + clusters_Discover[trial].length + 
						" min " + minClusterSizes_Discover[trial] + " max " + maxClusterSizes_Discover[trial] + 
						" mid " + midClusterSizes_Discover[trial] + " avg " + avgClusterSizes_Discover[trial]);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	public void loadDiscover(int trials)
			throws IOException, InterruptedException {
		
		// a program to load the discovered clustering
		
		// initialize the variables for loading the discovered clusters
		clusters_Discover = new ArrayList[trials][];
		vertex2Cluster_Discover = new int[trials][Hypergraph.getVertexSize()];
		classifiedNums_Discover = new int[trials];
		minClusterSizes_Discover = new double[trials];
		maxClusterSizes_Discover = new double[trials];
		midClusterSizes_Discover = new double[trials];
		avgClusterSizes_Discover = new double[trials];
		
		// load the discovered clusters in parallel
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			LoadDiscover = new LoadDiscover(trial);
		}
		latch.await();
		
		// discovered clusters have been loaded
		DiscoverLoaded = true;
	}

	public void loadComponent() throws IOException, InterruptedException {
		
		loadAllComponent();
		
		int n = Hypergraph.getVertexSize();
		component = new int[n];
		component_weight = new double[allComponents.size()];
		
		ArrayList<Integer> componentArr = new ArrayList<Integer>();
		// for each component
		for (int i = 0; i < allComponents.size(); i ++) {
			componentArr = allComponents.get(i);
			for (int node : componentArr) {
				component[node] = i;
				component_weight[i] += (Hypergraph.getSecondIdx_INC(node) - Hypergraph.INC_head[node]);
			}
		}
		
		/////////////////////////////////////////////////////////////////////////////////
		
		ComponentLoaded = true;
	}
	
	private void loadAllComponent() throws FileNotFoundException, IOException {
		
		if (!Constant.CONNECTED) {
			allComponents = new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer> component;
			try (FileReader reader = new FileReader(FilePath_Mon.filePathPre + "/allComponent.txt");
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				String line;
				String[] strs;
				while ((line = bufferedReader.readLine()) != null) {
					strs = line.split("\t");
					component = new ArrayList<Integer>(strs.length);
					
					for (int i = 0; i < strs.length; i++) {
						component.add(Integer.parseInt(strs[i]));
					}
					
					allComponents.add(component);
				}
			}
			
		} else {
			allComponents = new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer> component = new ArrayList<Integer>();
			for (int i = 0; i < Hypergraph.getVertexSize(); i++) {
				component.add(i);
			}
			allComponents.add(component);
		}
	}
		
	private class PRF_Prime extends Thread {
		private int trial;

		// reference: Metrics for evaluating 3D medical image segmentation: analysis, selection, and tool
		// refernece: https://hal.science/hal-00611319/document
		
		PRF_Prime(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			
			int groundTruthNum = clusters_Ground.length;
			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			
			HashSet<Integer>[] discoverIdxs = new HashSet[allComponents.size()];
			HashSet<Integer>[] groundIdxs = new HashSet[allComponents.size()];
			double[] N_comps = new double[allComponents.size()];
			double N_comp;
			for (int compID = 0; compID < allComponents.size(); compID++) {
				N_comp = 0;
				discoverIdxs[compID] = new HashSet<Integer>();
				groundIdxs[compID] = new HashSet<Integer>();
				for (int node : allComponents.get(compID)) {
					if (vertex2Cluster_Ground[node] != -1 && vertex2Cluster_Discover[trial][node] != -1) {
						N_comp++;
						discoverIdxs[compID].add(vertex2Cluster_Discover[trial][node]);
						groundIdxs[compID].add(vertex2Cluster_Ground[node]);
					}
				}
				N_comps[compID] = N_comp;
			}
			
			double avgLocTPFN = 0;
			double avgLocTPFP = 0;
			
			Iterator<Integer> itr_discover;
			Iterator<Integer> itr_ground;
			int[] counterArr;
			double N, TP, FP, FN, TN, TPFN, TPFP, FPTN;
			int clusterSize_Discover, clusterSize_Ground, intersect, discoverIdx, groundTruthIdx;
			double cnt = 0;
			for (int compID = 0; compID < allComponents.size(); compID++) {
						
				N_comp = N_comps[compID];
				N = binomialCoe(N_comp);
				
				if (N_comp == 0) continue;
						
				TPFP = 0;
				TP = 0;
				FP = 0;
				itr_discover = discoverIdxs[compID].iterator();
				int N_discover = 0;
				while (itr_discover.hasNext()) {
					discoverIdx = itr_discover.next();
					
					// counterArr[i] is the intersection of c1 with g_i
					counterArr = new int[groundTruthNum];
					
					clusterSize_Discover = 0;
					// for each element in c1
					for (int vID : clusters_Discover_pointer[discoverIdx]) {
						if (vertex2Cluster_Ground[vID] != -1) {
							clusterSize_Discover++;
							counterArr[vertex2Cluster_Ground[vID]]++;
						}
					}
					N_discover += clusterSize_Discover;
					
					TPFP += binomialCoe(clusterSize_Discover);
					
					itr_ground = groundIdxs[compID].iterator();
					while (itr_ground.hasNext()) {
						groundTruthIdx = itr_ground.next();
						intersect = counterArr[groundTruthIdx];
						TP += binomialCoe(intersect);
					}
				}
				
				FP = TPFP - TP;
				
				TPFN = 0;
				itr_ground = groundIdxs[compID].iterator();
				int N_ground = 0;
				while (itr_ground.hasNext()) {
					groundTruthIdx = itr_ground.next();
					clusterSize_Ground = 0;
					for (int vID : clusters_Ground[groundTruthIdx]) {
						if (vertex2Cluster_Discover[trial][vID] != -1) {
							clusterSize_Ground++;
						}
					}
//					System.out.println("ground " + groundTruthIdx + " " + clusterSize_Ground + " " + binomialCoe(clusterSize_Ground));
					TPFN += binomialCoe(clusterSize_Ground);
					N_ground += clusterSize_Ground;
				}
				
				FN = TPFN - TP;
				
				FPTN = N - TPFN;
				
				TN = FPTN - FP;
				
				if (TPFP == 0 || TPFN == 0) continue;
				
				double precision = TP / (TP + FP);
				
				double recall = TP / (TP + FN);
				
				double f1 = 0;
				if ((precision + recall) == 0) f1 = 0;
				else f1 = (2 * precision * recall) / (precision + recall);
				
				double ri = (TP + TN) / (TP + FP + FN + TN);
				
				double ari = 2 * ((TP * TN) - (FP * FN));
				ari = ari / (((TP + FN) * (FN + TN)) + ((TP + FP) * (FP + TN)));
				
				double jcc = TP / (TP + FP + FN);
				
				double error = (FP + FN) / (TP + FP + FN + TN);
				
				double fpr = FP / (FP + TN);
				
				double sensitivity = TP / (TP + FN);
				double specificity = TN / (TN + FP);
				double balri = (sensitivity + specificity) / 2.0;
//				double balri = 1.0 / (0.5 * ((1.0 / sensitivity) + (1.0 / specificity)));
				
				avgLocTPFN += TPFN;
				avgLocTPFP += TPFP;
				
				ConfusionMatrix[trial][0] += TP;
				ConfusionMatrix[trial][1] += FN;
				ConfusionMatrix[trial][2] += FP;
				ConfusionMatrix[trial][3] += TN;
				precisions_Prime[trial] += precision;
				recalls_Prime[trial] += recall;
				f1s_Prime[trial] += f1;
				aris_Prime[trial] += ari;
				jcc_Prime[trial] += jcc;
				balRi_Prime[trial] += balri;
				cnt++;
			}
			
			avgLocTPFN /= cnt;
			avgLocTPFP /= cnt;
			
			avgTPFNs[trial] = avgLocTPFN;
			avgTPFPs[trial] = avgLocTPFP;
			
			ConfusionMatrix[trial][0] /= cnt;
			ConfusionMatrix[trial][1] /= cnt;
			ConfusionMatrix[trial][2] /= cnt;
			ConfusionMatrix[trial][3] /= cnt;
			precisions_Prime[trial] /= cnt;
			recalls_Prime[trial] /= cnt;
			f1s_Prime[trial] /= cnt;
			aris_Prime[trial] /= cnt;
			jcc_Prime[trial] /= cnt;
			balRi_Prime[trial] /= cnt;

			latch.countDown();
		}
		
		private double binomialCoe(double m) {
			double m_2;
			
			if (m == 0 || m == 1) { 
				m_2 = 0;
			} else {
				m_2 = m * (m - 1) / 2;
			}

			return m_2;
		}
	}

	public void PRF_Prime(int trials)
			throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials);
		
		if (!ComponentLoaded) loadComponent();
		
		avgTPFNs = new double[trials];
		avgTPFPs = new double[trials];
		
		ConfusionMatrix = new double[trials][4];
		precisions_Prime = new double[trials];
		recalls_Prime = new double[trials];
		f1s_Prime = new double[trials];
		aris_Prime = new double[trials];
		jcc_Prime = new double[trials];
		balRi_Prime = new double[trials];
		
		// for discover clusters
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			PRF_Prime = new PRF_Prime(trial);
		}
		latch.await();
		
		double avgTPFN = 0;
		double avgTPFP = 0;
		double cntTPFN = 0;
		double cntTPFP = 0;
		for (int trial = 0; trial < trials; trial++) {
			if (!Double.isNaN(avgTPFNs[trial])) {
				avgTPFN += avgTPFNs[trial];
				cntTPFN++;
			}
			
			if (!Double.isNaN(avgTPFPs[trial])) {
				avgTPFP += avgTPFPs[trial];
				cntTPFP++;
			}
		}
		
		avgTPFN /= cntTPFN;
		avgTPFP /= cntTPFP;
		
		System.out.println(avgTPFN + " " + avgTPFP);
		
		try {
			// for discover clusters
			FileWriter fw_user = null;
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/PRF_Prime_" + method + ".txt",true);
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",precisions_Prime[trial]) + ","
						+ String.format("%.4f",recalls_Prime[trial]) + "," + String.format("%.4f",f1s_Prime[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/PRF_Prime.txt",true);
			double minDiscoverSize = 0;
			double maxDiscoverSize = 0;
			double midDiscoverSize = 0;
			double avgDiscoverSize = 0;
			
			double avgDiscoverNum = 0;
			double avg_TP = 0;
			double avg_FN = 0;
			double avg_FP = 0;
			double avg_TN = 0;
			double avgPrecision = 0;
			double avgRecall = 0;
			double avgF1 = 0;
			double avgARI = 0;
			double avgJCC = 0;
			double avgBALRI = 0;
			
			int cnt_TP = 0;
			int cnt_FN = 0;
			int cnt_FP = 0;
			int cnt_TN = 0;
			int cnt_precision = 0;
			int cnt_recall = 0;
			int cnt_f1 = 0;
			int cnt_ari = 0;
			int cnt_jcc = 0;
			int cnt_balri = 0;
			
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				minDiscoverSize += minClusterSizes_Discover[trial];
				maxDiscoverSize += maxClusterSizes_Discover[trial];
				midDiscoverSize += midClusterSizes_Discover[trial];
				avgDiscoverSize += avgClusterSizes_Discover[trial];
				
				avgDiscoverNum += clusters_Discover[trial].length;
				
				if (!Double.isNaN(ConfusionMatrix[trial][0])) {
					avg_TP += ConfusionMatrix[trial][0];
					cnt_TP++;
				}
				
				if (!Double.isNaN(ConfusionMatrix[trial][1])) {
					avg_FN += ConfusionMatrix[trial][1];
					cnt_FN++;
				}
				
				if (!Double.isNaN(ConfusionMatrix[trial][2])) {
					avg_FP += ConfusionMatrix[trial][2];
					cnt_FP++;
				}
				
				if (!Double.isNaN(ConfusionMatrix[trial][3])) {
					avg_TN += ConfusionMatrix[trial][3];
					cnt_TN++;
				}
				
				if (!Double.isNaN(precisions_Prime[trial])) {
					avgPrecision += precisions_Prime[trial];
					cnt_precision++;
				}
				
				if (!Double.isNaN(recalls_Prime[trial])) {
					avgRecall += recalls_Prime[trial];
					cnt_recall++;
				}
				
				if (!Double.isNaN(f1s_Prime[trial])) {
					avgF1 += f1s_Prime[trial];
					cnt_f1++;
				}
				
				if (!Double.isNaN(aris_Prime[trial])) {
					avgARI += aris_Prime[trial];
					cnt_ari++;
				}
				
				if (!Double.isNaN(jcc_Prime[trial])) {
					avgJCC += jcc_Prime[trial];
					cnt_jcc++;
				}
				
				if (!Double.isNaN(balRi_Prime[trial])) {
					avgBALRI += balRi_Prime[trial];
					cnt_balri++;
				}
			}
			minDiscoverSize /= trials;
			maxDiscoverSize /= trials;
			midDiscoverSize /= trials;
			avgDiscoverSize /= trials;
			
			avgDiscoverNum /= trials;
			avg_TP /= cnt_TP;
			avg_FN /= cnt_FN;
			avg_FP /= cnt_FP;
			avg_TN /= cnt_TN;
			avgPrecision /= cnt_precision;
			avgRecall /= cnt_recall;
			avgF1 /= cnt_f1;
			avgARI /= cnt_ari;
			avgJCC /= cnt_jcc;
			avgBALRI /= cnt_balri;
			
			fw_user.write("PRF_Prime_" + method + "\n");
			fw_user.write(String.format("%.4f",avgF1) + "," + String.format("%.4f",avgPrecision) + "," + String.format("%.4f",avgRecall) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Hypergraph.garbbageCollector.gc();
	}

	private class NMI extends Thread {
		private int trial;

		NMI(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			
			// https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
			// http://www.vldb.org/pvldb/vol2/vldb09-404.pdf
			
			int groundTruthNum = clusters_Ground.length;
			int discoverNum = clusters_Discover[trial].length;

			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];

			// calculate NMI
			double N = 0.0;
			int[] clusterSize_Ground = new int[groundTruthNum];
			for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
				for (int vID : clusters_Ground[groundTruthIdx]) {
					if (vertex2Cluster_Discover[trial][vID] != -1) {
						N++;
						clusterSize_Ground[groundTruthIdx]++;
					}
				}
			}
			
			double MI = 0.0;
			
			// Numerator
			// for each c1
			double n_c;  // size of c
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				
				// counterArr[i] is the intersection of c1 with each g_i,
				double[] counterArr = new double[groundTruthNum];

				// for each element in c1
				n_c = 0.0;
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						n_c++;
						counterArr[vertex2Cluster_Ground[vID]]++;
					}
				}

				double n_g, intersect;
				for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
					n_g = clusterSize_Ground[groundTruthIdx];
					intersect = counterArr[groundTruthIdx];

					if (intersect == 0.0) {
						MI += 0.0;
//						System.out.println(0.0);
					} else {
						MI += ((intersect / N) * Math.log((N * intersect) / (n_c * n_g)));
//						if (((intersect / N) * Math.log((N * intersect) / (n_c * n_g))) > 0) {
//							System.out.println(((intersect / N) * Math.log((N * intersect) / (n_c * n_g))));
//							System.out.println("N " + N + " intersect " + intersect + " N * intersect " + (N * intersect) + " n_c " + n_c + " n_g " + n_g + " n_c * n_g " + n_c * n_g);
//						}
					}
				}
			}

//			System.out.println("MI " + MI);
					
			double entrotyC = 0.0;
			// for each ci
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				n_c = 0.0;
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						n_c++;
					}
				}
				
				if (n_c == 0.0) {
					entrotyC += 0.0;
				} else {
					entrotyC += ((n_c / N) * Math.log(n_c / N));
				}
			}
			entrotyC = -1.0 * entrotyC;
			
//			System.out.println("entrotyC " + entrotyC);
			
			double entrotyG = 0.0;
			double n_g;
			// for each gj
			for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
				n_g = clusterSize_Ground[groundTruthIdx];
				
				if (n_g == 0.0) {
					entrotyG += 0.0;
				} else {
					entrotyG += ((n_g / N) * Math.log(n_g / N));
				}
			}
			entrotyG = -1.0 * entrotyG;
			
//			System.out.println("entrotyG " + entrotyG);
			
			double NMI = MI / ((entrotyC + entrotyG) / 2);
//			
//			System.out.println("NMI " + NMI);
 
			NMIs[trial] = NMI;

			latch.countDown();
		}
	}
	
	public void NMI(int trials)
			throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials);
		
		NMIs = new double[trials];
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			NMI = new NMI(trial);
		}

		latch.await();
		
		try {
			// for discover clusters
			FileWriter fw_user = null;
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/NMI_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",NMIs[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/NMI.txt",true);
			double avgDiscover = 0;
			double avgNMI = 0;
			int cnt = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				if (!Double.isNaN(NMIs[trial])) {
					avgNMI += NMIs[trial];
					cnt++;
				}
			}
			avgDiscover /= trials;
			avgNMI /= cnt;
			fw_user.write("NMI_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgNMI) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		NMIs = null;
		
		Hypergraph.garbbageCollector.gc();
	}

	private class Purity extends Thread {
		private int trial;

		Purity(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			// refer to https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
			
			int groundTruthNum = clusters_Ground.length;
			int discoverNum = clusters_Discover[trial].length;

			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			
			ArrayList<Double> purity_sort = new ArrayList<Double>(discoverNum);
			
			double N = 0;
			double[] counterArr;
			double overlap;
			// for each c1
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				
				// counterArr[i] is the intersection of c1 with g_i,
				counterArr = new double[groundTruthNum];

				// for each element in c1
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						N++;
						counterArr[vertex2Cluster_Ground[vID]]++;
					}
				}
				
				// find the maximum
				overlap = -1;
				for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
					if (overlap < counterArr[groundTruthIdx]) overlap = counterArr[groundTruthIdx];
				}
				
				purity_sort.add(overlap);
			}
 
			double left = (double) 1 / N;
			double right = 0;
			
			if (groundTruthNum <= discoverNum) {
				Collections.sort(purity_sort, Collections.reverseOrder());
				
				for (int i = 0; i < groundTruthNum; i ++) {
					right += purity_sort.get(i);
				}
				
				Purities[trial] = left * right;
			} else {
				for (int i = 0; i < discoverNum; i ++) {
					right += purity_sort.get(i);
				}
				
				Purities[trial] = left * right;
			}

			latch.countDown();
		}
	}
	
	public void Purity(int trials) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials);
		
		Purities = new double[trials];
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			Purity = new Purity(trial);
		}

		latch.await();
		
		try {
			// for discover clusters
			FileWriter fw_user = null;
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/Purity_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",Purities[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/Purity.txt",true);
			double avgDiscover = 0;
			double avgPurity = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				avgPurity += Purities[trial];
			}
			avgDiscover /= trials;
			avgPurity /= trials;
			fw_user.write("Purity_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgPurity) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Purities = null;
		Hypergraph.garbbageCollector.gc();
	}
	
	private class ARI_Prime extends Thread {
		private int trial;

		ARI_Prime(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			// refer to https://link.springer.com/content/pdf/10.1007/BF01908075.pdf
			// or https://dl.acm.org/doi/pdf/10.1145/1281192.1281280
			
			int groundTruthNum = clusters_Ground.length;
			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			double N_graph = Hypergraph.getVertexSize();
			
			HashSet<Integer>[] discoverIdxs = new HashSet[allComponents.size()];
			HashSet<Integer>[] groundIdxs = new HashSet[allComponents.size()];
			double[] N_comps = new double[allComponents.size()];
			double N_comp;
			
			for (int compID = 0; compID < allComponents.size(); compID++) {
				
//				if (!validCompID.contains(compID)) continue;
				
				N_comp = 0;
				discoverIdxs[compID] = new HashSet<Integer>();
				groundIdxs[compID] = new HashSet<Integer>();
				for (int node : allComponents.get(compID)) {
					if (vertex2Cluster_Ground[node] != -1 && vertex2Cluster_Discover[trial][node] != -1) {
						N_comp++;
						discoverIdxs[compID].add(vertex2Cluster_Discover[trial][node]);
						groundIdxs[compID].add(vertex2Cluster_Ground[node]);
					}
				}
				N_comps[compID] = N_comp;
			}
			
			double avgDiscoverCoverCnt = 0;
			double avgGroundCoverCnt = 0;
			double avgNoCoverCnt = 0;
			double coverCnt = 0;
			
			Iterator<Integer> itr_discover;
			Iterator<Integer> itr_ground;
			int[] counterArr;
			double cnt = 0;
			int discoverIdx, groundTruthIdx;
			double N, ni, nij, nj;
			double N_2, ni_2, nj_2, nij_2, denominator, nominator;
			for (int compID = 0; compID < allComponents.size(); compID++) {
				
//				if (!validCompID.contains(compID)) continue;
				
				N = N_comps[compID];
				N_2 = binomialCoe(N);
				
//				System.out.println("N " + N + " N_2 " + N_2);
				
				if (N == 0 || N_2 == 0) continue;
				
				ni_2 = 0;
				itr_discover = discoverIdxs[compID].iterator();
				// for each ci
				while(itr_discover.hasNext()) {
					discoverIdx = itr_discover.next();
					
					ni = 0;
					for (int vID : clusters_Discover_pointer[discoverIdx]) {
						if (vertex2Cluster_Ground[vID] != -1) {
							ni++;
						}
					}
					ni_2 += binomialCoe(ni);
					
//					System.out.println("\tni " + ni + " binomialCoe(ni) " + binomialCoe(ni));
				}
				
//				System.out.println("ni_2 " + ni_2);
				
				nj_2 = 0;
				itr_ground = groundIdxs[compID].iterator();
				// for each gi
				while(itr_ground.hasNext()) {
					groundTruthIdx = itr_ground.next();
					nj = 0;
					for (int vID : clusters_Ground[groundTruthIdx]) {
						if (vertex2Cluster_Discover[trial][vID] != -1) {
							nj++;
						}
					}
					nj_2 += binomialCoe(nj);
					
//					System.out.println("\tnj " + nj + " binomialCoe(nj) " + binomialCoe(nj));
				}
				
//				System.out.println("nj_2 " + nj_2);
				
				nij_2 = 0;
				nij = 0;
				
				itr_discover = discoverIdxs[compID].iterator();
				// for each ci
				while(itr_discover.hasNext()) {
					discoverIdx = itr_discover.next();
					
					// counterArr[i] is the intersection of ci with g_j,
					counterArr = new int[groundTruthNum];

					// for each element in ci
					for (int vID : clusters_Discover_pointer[discoverIdx]) {
						if (vertex2Cluster_Ground[vID] != -1) {
							counterArr[vertex2Cluster_Ground[vID]]++;
						}
					}
					
					for (int j = 0; j < counterArr.length; j++) {
						nij = counterArr[j];
						nij_2 += binomialCoe(nij);
						
//						if (binomialCoe(nij) > 0) System.out.println("\tnij " + nij + " binomialCoe(nij) " + binomialCoe(nij));
					}
				}
				
//				System.out.println("nij_2 " + nij_2);
				
				denominator = ((ni_2 + nj_2) / 2) - ((ni_2 * nj_2) / N_2);
				nominator = nij_2 - ((ni_2 * nj_2) / N_2);
				
				if (denominator == 0) continue;
				
				double ari = nominator / denominator;
				
//				System.out.println("nominator " + nominator + " denominator " + denominator + " ari " + ari);
				
				if (ni_2 == N_2) avgDiscoverCoverCnt ++;
				else if (nj_2 == N_2) avgGroundCoverCnt++;
				else avgNoCoverCnt++;
				coverCnt++;
				
				if (ni_2 == N_2) continue;
				else if (nj_2 == N_2) continue;
				if (ari <= 0) continue;
				
				if (Double.isNaN(ari)) {
					System.out.println(nominator + " " + denominator);
					System.out.println("ni_2 " + ni_2 + " nj_2 " + nj_2 + " nij_2 " + nij_2 + " N_2 " + N_2);
				}
				
				ARIs_Prime[trial] += ari;
				
				cnt++;
			}
			
			if (cnt == 0) ARIs_Prime[trial] = -2;
			else ARIs_Prime[trial] /= (double)cnt;
			
			avgDiscoverCoverCnt /= coverCnt;
			avgGroundCoverCnt /= coverCnt;
			avgNoCoverCnt /= coverCnt;
			
			avgDiscoverCoverCnts[trial] = avgDiscoverCoverCnt;
			avgGroundCoverCnts[trial] = avgGroundCoverCnt;
			avgNoCoverCnts[trial] = avgNoCoverCnt;

			latch.countDown();
		}
		
		private double binomialCoe(double m) {
			double m_2;
			if (m == 0 || m == 1) { 
				m_2 = 0;
			} else {
				m_2 = m * (m - 1) / 2;
			}
			return m_2;
		}
	}
	
	public void ARI_Prime(int trials) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials);
		
		if (!ComponentLoaded) loadComponent();
		
		ARIs_Prime = new double[trials];
		
		avgDiscoverCoverCnts = new double[trials];
		avgGroundCoverCnts = new double[trials];
		avgNoCoverCnts = new double[trials];
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			ARI_Prime = new ARI_Prime(trial);
		}
		latch.await();
		
		double avgDiscoverCoverCnt = 0;
		double avgGroundCoverCnt = 0;
		double avgNoCoverCnt = 0;
		double DiscoverCoverCnt = 0;
		double GroundCoverCnt = 0;
		double NoCoverCnt = 0;
		for (int trial = 0; trial < trials; trial++) {
			if (!Double.isNaN(avgDiscoverCoverCnts[trial])) {
				avgDiscoverCoverCnt += avgDiscoverCoverCnts[trial];
				DiscoverCoverCnt++;
			}
			
			if (!Double.isNaN(avgGroundCoverCnts[trial])) {
				avgGroundCoverCnt += avgGroundCoverCnts[trial];
				GroundCoverCnt++;
			}
			
			if (!Double.isNaN(avgNoCoverCnts[trial])) {
				avgNoCoverCnt += avgNoCoverCnts[trial];
				NoCoverCnt++;
			}
		}
		
		avgDiscoverCoverCnt /= DiscoverCoverCnt;
		avgGroundCoverCnt /= GroundCoverCnt;
		avgNoCoverCnt /= NoCoverCnt;
		
		System.out.println(avgDiscoverCoverCnt + " " + avgGroundCoverCnt + " " + avgNoCoverCnt);
		
		try {
			// for discover clusters
			FileWriter fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/ARI_Prime_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",ARIs_Prime[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/ARI_Prime.txt",true);
			double avgDiscover = 0;
			double avgARI = 0;
			double cnt = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				if (ARIs_Prime[trial] > -2) {
					avgARI += ARIs_Prime[trial];
					cnt++;
				}
			}
			avgDiscover /= trials;
			if (cnt > 0) avgARI /= cnt;
			else avgARI = -2;
			fw_user.write("ARI_Prime_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgARI) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ARIs_Prime = null;
		Hypergraph.garbbageCollector.gc();
	}

}
