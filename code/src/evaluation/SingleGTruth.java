package evaluation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
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
	
	int maxClusterSize;

	// Hypergraph
	int n;
	double gamma;
	double totalEdgeWeight;
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
	
	// for top k
	// [trial][cluster]
	int[][][] clusters_Ground_largestk;
	int[][] vertex2Cluster_Ground_largestk;
	int[] classifiedNums_Ground_largestk;
	// [trial][cluster]
	List<Integer>[][] clusters_Discover_largestk;	// cluster to nodes
	int[][] vertex2Cluster_Discover_largestk;
	int[] classifiedNums_Discover_largestk;
	
	LoadDiscover LoadDiscover;
	LoadDiscover_GroundTruthOnly LoadDiscover_GroundTruthOnly;

	CountDownLatch latch;
	
	public boolean DiscoverLoaded;
	public boolean GraphLoaded;

	/////////////////////////////////// for PRF ///////////////////////////////////
	// [trial]
	double[] precisions;
	double[] recalls;
	double[] f1s;
	private PRF PRF;
	
	/////////////////////////////////// for PRF_Prime ///////////////////////////////////
	// [trial]
	double[] precisions_Prime;
	double[] recalls_Prime;
	double[] f1s_Prime;
	double[] ris_Prime;
	double[] aris_Prime;
	PRF_Prime PRF_Prime;
	
	double[] avgTPFNs;
	double[] avgTPFPs;
	
	/////////////////////////////////// for NMI ///////////////////////////////////
	// [trials]
	double[] NMIs;
	NMI NMI;
	
	/////////////////////////////////// for Purity ///////////////////////////////////
	// [trials]
	double[] Purities;
	Purity Purity;
	
	/////////////////////////////////// for ARIs ///////////////////////////////////
	// [trials]
	double[] ARIs;
	ARI ARI;
	
	/////////////////////////////////// for ARIs_Prime ///////////////////////////////////
	// [trials]
	double[] ARIs_Prime;
	ARI_Prime ARI_Prime;
	double[] avgDiscoverCoverCnts;
	double[] avgGroundCoverCnts;
	double[] avgNoCoverCnts;
	
	///////////////////////////// for HModularity /////////////////////////////////
	// [trials]
	double[] HModularities;
	double HModularities_Ground;
	HModularity HModularity;
	
	///////////////////////////// for HConductance /////////////////////////////////
	// [trials]
	double[] HConductances;
	double HConductances_Ground;
	HConductance HConductance;

	// for BogLouvain, BogCNMOpt, BogCNMRan, BipartiteLouvain, IRMM, HMLL
	public SingleGTruth(boolean isBaseline, String algorithm, int maxClusterSize) throws IOException, InterruptedException {
		
		System.out.println(algorithm);
		
		if (!Constant.CONNECTED) fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		else fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		this.isBaseline = isBaseline;
		
		switch (algorithm) {
			
		case "BogLouvain":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "BogLouvain";
				else method = "BogLouvain_connect";
			}
			
			break;
			
		case "BogCNMOpt":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "BogCNMOpt";
				else method = "BogCNMOpt_connect";
			}
			
			break;
		
		case "BogCNMRan":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "BogCNMRan";
				else method = "BogCNMRan_connect";
			}
			
			break;
		case "BipartiteLouvain":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "BipartiteLouvain";
				else method = "BipartiteLouvain_connect";
			}
			break;
			
		case "IRMM":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "IRMM";
				else method = "IRMM_connect";
			}
			break;
		
		case "HMLL":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "HMLL";
				else method = "HMLL_connect";
			}
			break;
			
		case "GMLL":
			{
				fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
				if (!Constant.CONNECTED) method = "GMLL";
				else method = "GMLL_connect";
			}
			break;
			
		case "HPPR":
		{
			fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/baseline/";
			if (!Constant.CONNECTED) method = "HPPR";
			else method = "HPPR_connect";
		}
			break;
		}
		
		loadGroundTruth(maxClusterSize);
		
		DiscoverLoaded = false;
		GraphLoaded = false;
	}
	
	// for PIC
	public SingleGTruth(boolean isBaseline, String moveStrategy, int strategy, boolean toHigherOrder, String ordering, double ratio) throws IOException, InterruptedException {
		
		System.out.println(moveStrategy);
		
		if (!Constant.CONNECTED) fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		else fileInput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		this.isBaseline = isBaseline;
		
		fileInput_discover_pre = FilePath_Mon.filePathPre + "/clustering/pic/";
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		if (!Constant.CONNECTED) method = "pic_" + moveStrategy + "_ordered_" + toHigherOrder + "_order_" + ordering + "_ratio_" + ratio;
		else method = "pic_" + moveStrategy + "_ordered_" + toHigherOrder + "_order_" + ordering + "_ratio_" + ratio + "_connect";
		
		loadGroundTruth(strategy);
		
		DiscoverLoaded = false;
		GraphLoaded = false;
	}
	
	public void loadGroundTruth(int maxClusterSize) throws IOException {
		
		this.maxClusterSize = maxClusterSize;
		
		ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();

		ArrayList<Integer> nodes;
		String[] strs;
		try (FileReader reader = new FileReader(fileInput_groundTruth);
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
						
						if (strs.length <= maxClusterSize) continue;
						
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
	
	public void loadDiscover(int trials, int maxClusterSize)
			throws IOException, InterruptedException {
		
		this.maxClusterSize = maxClusterSize;
		clusters_Discover = new ArrayList[trials][];
		vertex2Cluster_Discover = new int[trials][Hypergraph.getVertexSize()];
		classifiedNums_Discover = new int[trials];
		minClusterSizes_Discover = new double[trials];
		maxClusterSizes_Discover = new double[trials];
		midClusterSizes_Discover = new double[trials];
		avgClusterSizes_Discover = new double[trials];
		
		clusters_Ground_largestk = new int[trials][][];
		vertex2Cluster_Ground_largestk = new int[trials][];
		classifiedNums_Ground_largestk = new int[trials];
		// [trial][cluster]
		clusters_Discover_largestk = new ArrayList[trials][];;
		vertex2Cluster_Discover_largestk = new int[trials][Hypergraph.getVertexSize()];;
		classifiedNums_Discover_largestk = new int[trials];
		
		latch = new CountDownLatch(trials);
		
		for (int trial = 0; trial < trials; trial++) {
			LoadDiscover = new LoadDiscover(trial);
		}
		
		latch.await();
		
		DiscoverLoaded = true;
	}

	private class LoadDiscover_GroundTruthOnly extends Thread{
		private int trial;
		
		LoadDiscover_GroundTruthOnly(int trial) throws InterruptedException {
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
						
						if (strs.length <= maxClusterSize) continue;
						
						nodes = new ArrayList<Integer>(strs.length);
						
						for (int i = 0; i < strs.length; i++) {
							v = Integer.parseInt(strs[i]);
							if (vertex2Cluster_Ground[v] == -1) continue;
							
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
					
					for (int node : clusters.get(clusterID)) {
						vertex2Cluster_Discover[trial][node] = clusterID;
					}
				}
				
				avgSize /= clusters.size();
				SingleMergeSort SSort = new SingleMergeSort();
				SSort.sort(sizes);
				
				classifiedNums_Discover[trial] = classifiedNum;
				minClusterSizes_Discover[trial] = sizes[0];
				maxClusterSizes_Discover[trial] = sizes[sizes.length - 1];
				midClusterSizes_Discover[trial] = sizes[(sizes.length-0)/2];
				avgClusterSizes_Discover[trial] = avgSize;

//				System.out.println("trial " + trial + " number of discovered cluster " + clusters_Discover[trial].length + " classified node size " + classifiedNum);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	public void loadDiscover_GroundTruthOnly(int trials, int maxClusterSize)
			throws IOException, InterruptedException {
		
		this.maxClusterSize = maxClusterSize;
		
		clusters_Discover = new ArrayList[trials][];
		vertex2Cluster_Discover = new int[trials][Hypergraph.getVertexSize()];
		classifiedNums_Discover = new int[trials];
		
		latch = new CountDownLatch(trials);
		
		for (int trial = 0; trial < trials; trial++) {
			LoadDiscover_GroundTruthOnly = new LoadDiscover_GroundTruthOnly(trial);
		}
		
		latch.await();
		
		DiscoverLoaded = true;
	}

	public void loadgraph() throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
		/////////////////////////////////////////////////////////////////////////////////
		double lambda = -1;
		gamma = -1;
		String[] strs = FilePath_Mon.filePathPre.split("/");
		String dataset = strs[strs.length - 1];
		Path path = Paths.get(System.getProperty("user.dir") + "/python/lambdas.txt");
		Scanner scanner = new Scanner(path.toAbsolutePath());
		
		while(scanner.hasNextLine()){
		    //process each line
		    strs = scanner.nextLine().split("\t");
		    if (strs[0].equals(dataset)) lambda = Double.parseDouble(strs[1]);
		}
		
		scanner.close();
		
		gamma = Math.exp(-1 * lambda);
		System.out.println(dataset + " lambda " + lambda + " gamma " + gamma);
		
		/////////////////////////////////////////////////////////////////////////////////
		
		n = Hypergraph.getVertexSize();
		node_weights = new double[n];
		
		double[] INC_weight = Hypergraph.INC_weight;
		int[] INC_head = Hypergraph.INC_head;
    	double weight;
    	int curID = 0;
    	int first_idx, second_idx;
    	
    	// for each node
    	for (curID = 0; curID < n - 1; curID++) {
    		first_idx = INC_head[curID];
			second_idx = INC_head[curID + 1];
			
			// for each neighbor
			for (int i = first_idx; i < second_idx; i++) {
				weight = INC_weight[i];
				node_weights[curID] += weight;
	        	totalEdgeWeight += weight;
			}
    	}
    	first_idx = INC_head[curID];
		second_idx = INC_weight.length;
		for (int i = first_idx; i < second_idx; i++) {
			weight = INC_weight[i];
			node_weights[curID] += weight;
        	totalEdgeWeight += weight;
		}
		
		/////////////////////////////////////////////////////////////////////////////////
		
		loadAllComponent();
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
		
		GraphLoaded = true;
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
		
	private class PRF extends Thread {
		private int trial;

		PRF(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			int groundTruthNum = clusters_Ground.length;
			int discoverNum = clusters_Discover[trial].length;
			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			
			double matchingNum = discoverNum;
			
			int[] clusterSize_Ground = new int[groundTruthNum];
			for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
				for (int vID : clusters_Ground[groundTruthIdx]) {
					if (vertex2Cluster_Discover[trial][vID] != -1) {
						clusterSize_Ground[groundTruthIdx]++;
					}
				}
			}
			
			ArrayList<Double> precision_sort = new ArrayList<Double>(discoverNum);
			ArrayList<Double> recall_sort = new ArrayList<Double>(discoverNum);
			ArrayList<Double> f1_sort = new ArrayList<Double>(discoverNum);
					
			// for each c1
			boolean hasMatching;
			double clusterSize_Discover;
			double[] counterArr;
			double intersect, tempPrecision, tempRecall, tempF1, maxPrecision, maxRecall, maxF1;
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {

				hasMatching = false;
				clusterSize_Discover = 0;
				
				// counterArr[i] is the intersection of c1 with g_i,
				counterArr = new double[groundTruthNum];

				// for each element in c1
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						clusterSize_Discover++;
						counterArr[vertex2Cluster_Ground[vID]]++;
						hasMatching = true;
					}
				}
				
				if (!hasMatching) {
					matchingNum--;
					continue;
				}
				
				maxPrecision = -1;
				maxRecall = -1;
				maxF1 = -1;
				
				for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
					intersect = counterArr[groundTruthIdx];
					
					tempPrecision = intersect / clusterSize_Discover;
					if (tempPrecision > maxPrecision) maxPrecision = tempPrecision;
					
					tempRecall = intersect / (double) clusterSize_Ground[groundTruthIdx];
					if (tempRecall > maxRecall) maxRecall = tempRecall;
					
					tempF1 = (2 * tempPrecision * tempRecall) / (tempPrecision + tempRecall);
					if (tempF1 > maxF1) maxF1 = tempF1;
				}
				
				precision_sort.add(maxPrecision);
				recall_sort.add(maxRecall);
				f1_sort.add(maxF1);
			}
			
			if (groundTruthNum <= matchingNum) {
				Collections.sort(precision_sort, Collections.reverseOrder());   
				Collections.sort(recall_sort, Collections.reverseOrder());   
				Collections.sort(f1_sort, Collections.reverseOrder());
				
				for (int i = 0; i < groundTruthNum; i ++) {
					precisions[trial] += precision_sort.get(i);
					recalls[trial] += recall_sort.get(i);
					f1s[trial] += f1_sort.get(i);
				}
				
				precisions[trial] /= groundTruthNum;
				recalls[trial] /= groundTruthNum;
				f1s[trial] /= groundTruthNum;
			} else {
				for (int i = 0; i < matchingNum; i ++) {
					precisions[trial] += precision_sort.get(i);
					recalls[trial] += recall_sort.get(i);
					f1s[trial] += f1_sort.get(i);
				}
				
				precisions[trial] /= matchingNum;
				recalls[trial] /= matchingNum;
				f1s[trial] /= matchingNum;
			}

			latch.countDown();
		}
	}

	public void PRF(int trials, int strategy)
			throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		precisions = new double[trials];
		recalls = new double[trials];
		f1s = new double[trials];
		
		// for discover clusters
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			PRF = new PRF(trial);
		}
		latch.await();
		
		try {
			// for discover clusters
			FileWriter fw_user = null;
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/PRF_" + method + ".txt",true);
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",precisions[trial]) + ","
						+ String.format("%.4f",recalls[trial]) + "," + String.format("%.4f",f1s[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/PRF.txt",true);
			double avgDiscoverSize = 0;
			double avgPrecision = 0;
			double avgRecall = 0;
			double avgF1 = 0;
			int cnt_precision = 0;
			int cnt_recall = 0;
			int cnt_f1 = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscoverSize += clusters_Discover[trial].length;
				
				if (!Double.isNaN(precisions[trial])) {
					avgPrecision += precisions[trial];
					cnt_precision++;
				}
				
				if (!Double.isNaN(recalls[trial])) {
					avgRecall += recalls[trial];
					cnt_recall++;
				}
				
				if (!Double.isNaN(f1s[trial])) {
					avgF1 += f1s[trial];
					cnt_f1++;
				}
				
			}
			avgDiscoverSize /= trials;
			avgPrecision /= cnt_precision;
			avgRecall /= cnt_recall;
			avgF1 /= cnt_f1;
			fw_user.write("PRF_" + method + "\n");
			fw_user.write(String.format("%.1f",avgDiscoverSize) + "," + clusters_Ground.length + "," + String.format("%.4f",avgPrecision) + ","
					+ String.format("%.4f",avgRecall) + "," + String.format("%.4f",avgF1) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Hypergraph.garbbageCollector.gc();
	}

	private class PRF_Prime extends Thread {
		private int trial;

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
				
				if (precision > 1 || recall > 1 || f1 > 1 || ri > 1 || precision < 0 || recall < 0 || f1 < 0 || ri < 0) {
					System.out.println("N_comp " + N_comp + " N_discover " + N_discover + " N_ground " + N_ground);
					System.out.println("N " + N + " TP " + TP + " FN " + FN + " FP " + FP + " TN " + TN);
					System.out.println("TPFP " + TPFP + " TPFN " + TPFN + " FPTN " + FPTN);
					System.out.println("precision " + precision + " recall " + recall + " f1 " + f1 + " ri " + ri);
				}
				
				avgLocTPFN += TPFN;
				avgLocTPFP += TPFP;
				
				precisions_Prime[trial] += precision;
				recalls_Prime[trial] += recall;
				f1s_Prime[trial] += f1;
				ris_Prime[trial] += ri;
				aris_Prime[trial] += ari;
				cnt++;
			}
			
			avgLocTPFN /= cnt;
			avgLocTPFP /= cnt;
			
			avgTPFNs[trial] = avgLocTPFN;
			avgTPFPs[trial] = avgLocTPFP;
			
			precisions_Prime[trial] /= cnt;
			recalls_Prime[trial] /= cnt;
			f1s_Prime[trial] /= cnt;
			ris_Prime[trial] /= cnt;
			aris_Prime[trial] /= cnt;
			
			if (precisions_Prime[trial] > 1 || recalls_Prime[trial] > 1 || f1s_Prime[trial] > 1 || ris_Prime[trial] > 1 || aris_Prime[trial] > 1) {
				System.out.println("trial precision " + precisions_Prime[trial] + 
						" recall " + recalls_Prime[trial] + " f1 " + f1s_Prime[trial] + " ri " + ris_Prime[trial] + " ri " + ris_Prime[trial]);
			}

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

	public void PRF_Prime(int trials, int strategy)
			throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		if (!GraphLoaded) {
			loadgraph();
		}
		
		avgTPFNs = new double[trials];
		avgTPFPs = new double[trials];
		
		precisions_Prime = new double[trials];
		recalls_Prime = new double[trials];
		f1s_Prime = new double[trials];
		ris_Prime = new double[trials];
		aris_Prime = new double[trials];
		
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
						+ String.format("%.4f",recalls_Prime[trial]) + "," + String.format("%.4f",f1s_Prime[trial]) + "," + String.format("%.4f",ris_Prime[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/PRF_Prime.txt",true);
			double minDiscoverSize = 0;
			double maxDiscoverSize = 0;
			double midDiscoverSize = 0;
			double avgDiscoverSize = 0;
			
			double avgDiscoverNum = 0;
			double avgPrecision = 0;
			double avgRecall = 0;
			double avgF1 = 0;
			double avgRI = 0;
			double avgARI = 0;
			int cnt_precision = 0;
			int cnt_recall = 0;
			int cnt_f1 = 0;
			int cnt_ri = 0;
			int cnt_ari = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				minDiscoverSize += minClusterSizes_Discover[trial];
				maxDiscoverSize += maxClusterSizes_Discover[trial];
				midDiscoverSize += midClusterSizes_Discover[trial];
				avgDiscoverSize += avgClusterSizes_Discover[trial];
				
				avgDiscoverNum += clusters_Discover[trial].length;
				
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
				
				if (!Double.isNaN(ris_Prime[trial])) {
					avgRI += ris_Prime[trial];
					cnt_ri++;
				}
				
				if (!Double.isNaN(aris_Prime[trial])) {
					avgARI += aris_Prime[trial];
					cnt_ari++;
				}
			}
			minDiscoverSize /= trials;
			maxDiscoverSize /= trials;
			midDiscoverSize /= trials;
			avgDiscoverSize /= trials;
			
			avgDiscoverNum /= trials;
			avgPrecision /= cnt_precision;
			avgRecall /= cnt_recall;
			avgF1 /= cnt_f1;
			avgRI /= cnt_ri;
			avgARI /= cnt_ari;
			
			fw_user.write("PRF_Prime_" + method + "\n");
			fw_user.write(minClusterSize_Ground + "," + maxClusterSize_Ground + "," + midClusterSize_Ground+ "," + avgClusterSize_Ground + "\n");
			fw_user.write(String.format("%.1f",minDiscoverSize) + "," + String.format("%.1f",maxDiscoverSize) + 
					"," + String.format("%.1f",midDiscoverSize) + "," + String.format("%.1f",avgDiscoverSize) + "\n");
			fw_user.write(String.format("%.1f",avgDiscoverNum) + "," + clusters_Ground.length + "," + String.format("%.4f",avgPrecision) + 
					"," + String.format("%.4f",avgRecall) +  "," + String.format("%.4f",avgF1) + "," + String.format("%.4f",avgRI) + "," + String.format("%.4f",avgARI) + "\n");
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
	
	public void NMI(int trials, int strategy)
			throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
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
	
	public void Purity(int trials, int strategy) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
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
	
	private class ARI extends Thread {
		private int trial;

		ARI(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}

		public void run() {
			// refer to https://link.springer.com/content/pdf/10.1007/BF01908075.pdf
			// or https://dl.acm.org/doi/pdf/10.1145/1281192.1281280
			
			int groundTruthNum = clusters_Ground.length;
			int discoverNum = clusters_Discover[trial].length;

			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			
			double N = 0;
			int[] clusterSize_Ground = new int[groundTruthNum];
			for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
				for (int vID : clusters_Ground[groundTruthIdx]) {
					if (vertex2Cluster_Discover[trial][vID] != -1) {
						N++;
						clusterSize_Ground[groundTruthIdx]++;
					}
				}
			}
			double N_2 = binomialCoe(N);
			
			double ni_2 = 0;
			double ni;
			// for each ci
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				ni = 0;
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						ni++;
					}
				}
				ni_2 += binomialCoe(ni);
			}
			
			double nj_2 = 0;
			double nj;
			// for each gi
			for (int groundTruthIdx = 0; groundTruthIdx < groundTruthNum; groundTruthIdx++) {
				nj = clusterSize_Ground[groundTruthIdx];
				nj_2 += binomialCoe(nj);
			}
			
			double nij_2 = 0;
			double nij = 0;
			
			double[] counterArr;
			// for each ci
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				
				// counterArr[i] is the intersection of ci with g_j,
				counterArr = new double[groundTruthNum];

				// for each element in ci
				for (int vID : clusters_Discover_pointer[discoverIdx]) {
					if (vertex2Cluster_Ground[vID] != -1) {
						counterArr[vertex2Cluster_Ground[vID]]++;
					}
				}
				
				for (int j = 0; j < counterArr.length; j++) {
					nij = counterArr[j];
					nij_2 += binomialCoe(nij);
				}
			}
 
			double denominator = ((ni_2 + nj_2) / 2) - ((ni_2 * nj_2) / N_2);
			double nominator = nij_2 - ((ni_2 * nj_2) / N_2);
			
			if (N == 0 || N_2 == 0) ARIs[trial] = -2;
			if (denominator == 0) ARIs[trial] = -2;
			
			double ari = nominator / denominator;
			
			if (ari <= 0) ARIs[trial] = -2;
			else ARIs[trial] = ari;

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
	
	public void ARI(int trials, int strategy) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		ARIs = new double[trials];
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			ARI = new ARI(trial);
		}

		latch.await();
		
		try {
			// for discover clusters
			FileWriter fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/ARI_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f",ARIs[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/ARI.txt",true);
			double avgDiscover = 0;
			double avgARI = 0;
			double cnt = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				if (ARIs[trial] > -2) {
					avgARI += ARIs[trial];
					cnt++;
				}
			}
			avgDiscover /= trials;
			if (cnt > 0) avgARI /= cnt;
			else avgARI = -2;
			fw_user.write("ARI_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgARI) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ARIs = null;
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
			
//			double[] baseArr = new double[allComponents.size()];
//			int[] followArr = new int[allComponents.size()];
//			int N_total = 0;
//			for (int compID = 0; compID < allComponents.size(); compID++) {
//				N_comp = 0;
//				discoverIdxs[compID] = new HashSet<Integer>();
//				groundIdxs[compID] = new HashSet<Integer>();
//				for (int node : allComponents.get(compID)) {
//					if (vertex2Cluster_Ground[node] != -1 && vertex2Cluster_Discover[trial][node] != -1) {
//						N_comp++;
//						N_total++;
//						discoverIdxs[compID].add(vertex2Cluster_Discover[trial][node]);
//						groundIdxs[compID].add(vertex2Cluster_Ground[node]);
//					}
//				}
//				baseArr[compID] = N_comp;
//				followArr[compID] = compID;
//			}
//			DoubleMergeSort Dsort = new DoubleMergeSort();
//			Dsort.sort(baseArr, followArr);
//			
//			int clusterNum = (int) Math.min(Math.ceil(Math.log(N_total) / Math.log(2)), followArr.length);
//			
//			HashSet<Integer> validCompID = new HashSet<Integer>(clusterNum);
//			for(int i = followArr.length - 1; i >= (followArr.length - clusterNum); i--) {
//				validCompID.add(followArr[i]);
//			}
			
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
	
	public void ARI_Prime(int trials, int strategy) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		if (!GraphLoaded) {
			loadgraph();
		}
		
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

	private class HModularity extends Thread {
		private int trial;

		HModularity(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}
		
		public void run() {
			// for discovered clusters
			int discoverNum = clusters_Discover[trial].length;
			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			
			// for each c1
			double modularity = 0;
			HashSet<Integer> nodes;
			Iterator<Integer> itr;
			double phi, edgeWeight, vol_c, gamma_c;
			int edgeSize, second_idx_neighbor;
			boolean contain;
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				
				nodes = new HashSet<Integer>(clusters_Discover_pointer[discoverIdx].size());
				for (int node : clusters_Discover_pointer[discoverIdx]) {
					nodes.add(node);
				}
				
				phi = 0;
			    edgeSize = Hypergraph.getEdgeSize();
			    for (int edgeID = 0; edgeID < edgeSize; edgeID++) {
					
			    	contain = true;
			    	
			    	second_idx_neighbor = Hypergraph.getSecondIdx_EINC(edgeID);
					edgeWeight = 0;
					for (int j = Hypergraph.EINC_head[edgeID]; j < second_idx_neighbor; j++) {
						if (!nodes.contains(Hypergraph.EINC_nID[j])) contain = false;
						edgeWeight += Hypergraph.EINC_weight[j];
					}
					
					if (contain) phi += edgeWeight;
			    }
			    
			    vol_c = 0;
			    itr = nodes.iterator();
			    while (itr.hasNext()) {
			    	vol_c += node_weights[itr.next()];
			    }
				
				gamma_c = gamma * vol_c;
				modularity += (((phi + gamma_c) / totalEdgeWeight) - (gamma_c / (totalEdgeWeight - gamma_c)));
			}
					
			HModularities[trial] = modularity;
			
			latch.countDown();
		}
	}
		
	public void HModularity(int trials, int strategy) throws IOException, InterruptedException {
		
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		HModularities = new double[trials];
		
		if (isBaseline) {
			if (!GraphLoaded) loadgraph();
			
			// for discover clusters
			latch = new CountDownLatch(trials);
			for (int trial = 0; trial < trials; trial++) {
				HModularity = new HModularity(trial);
			}
			latch.await();
			
		} else {
			String line;
			for (int trial = 0; trial < trials; trial++) {
				try (FileReader reader = new FileReader(fileInput_discover_pre + "modularity_" + method + "_trial_" + trial + ".txt");
						BufferedReader bufferedReader = new BufferedReader((reader))) {
					while ((line = bufferedReader.readLine()) != null) {
						HModularities[trial] = Double.parseDouble(line.split("\t")[1]);
					}
				}
			}
		}
		
		try {
			// for discover clusters
			FileWriter fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/HModularity_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f", HModularities[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/HModularity.txt",true);
			double avgDiscover = 0;
			double avgHMod = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				avgHMod += HModularities[trial];
			}
			avgDiscover /= trials;
			avgHMod /= trials;
			fw_user.write("HModularity_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgHMod) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		HModularities = null;
		
		Hypergraph.garbbageCollector.gc();
	}
	
	private class HConductance extends Thread {
		private int trial;

		HConductance(int trial) throws InterruptedException {
			this.trial = trial;
			start();
		}
		
		public void run() {
			// for discovered clusters
			int discoverNum = clusters_Discover[trial].length;
			List<Integer>[] clusters_Discover_pointer = clusters_Discover[trial];
			int[] vertex2Cluster_Discover_pointer = vertex2Cluster_Discover[trial];
			
			double[] cond_components = new double[component_weight.length];
			for (int i = 0; i < cond_components.length; i++) cond_components[i] = Constant.large;
			
			// for each c1
			int componentID, first_idx, second_idx, u;
			double c_s, m_s, vol_s, vol_sbar, vol_component, denominator, cond;
			HashSet<Integer> nodeInCluster;
			for (int discoverIdx = 0; discoverIdx < discoverNum; discoverIdx++) {
				
				c_s = 0;
				m_s = 0;
				vol_s = 0;
				vol_sbar = 0;
				vol_component = 0;
				
				componentID = component[clusters_Discover_pointer[discoverIdx].get(0)];
				vol_component = component_weight[componentID];
				
				nodeInCluster = new HashSet<Integer>();
				for (int v : clusters_Discover_pointer[discoverIdx]) {
					nodeInCluster.add(v);
				}
				
				// for each element in c1
				for (int v : clusters_Discover_pointer[discoverIdx]) {
							
					vol_s += (Hypergraph.getSecondIdx_INC(v) - Hypergraph.INC_head[v]);
					
					int edgeID;
					first_idx = Hypergraph.INC_head[v];
					second_idx = Hypergraph.getSecondIdx_INC(v);
					for (int j = first_idx; j < second_idx; j++) {
						edgeID = Hypergraph.INC_eID[j];
						
						int first_idx_edge = Hypergraph.EINC_head[edgeID];
						int second_idx_edge = Hypergraph.getSecondIdx_EINC(edgeID);
						boolean allEdge = true;
						for (int k = first_idx_edge; k < second_idx_edge; k++) {
							u = Hypergraph.EINC_nID[k];
							if (vertex2Cluster_Discover_pointer[u] != -1 && vertex2Cluster_Discover_pointer[u] != discoverIdx) {
								allEdge = false;
							}
						}
						
						if (allEdge) m_s += 1;
						else c_s += 1;
					}
				}
				
				vol_sbar = vol_component - vol_s;
				
				denominator = Math.min(vol_s, vol_sbar);
				
				if (denominator <= 0) {
					continue;
				}
				if (c_s <= 0) {
					continue;
				}
				
				cond = c_s / denominator;
				
				if (cond_components[componentID] > cond) cond_components[componentID] = cond;
			}
			
			double conductance = 0;
			int cnt = 0;
			for (int i = 0; i < cond_components.length; i++) {
				if (cond_components[i] >= Constant.large) continue;
				conductance += cond_components[i];
				cnt++;
			}
			if (conductance == 0 || cnt == 0) conductance = 0;
			else conductance /= cnt;
			
			HConductances[trial] = conductance;
			
			latch.countDown();
		}
	}
		
	public void HConductance(int trials, int strategy) throws IOException, InterruptedException {
				
		if (!DiscoverLoaded) loadDiscover(trials, strategy);
		
		if (!GraphLoaded) {
			loadgraph();
		}
		
		HConductances = new double[trials];
		
		// for discover clusters
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			HConductance = new HConductance(trial);
		}
		latch.await();
		
		try {
			// for discover clusters
			FileWriter fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/HConductance_" + method + ".txt",true);
			for (int trial = 0; trial < trials; trial++) {
				fw_user.write(trial + "," + clusters_Discover[trial].length + "," + clusters_Ground.length + "," + String.format("%.4f", HConductances[trial]) + "\n");
			}
			fw_user.write("\n");
			fw_user.flush();
			fw_user.close();
			
			fw_user = new FileWriter(FilePath_Mon.filePathPre + "/measures/HConductance.txt",true);
			double avgDiscover = 0;
			double avgHCond = 0;
			int cnt = 0;
			// iteration level precision recall f1
			for (int trial = 0; trial < trials; trial++) {
				avgDiscover += clusters_Discover[trial].length;
				if (HConductances[trial] > 0) {
					avgHCond += HConductances[trial];
					cnt++;
				}
			}
			avgDiscover /= trials;
			avgHCond /= cnt;
			fw_user.write("HModularity_" + method + "\n");
			fw_user.write(String.format("%.4f",avgDiscover) + "," + clusters_Ground.length + "," + String.format("%.4f",avgHCond) + "\n");
			fw_user.write("\n");
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		HConductances = null;
		
		Hypergraph.garbbageCollector.gc();
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
		int trial_Louvain = 1;
		String moveStrategy = "move";
		boolean toHigherOrder = false;
		String ordering = "randomOrder";
		double ratio = 0.5;
		int strategy = 0;
		SingleGTruth calculator = new SingleGTruth(false, moveStrategy, strategy, toHigherOrder, ordering, ratio);
		calculator.DiscoverLoaded = false;
		calculator.GraphLoaded = false;
//		calculator.PRF(trial_Louvain, strategy);
		calculator.PRF_Prime(trial_Louvain, strategy);
		calculator.NMI(trial_Louvain, strategy);
//		calculator.Purity(trial_Louvain, strategy);
//		calculator.ARI(trial_Louvain, strategy);
//		calculator.ARI_Prime(trial_Louvain, strategy);
//		calculator.HModularity(trial_Louvain, strategy);
//		calculator.HConductance(trial_Louvain, strategy);
		
		System.out.println();
		
		int trial_BogLouvain = 1;
		strategy = 0;
		calculator = new SingleGTruth(true, "BogLouvain", strategy);
		calculator.DiscoverLoaded = false;
		calculator.GraphLoaded = false;
//		calculator.PRF(trial_BogLouvain, strategy);
		calculator.PRF_Prime(trial_BogLouvain, strategy);
		calculator.NMI(trial_BogLouvain, strategy);
//		calculator.Purity(trial_BogLouvain, strategy);
//		calculator.ARI(trial_BogLouvain, strategy);
//		calculator.ARI_Prime(trial_BogLouvain, strategy);
//		calculator.HConductance(trial_BogLouvain, strategy);
		
		System.out.println();
		
		int trial_BogCNMRan = 1;
		strategy = 0;
		calculator = new SingleGTruth(true, "BogCNMRan", strategy);
		calculator.DiscoverLoaded = false;
		calculator.GraphLoaded = false;
//		calculator.PRF(trial_BogCNMRan, strategy);
		calculator.PRF_Prime(trial_BogCNMRan, strategy);
		calculator.NMI(trial_BogCNMRan, strategy);
//		calculator.Purity(trial_BogCNMRan, strategy);
//		calculator.ARI(trial_BogCNMRan, strategy);
//		calculator.ARI_Prime(trial_BogCNMRan, strategy);
//		calculator.HConductance(trial_BogCNMRan, strategy);
		
		
//		int trial_BogCNMRan = 30;
//		int strategy = 2;
//		SingleGTruth calculator = new SingleGTruth(true, "BogCNMRan", strategy);
//		calculator.DiscoverLoaded = false;
//		calculator.PRF(trial_BogCNMRan, strategy);
//		calculator.NMI(trial_BogCNMRan, strategy);
//		calculator.Purity(trial_BogCNMRan, strategy);
//		calculator.ARI(trial_BogCNMRan, strategy);
//		
//		strategy = 0;
//		calculator = new SingleGTruth(true, "BogCNMRan", strategy);
//		calculator.DiscoverLoaded = false;
//		calculator.HypergraphLoaded = false;
//		calculator.HModularity(trial_BogCNMRan, strategy);
		
//		int trial_BogCNMOpt = 26;
//		calculator = new SingleGTruth(true, "BogCNMOpt");
//		calculator.DiscoverLoaded = false;
//		calculator.PRF(trial_BogCNMOpt);
//		calculator.NMI(trial_BogCNMOpt);
//		calculator.Purity(trial_BogCNMOpt);
//		calculator.ARI(trial_BogCNMOpt);
	}
}
