package exp.sensitivity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.Functions;

public class ProcessDataSets {
	
	static double lambda_precision = 6;
	static double enlarger = Math.pow(10, lambda_precision);
	
	static double totalDegree;
	static double[] nodeDegrees;
	static int[] cardinalities;
	static List<Integer>[] EINC;
	static List<Integer>[] INC;
	
	public static void genRandomHypergraph_onLambda(double lambda) throws FileNotFoundException, IOException {
		
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		lambda = (double) Math.round(lambda * enlarger) / enlarger;
		
		// node degree distribution
		totalDegree = 0;
		nodeDegrees = new double[n];
		
		int node = 0;
		double degree;
    	for (node = 0; node < n - 1; node++) {
    		// for each incident hyperedge
			degree = Hypergraph.INC_head[node + 1] - Hypergraph.INC_head[node];
			nodeDegrees[node] = degree;
			totalDegree += degree;
    	}
		degree = Hypergraph.INC_weight.length - Hypergraph.INC_head[node];
		nodeDegrees[node] = degree;
		totalDegree += degree;
		
		System.out.println("generate random graph, n " + n + " m " + m + " total degree " + totalDegree);
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		
		readCardinalityDistribution(m, lambda);
		
		////////////////////////////////////////////////////////////////////////////////////////////
		
		generator(n, m);
		
		///////////////////////////////////////////////////////////////////////////////////////
		
		//save
		String fileOutput_edge_dir = FilePath_Mon.filePathPre + "/scalability/lambda_" + lambda;
		String fileOutput_inc_dir = FilePath_Mon.filePathPre + "/scalability/lambda_" + lambda;
		save(n, m, lambda, fileOutput_edge_dir, fileOutput_inc_dir);
	}
	
	public static void genRandomHypergraph_onNodeSize(int copies) throws IOException, InterruptedException {
		Random random = new Random();
		int n = Hypergraph.getVertexSize();
		int m = Hypergraph.getEdgeSize();
		int[] allNodes = new int[n];
		double[] allNodeDegrees = new double[n];
		
		double ratio = 1.0 / (double) copies;
		int nodeSizePerCopy = (int) Math.round(n * ratio);
		nodeSizePerCopy = Math.min(nodeSizePerCopy, n);
		System.out.println("copy number " + copies + " nodeSizePerCopy " + nodeSizePerCopy + " total node size " + n);
		
		// get lambda
		String dataset = Hypergraph.dataset;
		double lambda = getLambda(dataset);
		if (lambda == -999) {
			System.out.println("wrong lambda");
			return;
		}
		System.out.println("dataset " + dataset + " lambda " + lambda);
		
		int node, edge, newNode, newEdge;
		
		// randomize nodes
		for (int i = 0; i < n; i++) {
			allNodes[i] = i;
		}
		int j, temp;
		for (int i = 0; i < n; i++) {
			j = random.nextInt(n);
			temp = allNodes[i];
			allNodes[i] = allNodes[j];
			allNodes[j] = temp;
		}
		
		// get node degrees
		int first_idx, second_idx;
		for (int i = 0; i < n; i ++) {
			node = allNodes[i];
			first_idx = Hypergraph.INC_head[node];
			second_idx = Hypergraph.getSecondIdx_INC(node);
			
			allNodeDegrees[i] = (second_idx - first_idx);
		}
		
		EINC = new ArrayList[m];
		INC = new ArrayList[n];
		int newNodeID;
		int[] newNodeIDs = new int[n];
		int newEdgeID;
		int[] newEdgeIDs = new int[m];
		
		HashSet<Integer> selectedEdges = null;
		int[] selectedEdgesFrequency = null;
		
		double thisRatio;
		int thisNodeSize = 0, thisEdgeSize = 0;
		for (int copy = 1; copy <= copies; copy++) {
			thisRatio = ratio * copy;
			thisRatio = (double) Math.round(thisRatio * enlarger) / enlarger;
			thisNodeSize = (int) Math.min(n, (n * thisRatio));
			
			selectedEdges = new HashSet<Integer>();
			selectedEdgesFrequency = new int[m];
			
			// for each selected node
			for (int i = 0; i < thisNodeSize; i++) {
				node = allNodes[i];
				
				first_idx = Hypergraph.INC_head[node];
				second_idx = Hypergraph.getSecondIdx_INC(node);
				for (int k = first_idx; k < second_idx; k++) {
					edge = Hypergraph.INC_eID[k];
					selectedEdgesFrequency[edge]++;
					
					if (selectedEdgesFrequency[edge] >= 2) selectedEdges.add(Hypergraph.INC_eID[k]);
				}
			}
			thisEdgeSize = (int) Math.min(m, selectedEdges.size());
			System.out.println("this ratio " + thisRatio + " thisNodeSize " + thisNodeSize + " thisEdgeSize " + thisEdgeSize);
			
			for (int i = 0; i < n; i++) INC[i] = new ArrayList<Integer>();
			for (int i = 0; i < m; i++) EINC[i] = new ArrayList<Integer>();
			newNodeID = 0;
			for (int i = 0; i < n; i++) newNodeIDs[i] = -1;
			newEdgeID = 0;
			for (int i = 0; i < m; i++) newEdgeIDs[i] = -1;
			
			// for each selected node
			boolean nodeExist, edgeExist;
			for (int i = 0; i < thisNodeSize; i++) {
				node = allNodes[i];
				if (newNodeIDs[node] == -1) newNodeIDs[node] = newNodeID++;
				newNode = newNodeIDs[node];
				nodeExist = false;
				
				first_idx = Hypergraph.INC_head[node];
				second_idx = Hypergraph.getSecondIdx_INC(node);
				for (int k = first_idx; k < second_idx; k++) {
					edge = Hypergraph.INC_eID[k];
					if (newEdgeIDs[edge] == -1) newEdgeIDs[edge] = newEdgeID++;
					newEdge = newEdgeIDs[edge];
					edgeExist = false;
					
					if (selectedEdgesFrequency[edge] >= 2) {
						nodeExist = true;
						edgeExist = true;
						INC[newNode].add(newEdge);
						EINC[newEdge].add(newNode);
					}
					
					if (!edgeExist) {
						newEdgeIDs[edge] = -1;
						newEdgeID--;
					}
				}
				
				if (!nodeExist) {
					newNodeIDs[node] = -1;
					newNodeID--;
				}
			}
			
			System.out.println("newEdgeID " + newEdgeID + " newNodeID " + newNodeID);
			
			//////////////////////////////////////////////////////////////////////////
			//save
			
			String fileOutput_edge_dir = FilePath_Mon.filePathPre + "/scalability/divide_" + thisRatio;
			String fileOutput_inc_dir = FilePath_Mon.filePathPre + "/scalability/divide_" + thisRatio;
			
			File theDir = new File(fileOutput_edge_dir);
			if (!theDir.exists()) theDir.mkdirs();
			theDir = new File(fileOutput_inc_dir);
			if (!theDir.exists()) theDir.mkdirs();
			
			String fileOutput_edge;
			String fileOutput_inc;
			if (Constant.CONNECTED) {
				fileOutput_edge = fileOutput_edge_dir + "/edge_connect.txt";
				fileOutput_inc = fileOutput_inc_dir + "/inc_connect.txt";
			} else {
				fileOutput_edge = fileOutput_edge_dir + "/edge.txt";
				fileOutput_inc = fileOutput_inc_dir + "/inc.txt";
			}
			
			List<Integer> vs_sorted;
			String line_output;
			
			try {
				FileWriter fwCount = new FileWriter(fileOutput_edge);
				
				int cnt = 0;
				for (edge = 0; edge < m; edge++) {
					if (EINC[edge].size() > 0) {
						cnt++;
						vs_sorted = EINC[edge];
						Collections.sort(vs_sorted);
						
						line_output = vs_sorted.get(0) + "";
				        for (int i = 1; i < vs_sorted.size(); i++) {
				        	line_output = line_output + "\t" + vs_sorted.get(i);
				        }
				        fwCount.write(line_output + "\n");
					}
				}
				fwCount.close();
				
				System.out.println("final edge size " + cnt);
				
			} catch (IOException e) {
				e.printStackTrace();
				return;
			}
			
			
			////////////////////////////////////////////////////////////////////////////////////////////
			
			try {
				FileWriter fwCount = new FileWriter(fileOutput_inc);
				
				int cnt = 0;
				for (node = 0; node < n; node++) {
					if (INC[node].size() > 0) {
						cnt++;
						vs_sorted = INC[node];
						Collections.sort(vs_sorted);
						
						line_output = vs_sorted.stream().map(Object::toString).collect(Collectors.joining("\t"));
						fwCount.write(node + "\t" + line_output + "\n");
					}
				}
				fwCount.close();
				
				System.out.println("final node size " + cnt);
				
			} catch (IOException e) {
				e.printStackTrace();
				return;
			}
		}
	}
	
	public static void synCardinalityDistribution (int edgeSize, double lambda) {
		
		lambda = (double) Math.round(lambda * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		System.out.println("synthetic cardinality distribution, edge size " + edgeSize + " lambda " + lambda);
		
		int trials = 10;
		for (int trial = 0; trial < trials; trial++) {
			
			String saveFilePath;
			if (Constant.CONNECTED) {
				saveFilePath = "data/" + Hypergraph.dataset + "/scalability/synCardinality_connected_lambda_" + lambda + "_trial_" + trial + ".txt";
			} else {
				saveFilePath = "data/" + Hypergraph.dataset + "/scalability/synCardinality_lambda_" + lambda + "_trial_" + trial + ".txt";
			}
			
			Process mProcess = null;
			try {
				Process process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/generateRandom.py " 
						+ Hypergraph.dataset + " " + lambda + " " + edgeSize + " " + saveFilePath);
				mProcess = process;
			} catch (Exception e) {
				System.out.println("Exception Raised" + e.toString());
			}
			InputStream stdout = mProcess.getInputStream();
			BufferedReader reader = new BufferedReader(new InputStreamReader(stdout, StandardCharsets.UTF_8));
			String line;
			String prevLine = "";
			try {
				System.out.println("stdout: "+ reader.readLine());
				while((line = reader.readLine()) != null){
					prevLine = line;
				}
				System.out.println("stdout: "+ prevLine);
			} catch (IOException e) {
				System.out.println("Exception in reading output" + e.toString());
			}
		}
	}
	
	public static void synDegreeDistribution (int nodeSize, double alpha) {
		
		alpha = (double) Math.round(alpha * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		System.out.println("synthetic degree distribution, edge size " + nodeSize + " alpha " + alpha);
		
		int trials = 20;
		for (int trial = 0; trial < trials; trial++) {
			
			String saveFilePath;
			if (Constant.CONNECTED) {
				saveFilePath = "data/" + Hypergraph.dataset + "/scalability/synDegree_connected_alpha_" + alpha + "_trial_" + trial + ".txt";
			} else {
				saveFilePath = "data/" + Hypergraph.dataset + "/scalability/synDegree_alpha_" + alpha + "_trial_" + trial + ".txt";
			}
			
			Process mProcess = null;
			try {
				Process process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/generateRandomPowLaw.py " 
						+ Hypergraph.dataset + " " + alpha + " " + nodeSize + " " + saveFilePath);
				mProcess = process;
			} catch (Exception e) {
				System.out.println("Exception Raised" + e.toString());
			}
			InputStream stdout = mProcess.getInputStream();
			BufferedReader reader = new BufferedReader(new InputStreamReader(stdout, StandardCharsets.UTF_8));
			String line;
			String prevLine = "";
			try {
				System.out.println("stdout: "+ reader.readLine());
				while((line = reader.readLine()) != null){
					prevLine = line;
				}
				System.out.println("stdout: "+ prevLine);
			} catch (IOException e) {
				System.out.println("Exception in reading output" + e.toString());
			}
		}
	}
	
	private static void readCardinalityDistribution(int edgeSize, double lambda) throws FileNotFoundException, IOException {
		// hyperedge cardinality distribution
		cardinalities = new int[edgeSize];
		
		System.out.println("lambda " + lambda);
		
		String file_cardinality;
		if (Constant.CONNECTED) {
			file_cardinality = FilePath_Mon.filePathPre + "/scalability/synCardinality_connected_lambda_" + lambda + ".txt";
		} else {
			file_cardinality = FilePath_Mon.filePathPre + "/scalability/synCardinality_lambda_" + lambda + ".txt";
		}
		
		File f = new File(file_cardinality);
		if(!f.exists()) {
			System.out.println("synCardinalityDistribution");
			synCardinalityDistribution(Hypergraph.getEdgeSize(), lambda);
		}
		
		int edgeCnt = 0;
		try (FileReader reader = new FileReader(file_cardinality);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line;
			while ((line = bufferedReader.readLine()) != null) {
				cardinalities[edgeCnt++] = Integer.parseInt(line);
			}
		}
		
		if (edgeCnt < edgeSize) {
			System.out.println("wrong");
			return;
		}
	}
	
	private static void generator(int n, int m) {
		
		for (int i = 0; i < n; i++) {
			nodeDegrees[i] /= totalDegree;
		}
		
		int factor = (int) Math.pow(10, 9);
		
		int degree_int;
		int rollBinSize = 0;
		for (int i = 0; i < n; i++) {
			nodeDegrees[i] = Math.round(nodeDegrees[i] * factor);
			degree_int = (int) nodeDegrees[i];
			rollBinSize += degree_int;
		}
		
		int[] rollBin = new int[rollBinSize];
		int idx = 0;
		for (int i = 0; i < n; i++) {
			degree_int = (int) nodeDegrees[i];
			for (int j = 0; j < degree_int; j++) {
				rollBin[idx++] = i;
			}
		}
		
		System.out.println("rollBin size " + rollBinSize);
		
		EINC = new ArrayList[m];
		INC = new ArrayList[n];
		
		for (int i = 0; i < n; i++) {
			INC[i] = new ArrayList<Integer>();
		}
		
		Random random = new Random();
		
		// generate random hypergraph
		// for each hyperedge
		int cardinality, node;
		boolean[] visited = new boolean[n];
		int distinctNodeSize;
		List<Integer> distinctNodes;
		
		for (int edge = 0; edge < m; edge++) {
			cardinality = cardinalities[edge];
			EINC[edge] = new ArrayList<Integer>(cardinality);
			
			distinctNodeSize = 0;
			distinctNodes = new ArrayList<Integer>(cardinality);
			
			while (distinctNodeSize < cardinality) {
				// randomly pick a node
				node = rollBin[random.nextInt(rollBinSize)];
				
				if (!visited[node]) {
					visited[node] = true;
					distinctNodeSize++;
					distinctNodes.add(node);
				}
			}
			
			for (int distinctNode : distinctNodes) {
				EINC[edge].add(distinctNode);
				INC[distinctNode].add(edge);
				
				visited[distinctNode] = false;
			}
		}
	}
	
	private static void generator_largeGraph(int n, int m) {
		
		for (int i = 0; i < n; i++) {
			nodeDegrees[i] /= totalDegree;
		}
		System.out.println("totalDegree " + totalDegree);
		
		int maxValue = Integer.MAX_VALUE - 5;
		int split = 0;
		double totalRollBinSize = 0;
		int[] largeRollBinSize = new int[32];
		int[] largeSplit = new int[n];
		
		int factor = (int) Math.pow(10, 9);
		System.out.println("factor " + factor);
		int degree_int;
		for (int node = 0; node < n; node++) {
			nodeDegrees[node] = Math.round(nodeDegrees[node] * factor);
			degree_int = (int) nodeDegrees[node];
			
			if ((largeRollBinSize[split] + degree_int) < maxValue) {
				largeRollBinSize[split] += degree_int;
				totalRollBinSize += degree_int;
				largeSplit[node] = split;
			} else {
				split++;
				largeRollBinSize[split] += degree_int;
				totalRollBinSize += degree_int;
				largeSplit[node] = split;
			}
		}
		
		int[][] largeRollBin = new int[split+1][Integer.MAX_VALUE - 5];
		int idx = 0;
		for (int node = 0; node < n; node++) {
			degree_int = (int) nodeDegrees[node];
			split = largeSplit[node];
			
			for (int i = 0; i < degree_int; i++) {
				largeRollBin[split][idx++] = node;
			}
		}
		System.out.println("totalRollBinSize " + totalRollBinSize);
		
		int splitRollBinSize = 1000000;
		int[] splitRollBin = new int[splitRollBinSize];
		idx = 0;
		for (split = 0; split < largeRollBin.length; split++) {
			double size = largeRollBinSize[split];
			int adjustedSize = (int) (splitRollBinSize * size / totalRollBinSize);
			System.out.println("split " + split + " rollBin size " + size + " adjustedSize " + adjustedSize);
			for (int i = 0; i < adjustedSize; i++) {
				splitRollBin[idx++] = split;
			}
		}
		System.out.println("idx " + idx);
		
		EINC = new ArrayList[m];
		INC = new ArrayList[n];
		
		for (int i = 0; i < n; i++) {
			INC[i] = new ArrayList<Integer>();
		}
		
		Random random = new Random();
		
		// generate random hypergraph
		// for each hyperedge
		int cardinality, node, rollBinSize;
		boolean[] visited = new boolean[n];
		int distinctNodeSize;
		List<Integer> distinctNodes;
		
		for (int edge = 0; edge < m; edge++) {
			cardinality = cardinalities[edge];
			EINC[edge] = new ArrayList<Integer>(cardinality);
			
			distinctNodeSize = 0;
			distinctNodes = new ArrayList<Integer>(cardinality);
			
			while (distinctNodeSize < cardinality) {
				// randomly pick a split
				split = splitRollBin[random.nextInt(splitRollBinSize)];
				rollBinSize = largeRollBinSize[split];
				// randomly pick a node
				node = largeRollBin[split][random.nextInt(rollBinSize)];
				
				if (!visited[node]) {
					visited[node] = true;
					distinctNodeSize++;
					distinctNodes.add(node);
				}
			}
			
			for (int distinctNode : distinctNodes) {
				EINC[edge].add(distinctNode);
				INC[distinctNode].add(edge);
				
				visited[distinctNode] = false;
			}
		}
	}
	
	public static double getLambda(String dataset) {
		double lambda = -999;
		
		// get cardinality distribution
		try {
			FileWriter fwCount = new FileWriter(
					FilePath_Mon.filePathPre + "/stat/cardinality/cardinality_trial_0.txt");
			
			int m = Hypergraph.getEdgeSize();
			int endM = m - 1;
			int EINC_nID_length = Hypergraph.EINC_nID.length;
			
			int cardinality;
			for (int edgeID = 0; edgeID < m; edgeID++) {

				if (edgeID == endM) cardinality = (EINC_nID_length - Hypergraph.EINC_head[edgeID]);
				else cardinality = (Hypergraph.EINC_head[edgeID + 1] - Hypergraph.EINC_head[edgeID]);

				if (cardinality >= 2) fwCount.write(cardinality + "\n");
			}
			fwCount.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -999;
		}

		double extractClusterTime = -999;
		// ref:https://www.edureka.co/community/358/how-to-execute-a-python-file-with-few-arguments-in-java
		Process process = null;
		try {
			process = Runtime.getRuntime().exec("python " + System.getProperty("user.dir") + "/python/getLambda.py "
					+ "/data/" + dataset + "/stat/cardinality/cardinality_trial_0.txt" + " " + 2);
		} catch (Exception e) {
			System.out.println("Exception Raised" + e.toString());
		}
		InputStream stdout = process.getInputStream();
		BufferedReader reader = new BufferedReader(new InputStreamReader(stdout, StandardCharsets.UTF_8));
		String line;
		try {

			// read the first line
			line = reader.readLine();
			if (line != null) {
				extractClusterTime = Double.parseDouble(line.split(",")[1]);
			} else return -999;

			// read the second line
			line = reader.readLine();
			if (line != null && !line.equals("nan")) {
				lambda = Double.parseDouble(line);
			} else return -999;
			
		} catch (IOException e) {
			System.out.println("Exception in reading output" + e.toString());
		}

		return lambda;
	}
	
	private static void save(int n, int m, double lambda, String fileOutput_edge_dir, String fileOutput_inc_dir) {
		
		File theDir = new File(fileOutput_edge_dir);
		if (!theDir.exists()) theDir.mkdirs();
		theDir = new File(fileOutput_inc_dir);
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput_edge;
		String fileOutput_inc;
		if (Constant.CONNECTED) {
			fileOutput_edge = fileOutput_edge_dir + "/edge_connect.txt";
			fileOutput_inc = fileOutput_inc_dir + "/inc_connect.txt";
		} else {
			fileOutput_edge = fileOutput_edge_dir + "/edge.txt";
			fileOutput_inc = fileOutput_inc_dir + "/inc.txt";
		}
		
		List<Integer> vs_sorted;
		String line_output;
		
		double vol = 0;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_edge);
			
			for (int edgeID = 0; edgeID < m; edgeID++) {
				vs_sorted = EINC[edgeID];
				Collections.sort(vs_sorted);
				
				vol += vs_sorted.size();
				
				line_output = vs_sorted.get(0) + "";
		        for (int i = 1; i < vs_sorted.size(); i++) {
		        	line_output = line_output + "\t" + vs_sorted.get(i);
		        }
		        fwCount.write(line_output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		System.out.println("vol " + vol);
		
		////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_inc);
			
			for (int nodeID = 0; nodeID < n; nodeID++) {
				vs_sorted = INC[nodeID];
				Collections.sort(vs_sorted);
				
				line_output = vs_sorted.stream().map(Object::toString).collect(Collectors.joining("\t"));
				fwCount.write(nodeID + "\t" + line_output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws Exception {
		Hypergraph.loadGraph();
		
//		for (double lambda = 0.45; lambda <= 0.65; lambda += 0.05) {
//			synCardinalityDistribution(6, lambda);
//		}
		
		double alpha = 2.030879436;
		synDegreeDistribution(11, alpha);
		
//		double lambda = 0.1;
//		genRandomHypergraph_onLambda(lambda);
		
//		int divide = 5;
//		divideNodeSet(divide);
	}
}
