package hyperGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.Functions;

public class Hypergraph {

	public static String dataset;
	private static String fileInput_edge;
	private static String fileInput_inc;

	public static Runtime garbbageCollector = Runtime.getRuntime();
	// node - hyperedge matrix
	public static int[] INC_eID;
	public static double[] INC_weight;
	public static int[] INC_head;

	// hyperedge - node matrix
	public static int[] EINC_nID;
	public static double[] EINC_weight;
	public static int[] EINC_head;

	// others
	public static int[] array;
	public static int O_volG;
	
	//////////////////// temporal variables ////////////////////////

	public static List<List<Integer>> einc_nIDs;
	public static List<List<Double>> einc_weights;

	public static List<Integer>[] inc_eIDs; // temporal variable
	public static List<Double>[] inc_weights; // temporal variable

	// For loading a hypergraph
	public static void loadGraph() throws IOException, InterruptedException {
		
		// The paths to the data files
		if (!Constant.CONNECTED) {
			fileInput_edge = FilePath_Mon.filePathPre + "/edge.txt";
			fileInput_inc = FilePath_Mon.filePathPre + "/inc.txt";
			
		} else {
			fileInput_edge = FilePath_Mon.filePathPre + "/edge_connect.txt";
			fileInput_inc = FilePath_Mon.filePathPre + "/inc_connect.txt";
		}

		// extract the name of data set
		String[] strs = FilePath_Mon.filePathPre.split("/");
		dataset = strs[strs.length - 1];

		double start = System.currentTimeMillis();

		readFiles();	// read the data files
		storeInc();	// store the hypergraph in the form of incident matrix

		double end = System.currentTimeMillis();
		double runningTime = (end - start) / Constant.RUNNING_TIME_UNIT;
		String runningTime_update_print = String.format("%.2f", runningTime);

		System.out.println(fileInput_edge + " data set: " + dataset + " loaded. " + "vertices: " + getVertexSize() + " edge: " + getEdgeSize()
				+ " time: " + runningTime_update_print);
	}
	
	// For loading a synthetic hypergraph
	public static void loadsynGraph(double value, String type) throws IOException, InterruptedException {
		
		// The paths to the data files
		if (!Constant.CONNECTED) {
			fileInput_edge = FilePath_Mon.filePathPre + "/scalability/" + type + "_" + String.format("%.1f", value) + "/edge.txt";
			fileInput_inc = FilePath_Mon.filePathPre + "/scalability/" + type + "_" + String.format("%.1f", value) + "/inc.txt";
		} else {
			fileInput_edge = FilePath_Mon.filePathPre + "/scalability/" + type + "_" + String.format("%.1f", value) + "/edge_connect.txt";
			fileInput_inc = FilePath_Mon.filePathPre + "/scalability/" + type + "_" + String.format("%.1f", value) + "/inc_connect.txt";
		}

		System.out.println("fileInput_edge " + fileInput_edge);
		System.out.println("fileInput_inc " + fileInput_inc);
		
		// extract the name of data set
		String[] strs = FilePath_Mon.filePathPre.split("/");
		dataset = strs[strs.length - 1];

		double start = System.currentTimeMillis();

		readFiles();	// read the data files
		storeInc();	// store the hypergraph in the form of incident matrix

		double end = System.currentTimeMillis();
		double runningTime = (end - start) / Constant.RUNNING_TIME_UNIT;
		String runningTime_update_print = String.format("%.2f", runningTime);

		System.out.println(fileInput_edge + " data set: " + dataset + " loaded. " + "vertices: " + getVertexSize() + " edge: " + getEdgeSize()
				+ " time: " + runningTime_update_print);
	}

	// read the data files
	private static void readFiles() throws IOException {
		// vertex size
		Path path = Paths.get(fileInput_inc);
		int vertexSize = (int) Files.lines(path).count();

		einc_nIDs = new ArrayList<List<Integer>>();
		einc_weights = new ArrayList<List<Double>>();
		inc_eIDs = new ArrayList[vertexSize];
		inc_weights = new ArrayList[vertexSize];
		for (int i = 0; i < vertexSize; i++) {
			inc_eIDs[i] = new ArrayList<Integer>();
			inc_weights[i] = new ArrayList<Double>();
		}

		int edgeID = 0;
		String line;
		String[] strs;
		List<Integer> vs;
		List<Double> ws;
		int v;
		double weightPerNode;
		try (FileReader reader = new FileReader(fileInput_edge);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");

				if (strs.length < 2) continue;

				int cardinality = strs.length;
				vs = new ArrayList<Integer>(cardinality);
				ws = new ArrayList<Double>(cardinality);
				
				// 1
				weightPerNode = Constant.INITIAL_NODEONEDGE_WEIGHT;
				weightPerNode = weightPerNode / cardinality;

				for (int i = 0; i < cardinality; i++) {
					v = Integer.parseInt(strs[i]);
					vs.add(v);
					ws.add(weightPerNode);

					inc_eIDs[v].add(edgeID);
					inc_weights[v].add(weightPerNode);
				}

				// new an edge
				einc_nIDs.add(vs);
				einc_weights.add(ws);

				edgeID++;
			}
		}
	}

	// store the hypergraph in the form of incident matrix
	private static void storeInc() throws InterruptedException {
		int totalCardinality = 0;
		for (int i = 0; i < inc_eIDs.length; i++) {
			totalCardinality += inc_eIDs[i].size();
		}

		//////////////////////////////////////////////////////////////////
		// node - hyperedge matrix
		
		INC_head = new int[inc_eIDs.length];
		INC_eID = new int[totalCardinality];
		INC_weight = new double[totalCardinality];

		int head_ptr = 0;

		int vertexSize = Hypergraph.getVertexSize();
		int neighborSize;
		List<Integer> edgeIDs = new ArrayList<Integer>();
		List<Double> weights = new ArrayList<Double>();
		for (int nodeID = 0; nodeID < vertexSize; nodeID++) {

			neighborSize = inc_eIDs[nodeID].size();
			edgeIDs = inc_eIDs[nodeID];
			weights = inc_weights[nodeID];

			INC_head[nodeID] = head_ptr;

			for (int idx = 0; idx < neighborSize; idx++) {
				INC_eID[head_ptr] = edgeIDs.get(idx);
				INC_weight[head_ptr] = weights.get(idx);
				head_ptr++;
			}
		}

		//////////////////////////////////////////////////////////////////
		// hyperedge - node matrix
		
		EINC_head = new int[einc_nIDs.size()];

		int hyperedgeSize = Hypergraph.getEdgeSize();
		EINC_nID = new int[totalCardinality];
		EINC_weight = new double[totalCardinality];

		head_ptr = 0;

		List<Integer> vs;
		List<Double> ws;

		for (int edgeID = 0; edgeID < hyperedgeSize; edgeID++) {

			EINC_head[edgeID] = head_ptr;

			vs = einc_nIDs.get(edgeID);
			ws = einc_weights.get(edgeID);

			for (int idx = 0; idx < vs.size(); idx++) {
				EINC_nID[head_ptr] = vs.get(idx);
				EINC_weight[head_ptr] = ws.get(idx);
				head_ptr++;
			}
		}

		//////////////////////////////////////////////////////////////////
		// estimate the upper bound of volG, to alarm the overflow incurred by the large hypergraph size
		
		double sum = 0;
		int cardinality, comb;
		for (int edgeID = 0; edgeID < hyperedgeSize; edgeID++) {
			cardinality = einc_nIDs.get(edgeID).size();
			comb = cardinality * (cardinality-1);
			
			if ((sum + comb) >= (double) (Constant.maxArraySize - 1000)) {
				sum = Constant.maxArraySize;
				break;
			} else sum += comb;
		}
		
		sum = Math.max(sum, getVertexSize());
		sum = Math.max(sum, getEdgeSize());
		
		if (sum > (double) Constant.maxArraySize) {
			System.out.println("Caution, it is a very large hypergraph");
			O_volG = Constant.maxArraySize;
		} else O_volG = (int) sum;
		
		///////////////////////////////////////////////////////////////////////
		// clean the temporal variables
		
		inc_eIDs = null;
		inc_weights = null;
		einc_nIDs = null;
		einc_weights = null;

		garbbageCollector.gc();
	}

	public static int getVertexSize() {
		return INC_head.length;
	}

	public static int getEdgeSize() {
		return EINC_head.length;
	}

	public static int getCardinality(int edgeID) {
		if (edgeID == Hypergraph.getEdgeSize() - 1)
			return (EINC_nID.length - Hypergraph.EINC_head[edgeID]);
		else
			return (Hypergraph.EINC_head[edgeID + 1] - Hypergraph.EINC_head[edgeID]);
	}

	public static int getDegree(int nodeID) {
		if (nodeID == Hypergraph.getVertexSize() - 1)
			return (INC_eID.length - Hypergraph.INC_head[nodeID]);
		else
			return (Hypergraph.INC_head[nodeID + 1] - Hypergraph.INC_head[nodeID]);
	}

	public static double getWeightDegree(int nodeID) {
		double weight = 0;
		int first_idx = Hypergraph.INC_head[nodeID];
		int second_idx = getSecondIdx_INC(nodeID);

		for (int i = first_idx; i < second_idx; i++) {
			weight += INC_weight[i];
		}

		return weight;
	}

	public static int getSecondIdx_INC(int nodeID) {
		if (nodeID == Hypergraph.getVertexSize() - 1)
			return INC_eID.length;
		else
			return Hypergraph.INC_head[nodeID + 1];
	}

	public static int getSecondIdx_EINC(int edgeID) {
		if (edgeID == Hypergraph.getEdgeSize() - 1)
			return Hypergraph.EINC_nID.length;
		else
			return Hypergraph.EINC_head[edgeID + 1];
	}

	public static void main(String arg[]) throws IOException, InterruptedException {
//		Hypergraph.loadGraph();
//		
//		int nodeDegree;
//		int maxNodeDegree = -1;
//		for (int nodeID = 0; nodeID < Hypergraph.getVertexSize(); nodeID++) {
//			
//			int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[nodeID];
//			if (nodeID == Hypergraph.getVertexSize() - 1) secondIdx_INC = Hypergraph.INC_eID.length;
//			else secondIdx_INC = Hypergraph.INC_head[nodeID + 1];
//
//			nodeDegree = secondIdx_INC - firstIdx_INC;
//			if (maxNodeDegree < nodeDegree) maxNodeDegree = nodeDegree;
//		}
//		System.out.println("maxNodeDegree " + maxNodeDegree);
//		
//		double sumComb = 0;
//		int cardinality;
//		double comb;
//		for (int edgeID = 0; edgeID < Hypergraph.getEdgeSize(); edgeID++) {
//			
//			int secondIdx_EINC, firstIdx_EINC = Hypergraph.EINC_head[edgeID];
//			if (edgeID == Hypergraph.getEdgeSize() - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
//			else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];
//			
//			cardinality = secondIdx_EINC - firstIdx_EINC;
//			comb = Functions.combination(cardinality, 2);
//			sumComb += comb;
//		}
//		System.out.println("sumComb " + sumComb);
	}
}
