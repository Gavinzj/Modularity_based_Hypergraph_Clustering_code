package preprocess;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;

public class GetCoreGraph {
	static int K = 2;
	
	static int n;
	static int m;
	
	static List<Integer>[] INC_eID;
	static HashSet<Integer> deletedNodes;
	
	static List<Integer>[] EINC_nID;
	static HashSet<Integer> deletedEdges;
	
	static int[] clusters;
	
	public static void loadGraph() throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
		n = Hypergraph.getVertexSize();
		m = Hypergraph.getEdgeSize();
		
		INC_eID = new ArrayList[n];
		for (int node = 0; node < n; node++) {
			INC_eID[node] = new ArrayList<Integer>();
			
			// for each incident edge
			int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[node];
			if (node == n - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[node + 1];
			
//			System.out.print("node " + nodeID + ": ");
			for (int i = firstIdx_INC; i < secondIdx_INC; i++) {
				int edgeID = Hypergraph.INC_eID[i];
				
//				System.out.print(edgeID + ",");
				
				INC_eID[node].add(edgeID);
			}
//			System.out.println();
//			System.out.println(INC_eID[nodeID]);
		}
		
		EINC_nID = new ArrayList[m];
		for (int edge = 0; edge < m; edge++) {
			EINC_nID[edge] = new ArrayList<Integer>();
			
			// for each node in the edge
			int secondIdx_EINC, firstIdx_EINC = Hypergraph.EINC_head[edge];
			if (edge == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[edge + 1];
			
//			System.out.print("edge " + edgeID + ": ");
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				int nodeID = Hypergraph.EINC_nID[i];
				
//				System.out.print(nodeID + ",");
				
				EINC_nID[edge].add(nodeID);
			}
//			System.out.println();
//			System.out.println(EINC_nID[edgeID]);
		}
	}
	public static void KCoreGraph() throws IOException, InterruptedException {
		loadGraph();
		
		int minDegree = 999999999;
		// for each node
		for (int node = 0; node < n; node++) {
			int degree = INC_eID[node].size();
			if (minDegree > degree) minDegree = degree;
		}
		
		deletedNodes = new HashSet<Integer>();
		deletedEdges = new HashSet<Integer>();
		
//		System.out.println("min degree " + minDegree);
		
		int itr = 0;
		while (minDegree < K) {
			System.out.println("itr " + itr++);
			
			// get removable nodes
			List<Integer> removeNodes = new ArrayList<Integer>();
			
			// for each node
			for (int node = 0; node < n; node++) {
				if (deletedNodes.contains(node)) {
//					System.out.println("node " + node + " deleted");
					continue;
				}
				int degree = INC_eID[node].size();
				if (degree < K) {
//					System.out.println("node " + node + ": " + INC_eID[node] + " remove");
					removeNodes.add(node);
					deletedNodes.add(node);
				}
			}
			
//			System.out.println("# of node to remove " + removeNodes.size() + ": " + removeNodes);
			
			// for each node, remove node from edge
			for (int removeNode : removeNodes) {
				List<Integer> incidentEdges = INC_eID[removeNode];
				
//				System.out.println("remove node " + removeNode + " edges " + incidentEdges);
				
				for (int edge : incidentEdges) {
					
//					System.out.println("for edge " + edge + ": " + EINC_nID[edge]);
					
					int idx = EINC_nID[edge].indexOf(removeNode);
					
//					System.out.println("index " + idx);
					
					EINC_nID[edge].remove(idx);
					
//					System.out.println("after: " + EINC_nID[edge]);
				}
			}
			
			// get removable edges
			List<Integer> removeEdges = new ArrayList<Integer>();
			
			// for each edge
			for (int edge = 0; edge < m; edge++) {
				if (deletedEdges.contains(edge)) {
//					System.out.println("edge " + edge + " deleted");
					continue;
				}
				int cardinality = EINC_nID[edge].size();
				if (cardinality <= 1) {
//					System.out.println("edge " + edge + ": " + EINC_nID[edge] + " remove");
					removeEdges.add(edge);
					deletedEdges.add(edge);
				}
			}
			
//			System.out.println("# of edge to remove " + removeEdges.size() + ": " + removeEdges);
			
			// for each edge, remove edge from its nodes
			for (int removeEdge : removeEdges) {
				List<Integer> nodes = EINC_nID[removeEdge];
				
//				System.out.println("remove edge " + removeEdge + " nodes " + nodes);
				
				for (int node : nodes) {
					
//					System.out.println("for node " + node + ": " + INC_eID[node]);
					
					int idx = INC_eID[node].indexOf(removeEdge);
					
//					System.out.println("index " + idx);
					
					INC_eID[node].remove(idx);
					
//					System.out.println("after: " + INC_eID[node]);
				}
			}
			
			// calculate min degree
			minDegree = 999999999;
			// for each node
			for (int node = 0; node < n; node++) {
				if (deletedNodes.contains(node)) continue;
				
				int degree = INC_eID[node].size();
				if (minDegree > degree) minDegree = degree;
			}
			
//			System.out.println("min degree " + minDegree);
		}
	}
	
	public static void loadGroundTruth() throws FileNotFoundException, IOException {
		
		clusters = new int[Hypergraph.getVertexSize()];
		
		String fileInput_clustering = "";
		if (!Constant.CONNECTED) fileInput_clustering = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		else fileInput_clustering = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		int nodeCnt = 0;
		int clusterCnt = 0;
		try (FileReader reader = new FileReader(fileInput_clustering);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				// process each line
				strs = line.split("\t");
				
//				System.out.println("line " + line);
				
				for (int i = 0; i < strs.length; i++) {
					int node = Integer.parseInt(strs[i]);
					
//					System.out.print(node + ",");
							
					clusters[node] = clusterCnt;
					nodeCnt++;
				}
//				System.out.println();
				
				clusterCnt++;
			}
		}
		
		System.out.println("nodeCnt " + nodeCnt + " clusterCnt " + clusterCnt);
	}
	
	public static void save() {
		String[] strs = FilePath_Mon.filePathPre.split("/");
		String dataset = strs[strs.length - 1];
		System.out.println("dataset " + dataset);
		
		String fileOutput_edges = FilePath_Mon.filePathPre + "/" + dataset + "_" + K + "core.txt";
		String fileOutput_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth" + "_" + K + "core.txt";
		
		int nodeCnt = 0;
		int edgeCnt = 0;
		HashMap<Integer, Integer> getNodeIDs = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> getNodes = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> edges = new HashMap<Integer, Integer>();
		
		// for each node
		for (int node = 0; node < n; node++) {
			if (deletedNodes.contains(node)) continue;
			getNodes.put(nodeCnt, node);
			getNodeIDs.put(node, nodeCnt++);
			
//			System.out.println("node " + node + ": " + INC_eID[node]);
		}
		
		// for each edge
		for (int edge = 0; edge < m; edge++) {
			if (deletedEdges.contains(edge)) continue;
			edges.put(edge, edgeCnt++);
//			System.out.println("edge " + edge + ": " + EINC_nID[edge]);
		}
		
		System.out.println(K + "-core, # of nodes " + getNodeIDs.size() + " " + nodeCnt + " # of edges " + edges.size() + " " + edgeCnt);
		
		/////////////////////////////////////////////////////////
		
		// save edge
		try {
			FileWriter fwCount = new FileWriter(fileOutput_edges);

			Iterator<Integer> itr_edge = edges.keySet().iterator();
			while (itr_edge.hasNext()) {
				int edge = itr_edge.next();
				List<Integer> nodesInEdge = EINC_nID[edge];
//				System.out.println("edge " + edge + " " + edges.get(edge) + ": " + nodesInEdge);
				
				int node = nodesInEdge.get(0);
				int newID = getNodeIDs.get(node);
//				System.out.println("node " + node + " newID " + newID);
				String output = newID + "";
				for (int i = 1; i < nodesInEdge.size(); i++) {
					node = nodesInEdge.get(i);
					newID = getNodeIDs.get(node);
//					System.out.println("node " + node + " newID " + newID);
					output = output + "\t" + newID;
				}
				
//				System.out.println(output + "\n");
				fwCount.write(output + "\n");
			}
			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		/////////////////////////////////////////////////////////
		
		// save ground truth
		try {
			FileWriter fwCount = new FileWriter(fileOutput_groundTruth);

			for (int newID = 0; newID < nodeCnt; newID++) {
				int node = getNodes.get(newID);
				int label = clusters[node];
//				System.out.println("newID " + newID + " node " + node + " label " + label);
				fwCount.write(newID + "\t" + label + "\n");
			}
			
			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		KCoreGraph();
		loadGroundTruth();
		save();
	}
}
