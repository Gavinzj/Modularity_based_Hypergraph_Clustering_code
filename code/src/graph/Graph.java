package graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.Functions;

public class Graph {
	
	public static boolean loaded;
	
	public static int vertexSize;
	public static int edgeSize;
	
	// node - node matrix
	public static int[] ADJ_head;
	public static int[] ADJ_nID;
	public static double[] ADJ_weight;
	
	public static int getVertexSize() {
		return vertexSize;
	}
	
	public static int getEdgeSize() {
		return edgeSize;
	}
	
	public static void reduceGraph() {
		vertexSize = Hypergraph.getVertexSize();
		int hyperEdgeSize = Hypergraph.getEdgeSize();

		HashMap<String, Integer> distinctEdges = new HashMap<String, Integer>();
		List<Double> weights = new ArrayList<Double>();
		List<Integer>[] inc_eIDs = new ArrayList[vertexSize];
		List<Integer>[] adj_nIDs = new ArrayList[vertexSize];
		for (int i = 0; i < vertexSize; i++) {
			inc_eIDs[i] = new ArrayList<Integer>();
			adj_nIDs[i] = new ArrayList<Integer>();
		}
		
		int edgeID = 0;
		int u, v, head, end, secondIdx_EINC, firstIdx_EINC;
		double edgeWeight;
		for (int eid = 0; eid < hyperEdgeSize; eid++) {
		    
			firstIdx_EINC = Hypergraph.EINC_head[eid];
			if (eid == Hypergraph.getEdgeSize() - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[eid + 1];
			
			// 1
			edgeWeight = Constant.INITIAL_EDGE_WEIGHT;
			edgeWeight /= (double) Functions.combination((secondIdx_EINC - firstIdx_EINC), 2);
			
			for (int i = firstIdx_EINC; i < secondIdx_EINC; i++) {
				u = Hypergraph.EINC_nID[i];
				
				for (int j = i + 1; j < secondIdx_EINC; j++) {
					v = Hypergraph.EINC_nID[j];
					
					if (u < v) {
						head = u;
						end = v;
					} else {
						head = v;
						end = u;
					}
					
					if (!distinctEdges.containsKey(u + "-" + v)) {
						distinctEdges.put(u + "-" + v, edgeID);
						weights.add(edgeWeight);
						
						inc_eIDs[head].add(edgeID);
						inc_eIDs[end].add(edgeID);
						
						adj_nIDs[head].add(end);
						adj_nIDs[end].add(head);
						
						edgeID++;
						
					} else {
						int idx = distinctEdges.get(u + "-" + v);
						double newWeight = weights.get(idx) + edgeWeight;
						weights.set(idx, newWeight);
					}
				}
			}
		}
		
		//////////////////////////////////////////////////////////
		
		edgeSize = distinctEdges.size();
		ADJ_head = new int[adj_nIDs.length];
		ADJ_nID = new int[2 * Graph.edgeSize];
		ADJ_weight = new double[2 * Graph.edgeSize];
		
		int head_ptr = 0;
		
		int neighborSize;
		List<Integer> vs, es;
		for (int nodeID = 0; nodeID < vertexSize; nodeID++) {

			neighborSize = adj_nIDs[nodeID].size();
			vs = adj_nIDs[nodeID];
			es = inc_eIDs[nodeID];

			ADJ_head[nodeID] = head_ptr;

			for (int idx = 0; idx < neighborSize; idx++) {
				ADJ_nID[head_ptr] = vs.get(idx);
				ADJ_weight[head_ptr] = weights.get(es.get(idx)) / 2;
				head_ptr++;
			}
		}
		
		loaded = true;
		System.out.println("clique reduced graph. " + "vertices: " + getVertexSize() + " edge: " + getEdgeSize());
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		Graph.reduceGraph();
    }
}
