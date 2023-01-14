package statistic;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.Functions;

public class Statistic_Hypergraph {
	
	public static void basicInfo() {
		// data set name
		System.out.println("data set: " + Hypergraph.dataset);
		
		// number of nodes
		int n = Hypergraph.getVertexSize();
		System.out.println("number of nodes: " + n);
		
		// number of hyperedges
		int m = Hypergraph.getEdgeSize();
		System.out.println("number of hyperedges: " + m);
		
		// average/maximum/minimum cardinality
		double averageCardinality = 0;
		double maxCardinality = -1;
		double minCardinality = Constant.large;
		int cardinality;
		double sumComb = 0;
		double comb;
		for (int edge = 0; edge < m; edge++) {
			
			int secondIdx_EINC, firstIdx_EINC = Hypergraph.EINC_head[edge];
			if (edge == m - 1) secondIdx_EINC = Hypergraph.EINC_nID.length;
			else secondIdx_EINC = Hypergraph.EINC_head[edge + 1];
			
			cardinality = secondIdx_EINC - firstIdx_EINC;
			comb = Functions.combination(cardinality, 2);
			sumComb += comb;
			
			averageCardinality += cardinality;
			if (maxCardinality < cardinality) maxCardinality = cardinality;
			if (minCardinality > cardinality) minCardinality = cardinality;
		}
		System.out.println("total cardinality: " + String.format("%.5f",averageCardinality));
		System.out.println("average cardinality: " + String.format("%.5f",(averageCardinality / m)));
		System.out.println("minimum cardinality: " + String.format("%.5f",minCardinality));
		System.out.println("maximum cardinality: " + String.format("%.5f",maxCardinality));
		System.out.println("sumComb " + sumComb);
		
		// average/maximum/minimum node degree
		double averageDegree = 0;
		int maxNodeDegree = -1;
		int minDegree = (int) Constant.large;
		int nodeDegree;
		for (int node = 0; node < n; node++) {
			
			int secondIdx_INC, firstIdx_INC = Hypergraph.INC_head[node];
			if (node == Hypergraph.getVertexSize() - 1) secondIdx_INC = Hypergraph.INC_eID.length;
			else secondIdx_INC = Hypergraph.INC_head[node + 1];

			nodeDegree = secondIdx_INC - firstIdx_INC;
			averageDegree += nodeDegree;
			if (maxNodeDegree < nodeDegree) maxNodeDegree = nodeDegree;
			if (minDegree > nodeDegree) minDegree = nodeDegree;
		}
		System.out.println("total degree: " + String.format("%.5f",averageDegree));
		System.out.println("average node degree: " + String.format("%.5f",(averageDegree / n)));
		System.out.println("minimum node degree " + minDegree);
		System.out.println("maximum node degree " + maxNodeDegree);
		
		// for clique reduction graph
		forCliqueReduceGraph_largeGraph(n, m);
	}
	
	public static void forCliqueReduceGraph(int n, int m) {
		
		boolean[] visited = new boolean[n];
		int endN = n-1;
		int endM = m-1;
		int INC_eID_length = Hypergraph.INC_eID.length;
		int EINC_nID_length = Hypergraph.EINC_nID.length;
		
		int lastHead = 0;
		int head = 0;
		int arrayIndex = 0;

		int maxDegree = -1;
		int secondIdx_INC, secondIdx_EINC, edgeID, neighbor, degree;
		for (int thisNode = 0; thisNode < n; thisNode++) {

			lastHead = head;

			// for each incident edge
			if (thisNode == endN) secondIdx_INC = INC_eID_length;
			else secondIdx_INC = Hypergraph.INC_head[thisNode + 1];
			
			for (int i = Hypergraph.INC_head[thisNode]; i < secondIdx_INC; i++) {
				edgeID = Hypergraph.INC_eID[i];

				// for each node in the edge
				if (edgeID == endM) secondIdx_EINC = EINC_nID_length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];

				for (int j = Hypergraph.EINC_head[edgeID]; j < secondIdx_EINC; j++) {
					neighbor = Hypergraph.EINC_nID[j];
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
		
		System.out.println("clique-reduced graph, total vol " + arrayIndex);
		System.out.println("clique-reduced graph, maximum degree " + maxDegree);
	}
		
	public static void forCliqueReduceGraph_largeGraph(int n, int m) {
		
		boolean[] visited = new boolean[n];
		int endN = n-1;
		int endM = m-1;
		int INC_eID_length = Hypergraph.INC_eID.length;
		int EINC_nID_length = Hypergraph.EINC_nID.length;
		
		int lastHead = 0;
		int maxArraySize = Hypergraph.maxArraySize;
		int arrayIndex = 0;
		
		boolean full = false;
		int maxDegree = -1;
		long totVol = 0;
		int secondIdx_INC, secondIdx_EINC, edgeID, neighbor, degree;
		for (int thisNode = 0; thisNode < n; thisNode++) {
			
			lastHead = arrayIndex;
			
			// for each incident edge
			if (thisNode == endN) secondIdx_INC = INC_eID_length;
			else secondIdx_INC = Hypergraph.INC_head[thisNode + 1];
			
			for (int i = Hypergraph.INC_head[thisNode]; i < secondIdx_INC; i++) {
				edgeID = Hypergraph.INC_eID[i];

				// for each node in the edge
				if (edgeID == endM) secondIdx_EINC = EINC_nID_length;
				else secondIdx_EINC = Hypergraph.EINC_head[edgeID + 1];

				for (int j = Hypergraph.EINC_head[edgeID]; j < secondIdx_EINC; j++) {
					neighbor = Hypergraph.EINC_nID[j];
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
				
				for (int i = lastHead; i < arrayIndex; i++) {
					visited[Hypergraph.array[i]] = false;
				}
				totVol += lastHead;
				
				////////////////////////////////////////////////////////
				
				thisNode--;
				full = false;
				
				////////////////////////////////////////////////////////
				
				lastHead = 0;
				arrayIndex = 0;
				
			} else {
				
				degree = arrayIndex - lastHead;
				if (maxDegree < degree) maxDegree = degree;
				
				for (int i = lastHead; i < arrayIndex; i++) {
					visited[Hypergraph.array[i]] = false;
				}
			}
		}
		for (int i = lastHead; i < arrayIndex; i++) {
			visited[Hypergraph.array[i]] = false;
		}
		totVol += arrayIndex;
		
		System.out.println("clique-reduced graph, total vol " + totVol);
		System.out.println("clique-reduced graph, maximum degree " + maxDegree);
	}
	
	public static int[] getCardinality() {
		
		int[] cardinalities = new int[Hypergraph.getEdgeSize()];
		int[] cardinality_sort = new int[Hypergraph.getEdgeSize()];
		for (int i = 0; i < Hypergraph.getEdgeSize(); i++) {
			cardinalities[i] = Hypergraph.getCardinality(i);
		}
		
		Arrays.sort(cardinalities);
		
		for (int i = 0; i < cardinalities.length; i++) {
			cardinality_sort[i] = cardinalities[cardinalities.length - 1 - i];
		}
		
		return cardinality_sort;
	}
	
	public static int[] perCardinalityCnt() {
		
		int maxCardinality = -1;
		for (int i = 0; i < Hypergraph.getEdgeSize(); i++) {
			if (maxCardinality < Hypergraph.getCardinality(i)) maxCardinality = Hypergraph.getCardinality(i);
		}
		
		int[] cnt = new int[maxCardinality + 1];
		for (int i = 0; i < Hypergraph.getEdgeSize(); i++) {
			cnt[Hypergraph.getCardinality(i)]++;
		}
		
		return cnt;
	}
	
	public static int[] getDegree() {
		
		int n = Hypergraph.getVertexSize();
		int endN = n - 1;
		
		int[] degrees = new int[n];
		int[] degrees_sort = new int[n];
		
		for (int node = 0; node < n; node++) {
			// for each incident edge
			if (node == endN) degrees[node] = Hypergraph.INC_eID.length - Hypergraph.INC_head[node];
			else degrees[node] = Hypergraph.INC_head[node + 1] - Hypergraph.INC_head[node];
		}
		Arrays.sort(degrees);
		
		for (int i = 0; i < n; i++) degrees_sort[i] = degrees[degrees.length - 1 - i];
		
		return degrees_sort;
	}
	
	public static void saveCardinality() {
		int[] cardinality_sort = new int[Hypergraph.getEdgeSize()];
		cardinality_sort = getCardinality();
		String fileOutputPre = FilePath_Mon.filePathPre + "/stat/";
		
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "cardinality_" + Hypergraph.dataset + ".txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "cardinality_" + Hypergraph.dataset + "_connect.txt");
			}
			
			for (int i = 0; i < cardinality_sort.length; i++) {
				fwCount.write(cardinality_sort[i] + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		int[] cnt = perCardinalityCnt();
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) {
				fwCount = new FileWriter(fileOutputPre + "perCardinalityCnt_" + Hypergraph.dataset + ".txt");
			} else {
				fwCount = new FileWriter(fileOutputPre + "perCardinalityCnt_" + Hypergraph.dataset + "_connect.txt");
			}
			
			for (int i = 0; i < cnt.length; i++) {
				fwCount.write(cnt[i] + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void saveDegree() {
		
		int[] degrees_sort = new int[Hypergraph.getVertexSize()];
		degrees_sort = getDegree();
		
		String fileOutputPre = FilePath_Mon.filePathPre + "/stat/";
		try {
			FileWriter fwCount = null;
			if (!Constant.CONNECTED) fwCount = new FileWriter(fileOutputPre + "degree_" + Hypergraph.dataset + ".txt");
			else fwCount = new FileWriter(fileOutputPre + "degree_" + Hypergraph.dataset + "_connect.txt");
			
			for (int i = 0; i < degrees_sort.length; i++) fwCount.write(degrees_sort[i] + "\n");
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
//		basicInfo();
//		saveCardinality();
		saveDegree();
	}
}
