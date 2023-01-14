package preprocess;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;
import java.util.stream.Collectors;

import hyperGraph.Hypergraph;
import utilities.FilePath_Mon;

public class CheckConnected {
	
	public static ArrayList<ArrayList<Integer>> allComponents;
	
	public static void findAllComponents() {
		HashSet<Integer> checkNodes = new HashSet<Integer>();
		
		allComponents = new ArrayList<ArrayList<Integer>>();
		
		for (int o = 0; o < Hypergraph.getVertexSize(); o++) {
			if (checkNodes.contains(o)) continue;
			
			boolean[] visited = new boolean[Hypergraph.getVertexSize()];
			Queue<Integer> queue = new LinkedList<Integer>();
			queue.add(o);
			visited[o] = true;
			
			while (!queue.isEmpty()) {
				int u = queue.remove();
				
				int second_idx = Hypergraph.getSecondIdx_INC(u);
				for (int edgeIdx = Hypergraph.INC_head[u]; edgeIdx < second_idx; edgeIdx++) {
					
					int edgeID = Hypergraph.INC_eID[edgeIdx];
					
					int second_idx_neighbor = Hypergraph.getSecondIdx_EINC(edgeID);
					for (int j = Hypergraph.EINC_head[edgeID]; j < second_idx_neighbor; j++) {
						int neighbor = Hypergraph.EINC_nID[j];
						
						if (!visited[neighbor] && !checkNodes.contains(neighbor)) {
							queue.add(neighbor);
							visited[neighbor] = true;
						}
					}
				}
			}
			
			ArrayList<Integer> nodesInComponent = new ArrayList<Integer>();
			for (int i = 0; i < visited.length; i ++) {
				if (visited[i]) {
					nodesInComponent.add(i);
					checkNodes.add(i);
				}
			}
			allComponents.add(nodesInComponent);
		}
		
		System.out.println("allComponents number " + allComponents.size());
	}
	
	public static ArrayList<Integer> getLargestComponent() {
		findAllComponents();
		
		int maxComponentSize = -1;
		int maxIndex = -1;
		
		for (int i = 0; i < allComponents.size(); i ++) {
			if (maxComponentSize < allComponents.get(i).size()) {
				maxComponentSize = allComponents.get(i).size();
				maxIndex = i;
			}
		}
		
		ArrayList<Integer> maxComponent = new ArrayList<Integer>();
		maxComponent = allComponents.get(maxIndex);
		System.out.println("maxComponentSize " + maxComponentSize);
		
		return maxComponent;
	}
	
	public static void saveAll_LargestComponent() {
		findAllComponents();
		
		// save all components
		String fileOutput_allComponent = FilePath_Mon.filePathPre + "/allComponent.txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput_allComponent);
			for (int i = 0; i < allComponents.size(); i++) {
				fwCount.write(allComponents.get(i).stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
			}
			fwCount.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// save the largest componnet
		ArrayList<Integer> maxComponent = new ArrayList<Integer>();
		
		// find largest component
		int maxComponentSize = -1;
		int maxIndex = -1;
		int componentSize;
		for (int i = 0; i < allComponents.size(); i ++) {
			componentSize = allComponents.get(i).size();
			if (maxComponentSize < componentSize) {
				maxComponentSize = componentSize;
				maxIndex = i;
			}
		}
		maxComponent = allComponents.get(maxIndex);
		
		// store the nodes that are inside the largest componnet
		HashSet<Integer> inComponent = new HashSet<Integer>(maxComponent.size());
		for (int node : maxComponent) inComponent.add(node);
		
		System.out.println("maxComponentSize " + maxComponentSize + " " + maxComponent.size() + " " + inComponent.size());
		
		int[] EINC_head = Hypergraph.EINC_head;
		int[] EINC_nID = Hypergraph.EINC_nID;
		
		// save .txt
		String fileOutput_largestComponent = FilePath_Mon.filePathPre + "/" + Hypergraph.dataset + "_connect.txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput_largestComponent);
			
			int edgeSize = Hypergraph.getEdgeSize();
			int first_idx, second_idx;
			boolean isInComponent;
			String output;
			
			// for each hyperedge
			for (int edgeID = 0; edgeID < edgeSize; edgeID++) {
				
				first_idx = EINC_head[edgeID];
				second_idx = Hypergraph.getSecondIdx_EINC(edgeID);
				isInComponent = true;
				
				// for each node in the hyperedge
				for (int i = first_idx; i < second_idx; i++) {
					if (!inComponent.contains(EINC_nID[i])) {
						isInComponent = false;
						break;
					}
				}
				
				if (isInComponent) {
					output = EINC_nID[first_idx] + "";
					for (int i = first_idx + 1; i < second_idx; i++) {
						output = output + "\t" + EINC_nID[i];
					}
					fwCount.write(output + "\n");
				}
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void saveAllcomponent() {
		findAllComponents();
		
		String fileOutput = FilePath_Mon.filePathPre + "/allComponent.txt";
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			
			for (int i = 0; i < allComponents.size(); i++) {
				fwCount.write(allComponents.get(i).stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void saveLargestComponent() {
		ArrayList<Integer> maxComponent = getLargestComponent();
		
		HashSet<Integer> inComponent = new HashSet<Integer>(maxComponent.size());
		for (int node : maxComponent) inComponent.add(node);
		
		int[] EINC_head = Hypergraph.EINC_head;
		int[] EINC_nID = Hypergraph.EINC_nID;
		
		// save dataset.txt
		try {
			FileWriter fwCount = new FileWriter(FilePath_Mon.filePathPre + "/" + Hypergraph.dataset + "_connect.txt");
			
			int edgeSize = Hypergraph.getEdgeSize();
			int first_idx, second_idx, node;
			boolean isInComponent;
			String output;
			for (int edgeID = 0; edgeID < edgeSize; edgeID++) {
				first_idx = EINC_head[edgeID];
				second_idx = Hypergraph.getSecondIdx_EINC(edgeID);
				isInComponent = true;
				
				output = EINC_nID[first_idx] + "";
				// for each node in the hyperedge
				for (int i = first_idx + 1; i < second_idx; i++) {
					node = EINC_nID[i];
					output = output + "\t" + node;
					
					if (!inComponent.contains(node)) {
						isInComponent = false;
						break;
					}
				}
				
				if (isInComponent) fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static boolean isConnected() {
		
		boolean returnBool = true;
		
		for (int o = 0; o < Hypergraph.getVertexSize(); o++) {
			System.out.println("from " + o);
			boolean[] visited = new boolean[Hypergraph.getVertexSize()];
			Queue<Integer> queue = new LinkedList<Integer>();
			queue.add(o);
			visited[0] = true;
			
			while (!queue.isEmpty()) {
				int u = queue.remove();
				
				int second_idx = Hypergraph.getSecondIdx_INC(u);
				for (int edgeIdx = Hypergraph.INC_head[u]; edgeIdx < second_idx; edgeIdx++) {
					
					int edgeID = Hypergraph.INC_eID[edgeIdx];
					
					int second_idx_neighbor = Hypergraph.getSecondIdx_EINC(edgeID);
					for (int j = Hypergraph.EINC_head[edgeID]; j < second_idx_neighbor; j++) {
						int neighbor = Hypergraph.EINC_nID[j];
						
						if (!visited[neighbor]) {
							queue.add(neighbor);
							visited[neighbor] = true;
						}
					}
				}
			}
			
			for (int i = 0; i < visited.length; i ++) {
				if (visited[i] == false) {
					return false;
				}
			}
			
			if (returnBool == true) return returnBool;
		}
		
		return returnBool;
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		Hypergraph.loadGraph();
		
//		if (!isConnected()) {
//			findLargestComponent();
//		} else {
//			System.out.println("true");
//		}
//		
//		saveAllcomponent();
		
		saveLargestComponent();
	}
}
