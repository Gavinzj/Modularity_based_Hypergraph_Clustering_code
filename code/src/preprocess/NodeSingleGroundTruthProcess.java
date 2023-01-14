package preprocess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;

public class NodeSingleGroundTruthProcess {

	private static String fileInput_node_groundTruth;
	private static String fileInput_vID_Idx;
	private static String fileOutput_clustering;
	private static String fileOutput_clustering_disconnect;
	private static String fileOutput_connect;

	static HashMap<Integer, Integer> Vertices_ID;
	static HashMap<String, ArrayList<Integer>> labels_nodes;
	static String[] nodeLabel;

	public static void read() throws IOException {

		fileInput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
		fileInput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx.txt";
		fileOutput_clustering = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth.txt";
		fileOutput_clustering_disconnect = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
		fileOutput_connect = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_connect.txt";
		
		Vertices_ID = new HashMap<Integer, Integer>();

		String[] strs;
		String line;
		int vertex, nodeID;
		try (FileReader reader = new FileReader(fileInput_vID_Idx);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");

				vertex = Integer.parseInt(strs[0]);
				nodeID = Integer.parseInt(strs[1]);

				Vertices_ID.put(vertex, nodeID);
			}
		}

		nodeLabel = new String[Vertices_ID.size()];
		labels_nodes = new HashMap<String, ArrayList<Integer>>();
		String label;
		int nodeCnt = 0;
		try (FileReader reader = new FileReader(fileInput_node_groundTruth);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");

				if (strs.length == 2) {
					nodeID = Integer.parseInt(strs[0]);
					label = strs[1];
				} else if (strs.length == 1) {
					nodeID = nodeCnt++;
					label = strs[0];
				} else {
					System.out.println("something wrong_NodeGroundTruthProcess_read");
					continue;
				}

				if (!Vertices_ID.containsKey(nodeID)) continue;
				nodeID = Vertices_ID.get(nodeID);
				
				if (!labels_nodes.containsKey(label)) labels_nodes.put(label, new ArrayList<Integer>());

				labels_nodes.get(label).add(nodeID);
				nodeLabel[nodeID] = label;
			}
		}

		ArrayList<String> labels = new ArrayList<String>(labels_nodes.size());
		Iterator<String> itr = labels_nodes.keySet().iterator();
		while (itr.hasNext()) {
			label = itr.next();
			labels.add(label);
		}

		int nodeSize = 0;
		for (int i = 0; i < labels.size(); i++) {
			label = labels.get(i);
			if (labels_nodes.get(label).size() == 0) {
				labels_nodes.remove(label);
				continue;
			}
			nodeSize += labels_nodes.get(label).size();
		}

		System.out.println("# of node label: " + labels_nodes.size());
		System.out.println("# of node: " + nodeSize);
	}

	public static void read_connect() throws IOException {

		fileInput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_connect.txt";
		fileInput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx_connect.txt";
		fileOutput_clustering = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
		Vertices_ID = new HashMap<Integer, Integer>();

		String[] strs;
		String line;
		int vertex, nodeID;
		try (FileReader reader = new FileReader(fileInput_vID_Idx);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				vertex = Integer.parseInt(strs[0]);
				nodeID = Integer.parseInt(strs[1]);
				Vertices_ID.put(vertex, nodeID);
			}
		}

		int nodeCnt = 0;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_clustering);

			String output;
			try (FileReader reader = new FileReader(fileInput_node_groundTruth);
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				while ((line = bufferedReader.readLine()) != null) {
					strs = line.split("\t");
					nodeCnt += strs.length;
					
					output = Vertices_ID.get(Integer.parseInt(strs[0])) + "";
					for (int i = 1; i < strs.length; i++) {
						nodeID = Vertices_ID.get(Integer.parseInt(strs[i]));
						output = output + "\t" + nodeID;
					}
					fwCount.write(output + "\n");
				}
			}

			fwCount.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("ground truth connected node count " + nodeCnt);
	}

	public static void save() throws InterruptedException, IOException {

		try {
			FileWriter fwCount = new FileWriter(fileOutput_clustering);

			ArrayList<Integer> cluster = new ArrayList<Integer>();
			String label, output;
			
			Iterator<String> itr = labels_nodes.keySet().iterator();
			while (itr.hasNext()) {
				label = itr.next();
				cluster = labels_nodes.get(label);

				output = cluster.get(0) + "";
				for (int i = 1; i < cluster.size(); i++) {
					output = output + "\t" + cluster.get(i);
				}

				fwCount.write(output + "\n");
			}

			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}

		//////////////////////////////////////////////////////////////////////////

		Constant.CONNECTED = false;
		Hypergraph.loadGraph();

		CheckConnected.findAllComponents();
		ArrayList<ArrayList<Integer>> allComponents = CheckConnected.allComponents;
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_clustering_disconnect);

			ArrayList<Integer> component = new ArrayList<Integer>();
			ArrayList<Integer> cluster = new ArrayList<Integer>();
			Iterator<String> itr;
			String label, output;
			for (int i = 0; i < allComponents.size(); i++) {

				component = allComponents.get(i);
				labels_nodes = new HashMap<String, ArrayList<Integer>>();

				// for each node in component
				for (int node : component) {
					label = nodeLabel[node];
					if (!labels_nodes.containsKey(label)) labels_nodes.put(label, new ArrayList<Integer>());
					labels_nodes.get(label).add(node);
				}

				itr = labels_nodes.keySet().iterator();
				while (itr.hasNext()) {
					label = itr.next();
					cluster = labels_nodes.get(label);
					
					output = cluster.get(0) + "";
					for (int j = 1; j < cluster.size(); j++) {
						output = output + "\t" + cluster.get(j);
					}
					fwCount.write(output + "\n");
				}
			}

			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}

		//////////////////////////////////////////////////////////////////////////

		try {
			FileWriter fwCount = new FileWriter(fileOutput_connect);

			ArrayList<Integer> cluster = new ArrayList<Integer>();
			ArrayList<Integer> component = new ArrayList<Integer>();
			
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
			component = allComponents.get(maxIndex);
			System.out.println("maxComponentSize " + maxComponentSize + " " + component.size());
			
			String label;
			labels_nodes = new HashMap<String, ArrayList<Integer>>();
			for (int node : component) {
				label = nodeLabel[node];
				if (!labels_nodes.containsKey(label)) labels_nodes.put(label, new ArrayList<Integer>());
				labels_nodes.get(label).add(node);
			}

			String output;
			Iterator<String> it = labels_nodes.keySet().iterator();
			while (it.hasNext()) {
				label = it.next();
				cluster = labels_nodes.get(label);
				
				output = cluster.get(0) + "";
				for (int j = 1; j < cluster.size(); j++) {
					output = output + "\t" + cluster.get(j);
				}
				fwCount.write(output + "\n");
			}

			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}

	public static void main(String arg[]) throws IOException, InterruptedException {
//		read();
//		save();
	}
}
