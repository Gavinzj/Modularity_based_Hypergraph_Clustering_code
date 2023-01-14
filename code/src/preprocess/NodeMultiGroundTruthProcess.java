package preprocess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.DoubleMergeSortString;
import utilities.FilePath_Mon;

public class NodeMultiGroundTruthProcess {
	
	private static String fileInput_node_groundTruth;
	private static String fileInput_vID_Idx;
	private static String fileInput_venue;
	private static String fileOutput;
	private static String fileOutput_disconnect;
	private static String fileOutput_connect;
	
	static HashMap<Integer, Integer> Vertices_ID;
//	static HashMap<String, String> venues;
	static HashMap<String, ArrayList<Integer>> labels_nodes;
	static List<String>[] nodeLabel;
	static HashMap<String, Integer> labelRanks;
	
	public static void read() throws IOException {
		
		fileInput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
		fileInput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx.txt";
		fileInput_venue = FilePath_Mon.filePathPre + "/venues.txt";
		fileOutput = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth.txt";
		fileOutput_disconnect = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_disconnect.txt";
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
		
//		venues = new HashMap<String, String>();
//		try (FileReader reader = new FileReader(fileInput_venue);
//				BufferedReader bufferedReader = new BufferedReader((reader))) {
//			String venue;
//			String venueID;
//			while ((line = bufferedReader.readLine()) != null) {
//				strs = line.split("\t");
//				venue = strs[0];
//				venueID = strs[1];
//				venues.put(venueID, venue);
//			}
//		}
		
		String label;
		labels_nodes = new HashMap<String, ArrayList<Integer>>();
		try (FileReader reader = new FileReader(fileInput_node_groundTruth);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			HashSet<String> distinctLabel;
			Iterator<String> itr;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				
				nodeID = Integer.parseInt(strs[0]);
				if (!Vertices_ID.containsKey(nodeID)) continue;
				nodeID = Vertices_ID.get(nodeID);
				
				distinctLabel = new HashSet<String>();
				for (int i = 1; i < strs.length; i++) {
					distinctLabel.add(strs[i]);
				}
				
				itr = distinctLabel.iterator();
				while (itr.hasNext()) {
					label = itr.next();
					if (!labels_nodes.containsKey(label)) labels_nodes.put(label, new ArrayList<Integer>());
					labels_nodes.get(label).add(nodeID);
				}
			}
		}
		
		sortByFrequency();
		
		nodeLabel = new ArrayList[Vertices_ID.size()];
		for (int i = 0; i < nodeLabel.length; i++) {
			nodeLabel[i] = new ArrayList<String>();
		}
		
		labels_nodes = new HashMap<String, ArrayList<Integer>>();

		int labelCnt, maxFreq, minRank, rank, minI;
		HashMap<String, Integer> distinctLabel;
		ArrayList<String> labels, maxFreqLabels;
		int[] labelFreq;
		try (FileReader reader = new FileReader(fileInput_node_groundTruth);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
						
				nodeID = Integer.parseInt(strs[0]);
				if (!Vertices_ID.containsKey(nodeID)) continue;
				nodeID = Vertices_ID.get(nodeID);
				
				labelCnt = 0;
				distinctLabel = new HashMap<String, Integer>();
				labels = new ArrayList<String>();
				for (int i = 1; i < strs.length; i++) {
					if (!distinctLabel.containsKey(strs[i])) {
						distinctLabel.put(strs[i], labelCnt++);
						labels.add(strs[i]);
					}
				}
				
				labelFreq = new int[distinctLabel.size()];
				for (int i = 1; i < strs.length; i++) labelFreq[distinctLabel.get(strs[i])]++;
				
				maxFreq = -1;
				for (int i = 0; i < labelFreq.length; i++) if (maxFreq < labelFreq[i]) maxFreq = labelFreq[i];
				
				maxFreqLabels = new ArrayList<String>();
				for (int i = 0; i < labelFreq.length; i++) {
					if (labelFreq[i] == maxFreq) {
						maxFreqLabels.add(labels.get(i));
					}
				}
				
				minRank = 9999999;
				minI = 0;
				for (int i = 1; i < strs.length; i++) {
					rank = labelRanks.get(strs[i]);
					if (minRank > rank) {
						minRank = rank;
						minI = i;
					}
				}
				label = strs[minI];
				
				if (!labels_nodes.containsKey(label)) labels_nodes.put(label, new ArrayList<Integer>());
				labels_nodes.get(label).add(nodeID);
				
				nodeLabel[nodeID].add(label);
			}
		}
		
		int labelNode = 0;
		int nonLabelNode = 0;
		for (int i = 0; i < nodeLabel.length; i++) {
			if (nodeLabel[i].size() > 0) labelNode++;
			else nonLabelNode++;
		}
		
		labels = new ArrayList<String>(labels_nodes.size());
		Iterator<String> it = labels_nodes.keySet().iterator();
        while (it.hasNext()) labels.add(it.next());
        
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
		System.out.println("# of label node: " + labelNode);
		System.out.println("# of non label node: " + nonLabelNode);
	}
	
	public static void read_connect() throws IOException {
		
		fileInput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_connect.txt";
		fileInput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx_connect.txt";
		fileOutput = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
		
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
			FileWriter fwCount = new FileWriter(fileOutput);

			String output;
			try (FileReader reader = new FileReader(fileInput_node_groundTruth);
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				while ((line = bufferedReader.readLine()) != null) {
					strs = line.split("\t");

					output = Vertices_ID.get(Integer.parseInt(strs[0])) + "";
					for (int i = 1; i < strs.length; i++) {
						nodeID = Vertices_ID.get(Integer.parseInt(strs[i]));
						output = output + "\t" + nodeID;
					}
					
					nodeCnt += strs.length;
					
					fwCount.write(output + "\n");
				}
			}

			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		System.out.println("ground truth connected " + nodeCnt);
	}
	
	public static void sortByFrequency() {
		Iterator<String> itr = labels_nodes.keySet().iterator();
		double[] srbaseArr = new double[labels_nodes.size()];
		String[] srfollowArr = new String[labels_nodes.size()];
		int cnt = 0;
		String label;
		while (itr.hasNext()) {
			label = itr.next();
			srbaseArr[cnt] = labels_nodes.get(label).size();
			srfollowArr[cnt] = label;
			cnt++;
		}
		
		DoubleMergeSortString Dsort = new DoubleMergeSortString();
		Dsort.sort(srbaseArr, srfollowArr);
		
		labelRanks = new HashMap<String, Integer>(srfollowArr.length);
		
		int top = 5;
		int rank = 0;
		for (int i = srbaseArr.length - 1; i >= 0; i--) {
			labelRanks.put(srfollowArr[i], rank++);
//			if (top-- > 0) System.out.println(venues.get(srfollowArr[i]) + " ranked " + rank + " frequency " + srbaseArr[i]);
		}
	}
	
	public static void save() throws InterruptedException, IOException {
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			
			ArrayList<Integer> cluster = new ArrayList<Integer>();
			Iterator<String> it = labels_nodes.keySet().iterator();
	        while (it.hasNext()) {
	        	String label = it.next();
	        	cluster = labels_nodes.get(label);
	        	
	        	String output = cluster.get(0) + "";
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
			FileWriter fwCount = new FileWriter(fileOutput_disconnect);
			
			ArrayList<Integer> component = new ArrayList<Integer>();
			ArrayList<Integer> cluster = new ArrayList<Integer>();
			Iterator<String> itr;
			
			for (int i = 0; i < allComponents.size(); i ++) {
				
				component = allComponents.get(i);
				labels_nodes = new HashMap<String, ArrayList<Integer>>();
				
				for (int node : component) {
					for (String l : nodeLabel[node]) {
						if (!labels_nodes.containsKey(l)) labels_nodes.put(l, new ArrayList<Integer>());
						labels_nodes.get(l).add(node);
					}
				}
				
				itr = labels_nodes.keySet().iterator();
				while (itr.hasNext()) {
		        	cluster = labels_nodes.get(itr.next());
		        	if (cluster.size() <= 0) continue;
					fwCount.write(cluster.stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
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
			
			labels_nodes = new HashMap<String, ArrayList<Integer>>();
			for (int node : component) {
				for (String l : nodeLabel[node]) {
					if (!labels_nodes.containsKey(l)) labels_nodes.put(l, new ArrayList<Integer>());
					labels_nodes.get(l).add(node);
				}
			}
			
			Iterator<String> it = labels_nodes.keySet().iterator();
			while (it.hasNext()) {
	        	cluster = labels_nodes.get(it.next());
	        	if (cluster.size() <= 0) continue;
				fwCount.write(cluster.stream().map(Object::toString).collect(Collectors.joining("\t")) + "\n");
	        }
			fwCount.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		read();
		save();
	}
}
