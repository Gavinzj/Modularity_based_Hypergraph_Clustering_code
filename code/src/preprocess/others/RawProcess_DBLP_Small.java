package preprocess.others;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;

import utilities.FilePath_Mon;

public class RawProcess_DBLP_Small {
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	
	private static String fileInput_node_label = FilePath_Mon.filePathPre + "/coauthor_dblp_small-node-labels-raw.txt";
	
	static HashSet<String> nodeLabels;
	static HashMap<Integer, String> nodes_Labels;
	
	public static void read() throws IOException {
		
		Path path = Paths.get(fileInput_node_label);
		Scanner scanner = new Scanner(path.toAbsolutePath());
		
		nodes_Labels = new HashMap<Integer, String>();
		nodeLabels = new HashSet<String>();
		
		String[] strs;
		String label;
		int nodeID = 0;
		while (scanner.hasNextLine()) {
			strs = scanner.nextLine().split(" ");
			
			label = strs[strs.length - 1];
			nodeID = Integer.parseInt(strs[0]);
			
			nodeLabels.add(label);
			nodes_Labels.put(nodeID, label);
		}
		
		scanner.close();
		
		System.out.println("# of node: " + nodes_Labels.size());
		System.out.println("# of node label: " + nodeLabels.size());
	}
	
	public static void save() {
		
		// save node class ground truth
		try {
			FileWriter fwCount = new FileWriter(fileOutput_node_groundTruth);
			
			ArrayList<Integer> nodeIDs = new ArrayList<Integer>();
			
			Iterator<Integer> it = nodes_Labels.keySet().iterator();
	        while (it.hasNext()) {
	        	nodeIDs.add(it.next());
	        }
	        
	        Collections.sort(nodeIDs);
	        
			for (int i = 0; i < nodeIDs.size(); i++) {
				fwCount.write(nodeIDs.get(i) + "\t" + nodes_Labels.get(nodeIDs.get(i)) + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
		save();
	}
}
