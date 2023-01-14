package preprocess.others;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import utilities.FilePath_Mon;

// for processing the data sets cora and citeseer available in https://linqs.soe.ucsc.edu/data
public class RawProcess_Citation {
	private static String dataset = "cora-full";
	private static String fileOutput_hypergraph_word = FilePath_Mon.filePathPre + "/cora_word.txt";
	private static String fileOutput_hypergraph_cocited = FilePath_Mon.filePathPre + "/cora_cocited.txt";
	private static String fileOutput_hypergraph_cociting = FilePath_Mon.filePathPre + "/cora_cociting.txt";
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	private static String fileInput_content = fileInputPre + "cora_content.txt";
	private static String fileInput_cites = fileInputPre + "cora_cites.txt";
	
	static HashMap<String, Integer> nodes_IDs;
	static HashSet<String> labels;
	static HashMap<Integer, String> nodes_Labels;
	
	static HashMap<String, HashSet<Integer>> wordSet_papers;
	static HashMap<Integer, HashSet<Integer>> cocited_papers;
	static HashMap<Integer, HashSet<Integer>> cociting_papers;
	
	public static void read() throws IOException {
		
		int nodeCnt = 0;
		nodes_IDs = new HashMap<String, Integer>();
		labels = new HashSet<String>();
		nodes_Labels = new HashMap<Integer, String>();
		wordSet_papers = new HashMap<String, HashSet<Integer>>();
		
		int stop = 2;
		try (FileReader reader = new FileReader(fileInput_content);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, paper, label, words;
			int node;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
//				if (stop-- <= 0) continue;
				
//				System.out.println(line);
				paper = strs[0];
				label = strs[strs.length - 1];
				
				words = "";
				for (int word = 1; word < strs.length - 1; word++) {
					if (strs[word].equals("0")) continue;
					words = words + word + ",";
//					System.out.println("word " + word + " " + strs[word]);
				}
//				System.out.println("paper " + paper + " label " + label + " words " + words);
				
				if (!nodes_IDs.containsKey(paper)) nodes_IDs.put(paper, nodeCnt++);
				node = nodes_IDs.get(paper);
				
				labels.add(label);
				
				if (!nodes_Labels.containsKey(node)) nodes_Labels.put(node, label);
				
				if (!wordSet_papers.containsKey(words)) {
					wordSet_papers.put(words, new HashSet<Integer>());
				}
				wordSet_papers.get(words).add(node);
			}
		}
		
		System.out.println("# of nodes " + nodes_IDs.size() + " " + nodeCnt + " " + nodes_Labels.size());
		System.out.println("# of labels " + labels.size());
		System.out.println("# of word sets " + wordSet_papers.size());
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		
		cocited_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then B be key
		cociting_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then A be key
		
		try (FileReader reader = new FileReader(fileInput_cites);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, citedPaper, citingPaper;
			int citedNode, citingNode;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
//				if (stop-- <= 0) continue;
				
//				System.out.println(line);
				citedPaper = strs[0];
				citingPaper = strs[1];
				
				if (nodes_IDs.containsKey(citedPaper)) citedNode = nodes_IDs.get(citedPaper);
				else continue;
				
				if (nodes_IDs.containsKey(citingPaper)) citingNode = nodes_IDs.get(citingPaper);
				else continue;
				
				if (!cocited_papers.containsKey(citingNode)) cocited_papers.put(citingNode, new HashSet<Integer>());
				cocited_papers.get(citingNode).add(citedNode);
				
				if (!cociting_papers.containsKey(citedNode)) cociting_papers.put(citedNode, new HashSet<Integer>());
				cociting_papers.get(citedNode).add(citingNode);
			}
		}
		
		System.out.println("# of cocited " + cocited_papers.size());
		System.out.println("# of cociting_papers " + cociting_papers.size());
	}
	
	public static void save() {
		// save hyperedges for words
		int minHyperedgeSize = 2;
		HashSet<String> visHyperedges = new HashSet<String>(wordSet_papers.size());
		
		int counter = 0;
		List<Integer> nodes;
		Iterator<String> itr_words;
		Iterator<Integer> itr_nodes;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_word);
			
			String words, output;
			itr_words = wordSet_papers.keySet().iterator();
			while (itr_words.hasNext()) {
				words = itr_words.next();
				
				nodes = new ArrayList<Integer>();
				
				itr_nodes = wordSet_papers.get(words).iterator();
				while(itr_nodes.hasNext()) {
					nodes.add(itr_nodes.next());
				}
				
				if (nodes.size() < minHyperedgeSize) continue;
				
				Collections.sort(nodes);
				
				output = nodes.get(0) + "";
				for (int i = 1; i < nodes.size(); i ++) {
					output = output + "\t" + nodes.get(i);
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
				counter++;
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of cowords hyperedges " + counter);
		
		///////////////////////////////////////////////////////////////////////////////////////////
		
		// save hyperedges for cited
		minHyperedgeSize = 2;
		visHyperedges = new HashSet<String>(cocited_papers.size());
		
		counter = 0;
		Iterator<Integer> itr_citingNodes;
		Iterator<Integer> itr_citedNodes;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_cocited);
			
			int citingNode;
			String output;
			itr_citingNodes = cocited_papers.keySet().iterator();
			
			while (itr_citingNodes.hasNext()) {
				citingNode = itr_citingNodes.next();
				
				nodes = new ArrayList<Integer>();
				
				itr_citedNodes = cocited_papers.get(citingNode).iterator();
				while(itr_citedNodes.hasNext()) {
					nodes.add(itr_citedNodes.next());
				}
				
				if (nodes.size() < minHyperedgeSize) continue;
				
				Collections.sort(nodes);
				
				output = nodes.get(0) + "";
				for (int i = 1; i < nodes.size(); i ++) {
					output = output + "\t" + nodes.get(i);
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
				counter++;
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of cocited hyperedges " + counter);
				
		///////////////////////////////////////////////////////////////////////////////////////////////
		
		// save hyperedges for citing
		minHyperedgeSize = 2;
		visHyperedges = new HashSet<String>(cociting_papers.size());
		
		counter = 0;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_cociting);
			
			int citedNode;
			String output;
			itr_citedNodes = cociting_papers.keySet().iterator();
			
			while (itr_citedNodes.hasNext()) {
				citedNode = itr_citedNodes.next();
				
				nodes = new ArrayList<Integer>();
				
				itr_citingNodes = cociting_papers.get(citedNode).iterator();
				while(itr_citingNodes.hasNext()) {
					nodes.add(itr_citingNodes.next());
				}
				
				if (nodes.size() < minHyperedgeSize) continue;
				
				Collections.sort(nodes);
				
				output = nodes.get(0) + "";
				for (int i = 1; i < nodes.size(); i ++) {
					output = output + "\t" + nodes.get(i);
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
				counter++;
				fwCount.write(output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of cociting hyperedges " + counter);
				
		////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//save node_groundTruth
		try {
			FileWriter fwCount = new FileWriter(fileOutput_node_groundTruth);
			
			itr_nodes = nodes_Labels.keySet().iterator();
			int node;
			String label;
			while(itr_nodes.hasNext()) {
				node = itr_nodes.next();
				label = nodes_Labels.get(node);
				fwCount.write(node + "\t" + label + "\n");
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
