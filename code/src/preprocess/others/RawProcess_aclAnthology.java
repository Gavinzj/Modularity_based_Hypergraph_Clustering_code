package preprocess.others;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import utilities.FilePath_Mon;

// for processing the data sets available in https://aclanthology.org/
public class RawProcess_aclAnthology {
	static String dataset = "aclAnthology20-full";
	static String fileOutput_hypergraph = FilePath_Mon.filePathPre + "/aclAnthology20.txt";
	static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	static String fileOutput_venues = FilePath_Mon.filePathPre + "/venues.txt";
	static String fileOutput_authors = FilePath_Mon.filePathPre + "/authors.txt";
	
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	
	static HashMap<Integer, HashSet<Integer>> author_Labels;
	static HashMap<String, Integer> author_IDs;
	static HashMap<String, Integer> venue_IDs;
	static List<String> publications;
	
	public static void read() throws IOException {
		
		int authorCnt = 0;
		int authorID;
		author_IDs = new HashMap<String, Integer>();
		author_Labels = new HashMap<Integer, HashSet<Integer>>();
		
		int venueCnt = 0;
		int venueID;
		venue_IDs = new HashMap<String, Integer>();
	
		String line, venue, author;
		publications = new ArrayList<String>();
		
		File path = new File(fileInputPre);

	    File [] files = path.listFiles();
	    for (int i = 0; i < files.length; i++){
	        if (files[i].isFile()){ //this line weeds out other directories/folders
	            
//	        	System.out.println();
//	        	System.out.println(files[i]);
	        	String fileName = files[i].getName();
	        	venue = fileName.replace(".bib", "");
	        	if (!venue_IDs.containsKey(venue)) venue_IDs.put(venue, venueCnt++);
	        	venueID = venue_IDs.get(venue);
	        	
	            int stop = 30;
	    		try (FileReader reader = new FileReader(fileInputPre + fileName);
	    				BufferedReader bufferedReader = new BufferedReader((reader))) {
	    			
	    			String publication = "";
	    			boolean startOfAuthor = false;
	    			while ((line = bufferedReader.readLine()) != null) {
//	    				if (stop-- <= 0) break;
	    				
	    				line = line.toLowerCase();
//	    				System.out.println(line);
	    				
	    				if (line.contains("@inproceedings{")) {
	    					
//	    					System.out.println("add:" + publication);
	    					publications.add(publication);
	    					
	    					publication = venueID + "";
	    				}
	    				
	    				if (line.contains("author = \"")) {
	    					startOfAuthor = true;
	    					line = line.replace("author = \"", "");
	    					line = line.replace("  and", "");
	    					line = line.replace("\",", "");
	    					author = line.trim();
	    					
	    					if (!author_IDs.containsKey(author)) author_IDs.put(author, authorCnt++);
	    					authorID = author_IDs.get(author);
	    					
//	    					System.out.println("author:" + author + " " + authorID);
	    					
	    					if (!author_Labels.containsKey(authorID)) author_Labels.put(authorID, new HashSet<Integer>());
	    					author_Labels.get(authorID).add(venueID);
	    					
	    					publication = publication + "\t" + authorID;
	    					
	    					continue;
	    				}
	    				
	    				if (startOfAuthor && !line.contains(" = \"")) {
	    					
	    					if (line.endsWith("  and")) {
		    					line = line.replace("  and", "");
		    					author = line.trim();
		    					
		    					if (!author_IDs.containsKey(author)) author_IDs.put(author, authorCnt++);
		    					authorID = author_IDs.get(author);
		    					
//		    					System.out.println("author:" + author + " " + authorID);
		    					
		    					if (!author_Labels.containsKey(authorID)) author_Labels.put(authorID, new HashSet<Integer>());
		    					author_Labels.get(authorID).add(venueID);
		    					
		    					publication = publication + "\t" + authorID;
		    					
	    					} else if (line.endsWith("\",")) {
	    						line = line.replace("\",", "");
		    					author = line.trim();
		    					
		    					if (!author_IDs.containsKey(author)) author_IDs.put(author, authorCnt++);
		    					authorID = author_IDs.get(author);
		    					
//		    					System.out.println("author:" + author + " " + authorID);
		    					
		    					if (!author_Labels.containsKey(authorID)) author_Labels.put(authorID, new HashSet<Integer>());
		    					author_Labels.get(authorID).add(venueID);
		    					
		    					publication = publication + "\t" + authorID;
		    					
	    					} else {
	    						System.out.println("wrong");
	    					}
	    					
	    				} else {
	    					startOfAuthor = false;
	    				}
	    			}
	    		}
	        }
	    }
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////

		System.out.println("number of venues " + venue_IDs.size());
		System.out.println("number of all publication " + publications.size());
		System.out.println("number of authors " + author_IDs.size());
	}
	
	public static void save() {
		
		// save hyperedges for words
		int minHyperedgeSize = 2;
		HashSet<String> visHyperedges = new HashSet<String>(publications.size());
		
		int counter = 0;
		try {
			FileWriter fwCount = new FileWriter(fileOutput_hypergraph);
			
			String[] strs;
			String output;
			for (String publication : publications) {
				strs = publication.split("\t");
				
				if ((strs.length - 1) < minHyperedgeSize) continue;
				
				output = strs[1];
				for (int i = 2; i < strs.length; i++) {
					output = output + "\t" + strs[i];
				}
				
				if (visHyperedges.contains(output)) continue;
				visHyperedges.add(output);
				
//				System.out.println(publication);
//				System.out.println(output);
//				System.out.println();
				
				fwCount.write(output + "\n");
				
				counter++;
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		System.out.println("# of hyperedges " + counter);
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		String venue;
		try {
			FileWriter fw_user = new FileWriter(fileOutput_venues);
			Iterator<String> it = venue_IDs.keySet().iterator();
			while(it.hasNext()) {
				venue = it.next();
				fw_user.write(venue + "\t" + venue_IDs.get(venue) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		String author;
		try {
			FileWriter fw_user = new FileWriter(fileOutput_authors);
			Iterator<String> it = author_IDs.keySet().iterator();
			while(it.hasNext()) {
				author = it.next();
				fw_user.write(author + "\t" + author_IDs.get(author) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		//save node_groundTruth
		try {
			FileWriter fwCount = new FileWriter(fileOutput_node_groundTruth);
			
			Iterator<Integer> itr_authorID = author_Labels.keySet().iterator();
			Iterator<Integer> itr_labels;
			int authorID, labelID;
			String output;
			while(itr_authorID.hasNext()) {
				authorID = itr_authorID.next();
				itr_labels = author_Labels.get(authorID).iterator();
				
				output = authorID + "";
				while (itr_labels.hasNext()) {
					labelID = itr_labels.next();
					output = output + "\t" + labelID;
//					System.out.print(labelID + "\t");
				}
				
//				System.out.println();
//				System.out.println(output);
//				System.out.println();
				
				fwCount.write(authorID + "\t" + output + "\n");
			}
			
			fwCount.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		
//		String fileOutput_publications = FilePath_Mon.filePathPre + "/publications_" + lowYear + "_" + topYear + ".txt";
//		String fileOutput_hypergraph_cocited = FilePath_Mon.filePathPre + "/cocited_" + lowYear + "_" + topYear + ".txt";
//		String fileOutput_hypergraph_cociting = FilePath_Mon.filePathPre + "/cociting_" + lowYear + "_" + topYear + ".txt";
//		String fileOutput_author_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/author_groundTruth_" + lowYear + "_" + topYear + ".txt";
//		String fileOutput_paper_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/paper_groundTruth_" + lowYear + "_" + topYear + ".txt";
//		
//		HashSet<Integer> authors = new HashSet<Integer>();
//		HashMap<Integer, List<String>> authorLabels = new HashMap<Integer, List<String>>();
//		ArrayList<String> hyperedges = new ArrayList<String>();
//		
//		HashMap<Integer, String> paperLabels = new HashMap<Integer, String>();
//		HashMap<Integer, HashSet<Integer>> cocited_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then B be key
//		HashMap<Integer, HashSet<Integer>> cociting_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then A be key
//		
//		List<Integer> authorsOfPublication;
//		String venue, hyperedge;
//		int index, year, idx, authorID, citedIndex;
//		String[] strs;
//		for (int i = 0; i < publications.size(); i++) {
//			
////			System.out.println(publications.get(i));
//			
//			strs = publications.get(i).split(",");
//			
//			index = Integer.parseInt(strs[0]);
//			
//			year = Integer.parseInt(strs[1].replace("y", ""));
//			
//			venue = strs[2].replace("v", "");
//			
//			if (!paperLabels.containsKey(index)) paperLabels.put(index, venue);
//			
//			idx = 3;
//			authorsOfPublication = new ArrayList<Integer>();
//			while (idx < strs.length && strs[idx].startsWith("a")) {
//				authorID = Integer.parseInt(strs[idx].replace("a", ""));
//				authorsOfPublication.add(authorID);
//				authors.add(authorID);
//				
//				if (!authorLabels.containsKey(authorID)) authorLabels.put(authorID, new ArrayList<String>());
//				authorLabels.get(authorID).add(venue);
//				
////				System.out.println("author " + strs[idx].replace("a", ""));
//				idx++;
//			}
//			
//			while (idx < strs.length && strs[idx].startsWith("c")) {
//				citedIndex = Integer.parseInt(strs[idx].replace("c", ""));
//					
//				if (!cocited_papers.containsKey(index)) cocited_papers.put(index, new HashSet<Integer>());
//				cocited_papers.get(index).add(citedIndex);
//				
//				if (!cociting_papers.containsKey(citedIndex)) cociting_papers.put(citedIndex, new HashSet<Integer>());
//				cociting_papers.get(citedIndex).add(index);
//				
////				System.out.println("cites " + strs[idx].replace("c", ""));
//				idx++;
//			}
//			
//			///////////////////////////////////////////////////////////////////////////
//			
//			int minHyperedgeSize = 2;
//			if (authorsOfPublication.size() >= minHyperedgeSize) {
//				hyperedge = authorsOfPublication.get(0) + "";
//				for (int j = 1; j < authorsOfPublication.size(); j++) {
//					hyperedge = hyperedge + "\t" + authorsOfPublication.get(j);
//				}
//				hyperedges.add(hyperedge);
//			}
//		}
//		
//		System.out.println("authors " + authors.size() + " labeled author " + authorLabels.keySet().size());
//		System.out.println("# of cocited " + cocited_papers.size());
//		System.out.println("# of cociting_papers " + cociting_papers.size());
//		
//		/////////////////////////////////////////////////////////////////////////////////////////////////////
//		
//		int counter = 0;
//		HashSet<String> visHyperedges = new HashSet<String>();
//		try {
//			FileWriter fw_user = new FileWriter(fileOutput_publications);
//			for (int i = 0; i < hyperedges.size(); i++) {
//				hyperedge = hyperedges.get(i);
//				if (!visHyperedges.add(hyperedge)) continue;
//				
//				fw_user.write(hyperedges.get(i) + "\n");
//				counter++;
//			}
//			fw_user.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//		System.out.println("# of coauthor hyperedges " + counter);
//		
//		/////////////////////////////////////////////////////////////////////////////////////////////////////
//		
//		// save hyperedges for cited
//		int minHyperedgeSize = 2;
//		visHyperedges = new HashSet<String>(cocited_papers.size());
//		
//		counter = 0;
//		List<Integer> nodes;
//		Iterator<Integer> itr_citingNodes;
//		Iterator<Integer> itr_citedNodes;
//		try {
//			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_cocited);
//			
//			int citingNode;
//			String output;
//			itr_citingNodes = cocited_papers.keySet().iterator();
//			
//			while (itr_citingNodes.hasNext()) {
//				citingNode = itr_citingNodes.next();
//				
//				nodes = new ArrayList<Integer>();
//				
//				itr_citedNodes = cocited_papers.get(citingNode).iterator();
//				while(itr_citedNodes.hasNext()) {
//					nodes.add(itr_citedNodes.next());
//				}
//				
//				if (nodes.size() < minHyperedgeSize) continue;
//				
//				Collections.sort(nodes);
//				
//				output = nodes.get(0) + "";
//				for (int i = 1; i < nodes.size(); i ++) {
//					output = output + "\t" + nodes.get(i);
//				}
//				
//				if (visHyperedges.contains(output)) continue;
//				visHyperedges.add(output);
//				
//				counter++;
//				fwCount.write(output + "\n");
//			}
//			
//			fwCount.close();
//			
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			return;
//		}
//		
//		System.out.println("# of cocited hyperedges " + counter);
//				
//		///////////////////////////////////////////////////////////////////////////////////////////////
//		
//		// save hyperedges for citing
//		minHyperedgeSize = 2;
//		visHyperedges = new HashSet<String>(cociting_papers.size());
//		
//		counter = 0;
//		try {
//			FileWriter fwCount = new FileWriter(fileOutput_hypergraph_cociting);
//			
//			int citedNode;
//			String output;
//			itr_citedNodes = cociting_papers.keySet().iterator();
//			
//			while (itr_citedNodes.hasNext()) {
//				citedNode = itr_citedNodes.next();
//				
//				nodes = new ArrayList<Integer>();
//				
//				itr_citingNodes = cociting_papers.get(citedNode).iterator();
//				while(itr_citingNodes.hasNext()) {
//					nodes.add(itr_citingNodes.next());
//				}
//				
//				if (nodes.size() < minHyperedgeSize) continue;
//				
//				Collections.sort(nodes);
//				
//				output = nodes.get(0) + "";
//				for (int i = 1; i < nodes.size(); i ++) {
//					output = output + "\t" + nodes.get(i);
//				}
//				
//				if (visHyperedges.contains(output)) continue;
//				visHyperedges.add(output);
//				
//				counter++;
//				fwCount.write(output + "\n");
//			}
//			
//			fwCount.close();
//			
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			return;
//		}
//		
//		System.out.println("# of cociting hyperedges " + counter);
//				
//		////////////////////////////////////////////////////////////////////////////////////////////
//		
//		// save author ground truth
//		try {
//			FileWriter fw_user = new FileWriter(fileOutput_author_groundTruth);
//			
//			String output = "";
//			Iterator<Integer> itr_author = authorLabels.keySet().iterator();
//			List<String> itr_label;
//			while(itr_author.hasNext()) {
//				authorID = itr_author.next();
//				
//				output = authorID + "";
//				itr_label = authorLabels.get(authorID);
//				for (String label : itr_label) {
//					output = output + "\t" + label;
//				}
//				fw_user.write(output + "\n");
//			}
//			fw_user.close();
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//		
//		// save paper ground truth
//		try {
//			FileWriter fw_user = new FileWriter(fileOutput_paper_groundTruth);
//			
//			String output = "";
//			Iterator<Integer> itr_paper = paperLabels.keySet().iterator();
//			int paperID;
//			while(itr_paper.hasNext()) {
//				paperID = itr_paper.next();
//				
//				if (!paperLabels.containsKey(paperID)) continue;
//				
//				output = paperID + "\t" + paperLabels.get(paperID);
//				
//				fw_user.write(output + "\n");
//			}
//			fw_user.close();
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
		save();
	}
}
