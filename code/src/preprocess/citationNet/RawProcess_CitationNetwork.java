package preprocess.citationNet;

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

// for processing the data sets available in https://www.aminer.cn/citation
public class RawProcess_CitationNetwork {
	private static String dataset = "acm9-full";
	static String fileOutput_venues = FilePath_Mon.filePathPre + "/venues.txt";
	static String fileOutput_authors = FilePath_Mon.filePathPre + "/authors.txt";
	
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	private static String fileInput = fileInputPre + "acm.txt";
	
	static int yearInterval = 1;
	static int minYear, maxYear;
	static ArrayList<Integer> yearsSort;
	
	static HashMap<String, Integer> author_IDs;
	static HashMap<String, Integer> venue_IDs;
	static List<String> allPublications;
	static ArrayList<String> publications;
	
	public static void read() throws IOException {
		
		int authorCnt = 0;
		author_IDs = new HashMap<String, Integer>();
		
		int venueCnt = 0;
		venue_IDs = new HashMap<String, Integer>();
		
		allPublications = new ArrayList<String>();
		
		HashSet<Integer> years = new HashSet<Integer>();
		
		String[] strs;
		String line, publication;
		int citedIndex;
		
		String title = "";
		String venue = "";
		int year = -1;
		int index = -1;
		List<Integer> citings = new ArrayList<Integer>();
		List<String> authors = new ArrayList<String>();
		
		int stop = 100;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			
			while ((line = bufferedReader.readLine()) != null) {
//				if (stop-- <= 0) break;
				
				line = line.toLowerCase();
//				System.out.println(line);
				
				if (line.startsWith("#*")) {
					// title
					title = line.replace("#*", "");
//					System.out.println("title " + title);
					
				} else if (line.startsWith("#@")) {
					// author
					line = line.replace("#@", "");
//					System.out.println("authors " + line);
					strs = line.split(",");
					
					for (int i = 0; i < strs.length; i++) {
						strs[i] = strs[i].trim();
						authors.add(strs[i]);
//						System.out.println("author " + strs[i]);
					}
					
				} else if (line.startsWith("#t")) {
					// year
					line = line.replace("#t", "");
					year = Integer.parseInt(line);
//					System.out.println("year " + year);
					
				} else if (line.startsWith("#c")) {
					// venue
					venue = line.replace("#c", "");
//					System.out.println("venue " + venue);
					
				} else if (line.startsWith("#index")) {
					// index
					line = line.replace("#index", "");
					index = Integer.parseInt(line);
//					System.out.println("index " + index);
					
				} else if (line.startsWith("#%")) {
					// citing
					line = line.replace("#%", "");
					citedIndex = Integer.parseInt(line);
					citings.add(citedIndex);
//					System.out.println("citing " + citing);
					
				} else if (line.startsWith("#!")) {
					// Abstract
				} else if (line.equals("")) {
//					System.out.println("empty");
					
					
					if (!title.equals("") && (authors.size() != 0) && (year != -1) && !venue.equals("") && (index != -1)) {
						
						publication = index + "";
						
						// process year
						publication = publication + ",y" + year;
						years.add(year);
						
						// process venue
						if (!venue_IDs.containsKey(venue)) venue_IDs.put(venue, venueCnt++);
						publication = publication + ",v" + venue_IDs.get(venue);
						
						// process authors
						for (String author : authors) {
							if (!author_IDs.containsKey(author)) author_IDs.put(author, authorCnt++);
							publication = publication + ",a" + author_IDs.get(author);
						}
						
						// process citing
						for (int i = 0; i < citings.size(); i++) {
							citedIndex = citings.get(i);
							publication = publication + ",c" + citedIndex;
						}
						
//						System.out.println(publication);
						allPublications.add(publication);
					}
					
					title = "";
					venue = "";
					year = -1;
					index = -1;
					citings = new ArrayList<Integer>();
					authors = new ArrayList<String>();
				}
			}
		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		
		yearsSort = new ArrayList<Integer>();
		Iterator<Integer> itr = years.iterator();
		while (itr.hasNext()) yearsSort.add(itr.next());
		Collections.sort(yearsSort);
		
		minYear = yearsSort.get(0);
		maxYear = yearsSort.get(yearsSort.size() - 1);
		
		System.out.println("number of venues " + venue_IDs.size());
		System.out.println("number of all publication " + allPublications.size());
		System.out.println("number of authors " + author_IDs.size());
		System.out.println("minYear " + minYear + " maxYear " + maxYear);
		
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
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_authors);
			String author;
			Iterator<String> it = author_IDs.keySet().iterator();
			while(it.hasNext()) {
				author = it.next();
				fw_user.write(author + "\t" + author_IDs.get(author) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//////////////////////////////////////////////////////////////////////////////////////
		
		int stopYear = 2010;
		for (int topYear = maxYear; topYear >= minYear; topYear = topYear - yearInterval) {
			int lowYear = topYear - yearInterval + 1;
			System.out.println("lowYear " + lowYear + " topYear " + topYear);
			
			publications = new ArrayList<String>();
			
			for (int i = 0; i < allPublications.size(); i++) {
				publication = allPublications.get(i);
				strs = publication.split(",");
				
				year = Integer.parseInt(strs[1].replace("y", ""));
				if (lowYear <= year && year <= topYear) {
					publications.add(publication);
				}
			}
			
			constructGraph(lowYear, topYear);
			
			if (lowYear <= stopYear) break;
		}
	}
	
	public static void constructGraph(int lowYear, int topYear) {
		
		String fileOutput_publications = FilePath_Mon.filePathPre + "/publications_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_hypergraph_cocited = FilePath_Mon.filePathPre + "/cocited_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_hypergraph_cociting = FilePath_Mon.filePathPre + "/cociting_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_author_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/author_groundTruth_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_paper_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/paper_groundTruth_" + lowYear + "_" + topYear + ".txt";
		
		HashSet<Integer> authors = new HashSet<Integer>();
		HashMap<Integer, List<String>> authorLabels = new HashMap<Integer, List<String>>();
		ArrayList<String> hyperedges = new ArrayList<String>();
		
		HashMap<Integer, String> paperLabels = new HashMap<Integer, String>();
		HashMap<Integer, HashSet<Integer>> cocited_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then B be key
		HashMap<Integer, HashSet<Integer>> cociting_papers = new HashMap<Integer, HashSet<Integer>>();	// A cited by B, then A be key
		
		List<Integer> authorsOfPublication;
		String venue, hyperedge;
		int index, year, idx, authorID, citedIndex;
		String[] strs;
		for (int i = 0; i < publications.size(); i++) {
			
//			System.out.println(publications.get(i));
			
			strs = publications.get(i).split(",");
			
			index = Integer.parseInt(strs[0]);
			
			year = Integer.parseInt(strs[1].replace("y", ""));
			
			venue = strs[2].replace("v", "");
			
			if (!paperLabels.containsKey(index)) paperLabels.put(index, venue);
			
			idx = 3;
			authorsOfPublication = new ArrayList<Integer>();
			while (idx < strs.length && strs[idx].startsWith("a")) {
				authorID = Integer.parseInt(strs[idx].replace("a", ""));
				authorsOfPublication.add(authorID);
				authors.add(authorID);
				
				if (!authorLabels.containsKey(authorID)) authorLabels.put(authorID, new ArrayList<String>());
				authorLabels.get(authorID).add(venue);
				
//				System.out.println("author " + strs[idx].replace("a", ""));
				idx++;
			}
			
			while (idx < strs.length && strs[idx].startsWith("c")) {
				citedIndex = Integer.parseInt(strs[idx].replace("c", ""));
					
				if (!cocited_papers.containsKey(index)) cocited_papers.put(index, new HashSet<Integer>());
				cocited_papers.get(index).add(citedIndex);
				
				if (!cociting_papers.containsKey(citedIndex)) cociting_papers.put(citedIndex, new HashSet<Integer>());
				cociting_papers.get(citedIndex).add(index);
				
//				System.out.println("cites " + strs[idx].replace("c", ""));
				idx++;
			}
			
			///////////////////////////////////////////////////////////////////////////
			
			int minHyperedgeSize = 2;
			if (authorsOfPublication.size() >= minHyperedgeSize) {
				hyperedge = authorsOfPublication.get(0) + "";
				for (int j = 1; j < authorsOfPublication.size(); j++) {
					hyperedge = hyperedge + "\t" + authorsOfPublication.get(j);
				}
				hyperedges.add(hyperedge);
			}
		}
		
		System.out.println("authors " + authors.size() + " labeled author " + authorLabels.keySet().size());
		System.out.println("# of cocited " + cocited_papers.size());
		System.out.println("# of cociting_papers " + cociting_papers.size());
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		int counter = 0;
		HashSet<String> visHyperedges = new HashSet<String>();
		try {
			FileWriter fw_user = new FileWriter(fileOutput_publications);
			for (int i = 0; i < hyperedges.size(); i++) {
				hyperedge = hyperedges.get(i);
				if (!visHyperedges.add(hyperedge)) continue;
				
				fw_user.write(hyperedges.get(i) + "\n");
				counter++;
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("# of coauthor hyperedges " + counter);
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// save hyperedges for cited
		int minHyperedgeSize = 2;
		visHyperedges = new HashSet<String>(cocited_papers.size());
		
		counter = 0;
		List<Integer> nodes;
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
				
		////////////////////////////////////////////////////////////////////////////////////////////
		
		// save author ground truth
		try {
			FileWriter fw_user = new FileWriter(fileOutput_author_groundTruth);
			
			String output = "";
			Iterator<Integer> itr_author = authorLabels.keySet().iterator();
			List<String> itr_label;
			while(itr_author.hasNext()) {
				authorID = itr_author.next();
				
				output = authorID + "";
				itr_label = authorLabels.get(authorID);
				for (String label : itr_label) {
					output = output + "\t" + label;
				}
				fw_user.write(output + "\n");
			}
			fw_user.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		// save paper ground truth
		try {
			FileWriter fw_user = new FileWriter(fileOutput_paper_groundTruth);
			
			String output = "";
			Iterator<Integer> itr_paper = paperLabels.keySet().iterator();
			int paperID;
			while(itr_paper.hasNext()) {
				paperID = itr_paper.next();
				
				if (!paperLabels.containsKey(paperID)) continue;
				
				output = paperID + "\t" + paperLabels.get(paperID);
				
				fw_user.write(output + "\n");
			}
			fw_user.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
	}
}
