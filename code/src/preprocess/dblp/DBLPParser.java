package preprocess.dblp;

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

public class DBLPParser {
	
	static String fileInput = FilePath_Mon.filePathPre + "/publications.txt";
	static String fileOutput_venues = FilePath_Mon.filePathPre + "/venues.txt";
	static String fileOutput_authors = FilePath_Mon.filePathPre + "/authors.txt";
	
	static int yearInterval = 1;
	static int minYear, maxYear;
	static ArrayList<Integer> yearsSort = new ArrayList<Integer>();
	static HashMap<String, Integer> venues = new HashMap<String, Integer>();
	static HashMap<String, Integer> authors = new HashMap<String, Integer>();
	static ArrayList<String> publications = new ArrayList<String>();
	
	public static void read() throws IOException {
		
		HashSet<Integer> years = new HashSet<Integer>();
		HashSet<String> skipTypes = new HashSet<String>();
		skipTypes.add("homepages");
		skipTypes.add("www");
		skipTypes.add("tr");
		skipTypes.add("reference");
		skipTypes.add("series");
		skipTypes.add("phd");
		
		HashSet<String> skipVenues = new HashSet<String>();
		skipVenues.add("corr");
		
		double hasfieldCnt = 0;
		int venueID = 0;
		int authorID = 0;
		int cnt = 0;
		String[] strs_t;
		String[] strs_y;
		String[] strs_a;
		String[] strs_conf;
		String conf, type, venue, line;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {

				strs_t = line.split("<<t>>");
				if (strs_t.length < 2) continue;
				
				conf = strs_t[0];
				strs_conf = conf.split("/");
				if (strs_conf.length < 3) continue;
				type = strs_conf[0];
				venue = strs_conf[1].toLowerCase();
				
				if (skipTypes.contains(type)) continue;
				if (skipVenues.contains(venue)) continue;
				
				strs_y = strs_t[1].split("<<y>>");
				if (strs_y.length < 2) continue;
				
				strs_a = strs_y[1].split("<<a>>");
				if (strs_a.length < 2) continue;
				
				int year = Integer.parseInt(strs_a[0]);
				
				if (!venues.containsKey(venue)) venues.put(venue, venueID++);
				years.add(year);
				for (int i = 1; i < strs_a.length; i ++) {
					if (!authors.containsKey(strs_a[i])) authors.put(strs_a[i], authorID++);
				}
						
				cnt++;
			}
			
			System.out.println(hasfieldCnt / cnt);
		}
		
		Iterator<Integer> itr = years.iterator();
		while (itr.hasNext()) yearsSort.add(itr.next());
		Collections.sort(yearsSort);
		
		minYear = yearsSort.get(0);
		maxYear = yearsSort.get(yearsSort.size() - 1);
		
		System.out.println("number of venues " + venues.size());
		System.out.println("number of publication " + cnt);
		System.out.println("number of authors " + authors.size());
		System.out.println("minYear " + minYear + " maxYear " + maxYear);
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_venues);
			Iterator<String> it = venues.keySet().iterator();
			while(it.hasNext()) {
				venue = it.next();
				fw_user.write(venue + "\t" + venues.get(venue) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_authors);
			String author;
			Iterator<String> it = authors.keySet().iterator();
			while(it.hasNext()) {
				author = it.next();
				fw_user.write(author + "\t" + authors.get(author) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//////////////////////////////////////////////////////////////////////////////////////
		int stopYear = 2019;
		for (int topYear = maxYear; topYear >= minYear; topYear = topYear - yearInterval) {
			int lowYear = topYear - yearInterval + 1;
			System.out.println("lowYear " + lowYear + " topYear " + topYear);
			
			publications = new ArrayList<String>();
			
			try (FileReader reader = new FileReader(fileInput);
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				while ((line = bufferedReader.readLine()) != null) {
					
					strs_t = line.split("<<t>>");
					if (strs_t.length < 2) continue;
					
					conf = strs_t[0];
					strs_conf = conf.split("/");
					if (strs_conf.length < 3) continue;
					type = strs_conf[0];
					venue = strs_conf[1].toLowerCase();
					
					if (skipTypes.contains(type)) continue;
					if (skipVenues.contains(venue)) continue;
					
					strs_y = strs_t[1].split("<<y>>");
					if (strs_y.length < 2) continue;
					
					strs_a = strs_y[1].split("<<a>>");
					if (strs_a.length < 2) continue;
					
					int year = Integer.parseInt(strs_a[0]);
					
					if (lowYear <= year && year <= topYear) {
						String publication = venues.get(venue) + "\t" + year;
						for (int i = 1; i < strs_a.length; i ++) {
							publication = publication + "\t" + authors.get(strs_a[i]);
						}
						publications.add(publication);
					}
				}
			}
			
			constructGraph(lowYear, topYear);
			
			if (lowYear <= stopYear) break;
		}
	}
	
	public static void constructGraph(int lowYear, int topYear) {
		String fileOutput_publications = FilePath_Mon.filePathPre + "/publications_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_" + lowYear + "_" + topYear + ".txt";
		
		HashSet<Integer> authors = new HashSet<Integer>();
		HashMap<Integer, List<Integer>> nodeLabels = new HashMap<Integer, List<Integer>>();
		ArrayList<String> hyperedges = new ArrayList<String>();
		
		String[] strs;
		int venue, author;
		for (String publication : publications) {
			strs = publication.split("\t");
			venue = Integer.parseInt(strs[0]);
			
//			if (year == topYear) {
			for (int i = 2; i < strs.length; i++) {
				author = Integer.parseInt(strs[i]);
				authors.add(author);
				if (!nodeLabels.containsKey(author)) nodeLabels.put(author, new ArrayList<Integer>());
				nodeLabels.get(author).add(venue);
			}
//			}
			
			String hyperedge = strs[2] + "";
			for (int i = 3; i < strs.length; i++) {
				hyperedge = hyperedge + "\t" + strs[i];
			}
			hyperedges.add(hyperedge);
		}
		
		System.out.println("authors " + authors.size() + " labeled author " + nodeLabels.keySet().size());
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_publications);
			HashSet<String> visHyperedges = new HashSet<String>();
			String hyperedge;
			for (int i = 0; i < hyperedges.size(); i++) {
				hyperedge = hyperedges.get(i);
				if (!visHyperedges.add(hyperedge)) continue;
				
				fw_user.write(hyperedges.get(i) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	////////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_node_groundTruth);
			
			String output = "";
			Iterator<Integer> itr_author = nodeLabels.keySet().iterator();
			List<Integer> itr_label;
			while(itr_author.hasNext()) {
				author = itr_author.next();
				
				output = author + "";
				itr_label = nodeLabels.get(author);
				for (int label : itr_label) {
					output = output + "\t" + label;
				}
				fw_user.write(output + "\n");
			}
			fw_user.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws Exception {
		read();
	}
}