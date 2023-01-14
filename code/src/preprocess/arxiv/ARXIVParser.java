package preprocess.arxiv;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import utilities.DoubleMergeSortString;
import utilities.FilePath_Mon;

public class ARXIVParser {
	
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
		
		HashSet<String> skipVenues = new HashSet<String>();
		skipVenues.add("corr");
		
		HashSet<Integer> years = new HashSet<Integer>();

		int stop = 100;
		
		int venueID = 0;
		int authorID = 0;
		int cnt = 0;
		String[] strs_t;
		String[] strs_y;
		String[] strs_a;
		String[] strs_conf;
		String line, conf, type, venue;
//		HashMap<String, Integer> venueFreqs = new HashMap<String, Integer>();
		
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {

//				if (stop-- < 0) break;
				
				strs_t = line.split("<<t>>");
				if (strs_t.length < 2) continue;
				
				conf = strs_t[0];
				strs_conf = conf.split("<<c>>");
				
				strs_y = strs_t[1].split("<<y>>");
				if (strs_y.length < 2) continue;
				
				strs_a = strs_y[1].split("<<a>>");
				if (strs_a.length < 2) continue;
				
				int year = Integer.parseInt(strs_a[0]);
				
//				int freq = 0;
				for (int i = 0; i < strs_conf.length; i++) {
					venue = strs_conf[i];
					if (skipVenues.contains(venue)) continue;
					
//					if (!venueFreqs.containsKey(venue)) {
//						venueFreqs.put(venue, 1);
//					} else {
//						freq = venueFreqs.get(venue) + 1;
//						venueFreqs.put(venue, freq);
//					}
					
					if (!venues.containsKey(venue)) venues.put(venue, venueID++);
				}
				
				years.add(year);
				
				for (int j = 1; j < strs_a.length; j ++) {
					if (!authors.containsKey(strs_a[j])) authors.put(strs_a[j], authorID++);
				}
						
				cnt++;
			}
		}
		
//		Iterator<String> itr_str = venueFreqs.keySet().iterator();
//		double[] baseArr = new double[venueFreqs.size()];
//		String[] followArr = new String[venueFreqs.size()];
//		int idx = 0;
//		while (itr_str.hasNext()) {
//			venue = itr_str.next();
//			baseArr[idx] = venueFreqs.get(venue);
//			followArr[idx] = venue;
//			idx++;
//		}
//		DoubleMergeSortString Dsort = new DoubleMergeSortString();
//		Dsort.sort(baseArr, followArr);
//		
//		for (int i = venueFreqs.size() - 1; i >= venueFreqs.size() - 10; i--) {
//			System.out.println(followArr[i] + " " + baseArr[i]);
//		}
		
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
			// TODO Auto-generated catch block
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//////////////////////////////////////////////////////////////////////////////////////
		
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
					strs_conf = conf.split("<<c>>");
					
					strs_y = strs_t[1].split("<<y>>");
					if (strs_y.length < 2) continue;
					
					strs_a = strs_y[1].split("<<a>>");
					if (strs_a.length < 2) continue;
					
					int year = Integer.parseInt(strs_a[0]);
					
					if (lowYear <= year && year <= topYear) {
						
						String publication = "";
						for (int i = 0; i < strs_conf.length; i++) {
							venue = strs_conf[i];
							if (skipVenues.contains(venue)) continue;
							
							venueID = venues.get(venue);
							publication = publication + venueID + "\t";
						}
						
						publication = publication + "y" + year;
						
						for (int i = 1; i < strs_a.length; i ++) {
							authorID = authors.get(strs_a[i]);
							publication = publication + "\t" + authorID;
						}
						publications.add(publication);
					}
				}
			}
			
			constructGraph(lowYear, topYear);
		}
	}
	
	public static void constructGraph(int lowYear, int topYear) {
		String fileOutput_publications = FilePath_Mon.filePathPre + "/publications_" + lowYear + "_" + topYear + ".txt";
		String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth_" + lowYear + "_" + topYear + ".txt";
		
		HashSet<Integer> authors = new HashSet<Integer>();
		HashMap<Integer, HashSet<Integer>> nodeLabels = new HashMap<Integer, HashSet<Integer>>();
		ArrayList<String> hyperedges = new ArrayList<String>();
		
		String[] strs;
		int year, author;
		ArrayList<Integer> venues;
		for (String publication : publications) {
			
			strs = publication.split("\t");
			
			venues = new ArrayList<Integer>();
			year = -1;
			
			int idx = -1;
			for (int i = 0; i < strs.length; i++) {
				if (strs[i].contains("y")) {
					year = Integer.parseInt(strs[i].replace("y", ""));
					idx = i + 1;
					break;
				} else {
					venues.add(Integer.parseInt(strs[i]));
				}
			}
			
//			if (year == topYear) {
			for (int i = idx; i < strs.length; i++) {
				author = Integer.parseInt(strs[i]);
				authors.add(author);
				if (!nodeLabels.containsKey(author)) nodeLabels.put(author, new HashSet<Integer>());
				
				for (int venue : venues) {
					nodeLabels.get(author).add(venue);
				}
			}
//			}
			
			String hyperedge = strs[idx] + "";
			for (int i = idx+1; i < strs.length; i++) {
				hyperedge = hyperedge + "\t" + strs[i];
			}
			hyperedges.add(hyperedge);
		}
		
		System.out.println("authors " + authors.size() + " labeled author " + nodeLabels.keySet().size());
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_publications);
			
			for (int i = 0; i < hyperedges.size(); i++) {
				fw_user.write(hyperedges.get(i) + "\n");
			}
			
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_node_groundTruth);
			
			String output = "";
			Iterator<Integer> itr_author = nodeLabels.keySet().iterator();
			Iterator<Integer> itr_label;
			while(itr_author.hasNext()) {
				author = itr_author.next();
				
				output = author + "";
				itr_label = nodeLabels.get(author).iterator();
				while(itr_label.hasNext()) {
					output = output + "\t" + itr_label.next();
				}
				
				fw_user.write(output + "\n");
			}
			
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws Exception {
		read();
	}
}