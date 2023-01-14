package preprocess.nbapbp;

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

public class NBAPBPParser {
	
	static String fileInput = FilePath_Mon.filePathPre + "/nba_pbp_2001-full/nba_pbp_2001.csv";
	static String fileOutput_venues = FilePath_Mon.filePathPre + "/nba_pbp_2001-full/venues.txt";
	static String fileOutput_authors = FilePath_Mon.filePathPre + "/nba_pbp_2001-full/authors.txt";
	
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
		
		HashMap<String, HashSet<Integer>> matchPlayers = new HashMap<String, HashSet<Integer>>();
		
		String[] strs;
		String line;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			
			bufferedReader.readLine();
			while ((line = bufferedReader.readLine()) != null) {
				
//				System.out.println("read " + line);
				line = line.replaceAll(",0,", ",,");
//				System.out.println("proc " + line);
				strs = line.split(",");
				int length = strs.length;
				
				int gameID = Integer.parseInt(strs[0]);
				int period = Integer.parseInt(strs[1]);
				
				int player1ID = -1;
				int player1TeamID = -1;
				if (length >= 4) {
					if (!strs[2].equals("") && !strs[3].equals("")) {
						player1ID = Integer.parseInt(strs[2]);
						player1TeamID = Integer.parseInt(strs[3]);
					}
				} else {
					continue;
				}
				
				int player2ID = -1;
				int player2TeamID = -1;
				if (length >= 6) {
					if (!strs[4].equals("") && !strs[5].equals("")) {
						player2ID = Integer.parseInt(strs[4]);
						player2TeamID = Integer.parseInt(strs[5]);
					}
				} else {
//					System.out.println(gameID + "," + period + "," + player1ID + "," + player1TeamID);
					if (!matchPlayers.containsKey(gameID + "-" + period)) {
						matchPlayers.put(gameID + "-" + period, new HashSet<Integer>());
					}
					matchPlayers.get(gameID + "-" + period).add(player1ID);
					continue;
				}
				
				int player3ID = -1;
				int player3TeamID = -1;
				if (length >= 8) {
					if (!strs[6].equals("") && !strs[7].equals("")) {
						player3ID = Integer.parseInt(strs[6]);
						player3TeamID = Integer.parseInt(strs[7]);
					}
				} else {
//					System.out.println(gameID + "," + period + "," + player1ID + "," + player1TeamID+ 
//							"," + player2ID + "," + player2TeamID);
					if (!matchPlayers.containsKey(gameID + "-" + period)) {
						matchPlayers.put(gameID + "-" + period, new HashSet<Integer>());
					}
					matchPlayers.get(gameID + "-" + period).add(player1ID);
					matchPlayers.get(gameID + "-" + period).add(player2ID);
					continue;
				}
				
//				System.out.println(gameID + "," + period + "," + player1ID + "," + player1TeamID+ 
//						"," + player2ID + "," + player2TeamID + "," + player3ID + "," + player3TeamID);
				if (!matchPlayers.containsKey(gameID + "-" + period)) {
					matchPlayers.put(gameID + "-" + period, new HashSet<Integer>());
				}
				matchPlayers.get(gameID + "-" + period).add(player1ID);
				matchPlayers.get(gameID + "-" + period).add(player2ID);
				matchPlayers.get(gameID + "-" + period).add(player3ID);
			}
		}
		
		System.out.println("# of matches " + matchPlayers.size());
		Iterator<String> itr_match = matchPlayers.keySet().iterator();
		while (itr_match.hasNext()) {
			String match = itr_match.next();
			System.out.println("match " + match + ": " + matchPlayers.get(match).size());
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