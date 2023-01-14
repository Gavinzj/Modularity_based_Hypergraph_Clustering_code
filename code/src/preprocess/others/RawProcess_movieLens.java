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

// for processing the data sets movieLens available in https://grouplens.org/datasets/movielens/
public class RawProcess_movieLens {
	private static String dataset = "movieLens-full";
	private static String fileOutput_hypergraph = FilePath_Mon.filePathPre + "/movieLens.txt";
	private static String fileOutput_node_groundTruth = FilePath_Mon.filePathPre + "/groundTruth/node_groundTruth.txt";
	
	private static String fileInputPre = FilePath_Mon.filePathPre + "/" + dataset + "/";
	private static String fileInput_metaData = fileInputPre + "metadata_updated.json";
	private static String fileInput_labels = fileInputPre + "movies.csv";
	
	static HashMap<String, Integer> movies_IDs;
	static HashMap<String, Integer> labels_IDs;
	static HashMap<Integer, HashSet<Integer>> movies_Labels;
	static HashMap<String, HashSet<Integer>> director_movieIDs;
	
	public static void read() throws IOException {
		
		int movieCnt = 0;
		movies_IDs = new HashMap<String, Integer>();
		int labelCnt = 0;
		labels_IDs = new HashMap<String, Integer>();
		
		movies_Labels = new HashMap<Integer, HashSet<Integer>>();
		
		director_movieIDs = new HashMap<String, HashSet<Integer>>();
		
		int stop = 2;
		try (FileReader reader = new FileReader(fileInput_labels);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, movie, label;
			int movieID, labelID;
			String[] strs, labelStrs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split(",");
//				if (stop-- <= 0) continue;
				
//				System.out.println(line);
				movie = strs[1];
				labelStrs = strs[2].split("\\|");
				
//				System.out.println(movie);
				
				if (!movies_IDs.containsKey(movie)) movies_IDs.put(movie, movieCnt++);
				movieID = movies_IDs.get(movie);
				
				for (int i = 0; i < labelStrs.length; i++) {
					label = labelStrs[i];
//					System.out.println(label);
					if (!labels_IDs.containsKey(label)) labels_IDs.put(label, labelCnt++);
				}
				
				for (int i = 0; i < labelStrs.length; i++) {
					labelID = labels_IDs.get(labelStrs[i]);
					
					if (!movies_Labels.containsKey(movieID)) movies_Labels.put(movieID, new HashSet<Integer>());
					movies_Labels.get(movieID).add(labelID);
				}
			}
		}
		
		System.out.println("# of movies " + movies_IDs.size() + " " + movieCnt);
		System.out.println("# of labels " + labels_IDs.size() + " " + labelCnt);
		
		stop = 100;
		try (FileReader reader = new FileReader(fileInput_metaData);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line, movie, director;
			int movieID;
			String[] strs;
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split(",");
//				if (stop-- <= 0) continue;
				
				movie = strs[0].trim().replace("{\"title\": \"", "");
				movie = movie.replace("\"", "");
				if (!movies_IDs.containsKey(movie)) continue;
				movieID = movies_IDs.get(movie);
				
				director = strs[1].trim();
				director = director.replace("\"directedBy\": \"", "");
				director = director.replace("\"", "");
				
				if (director.equals("") || director.equals(" ")) continue;
				
//				System.out.println(line);
//				System.out.println(movie + " " + movieID);
//				System.out.println(director);
				
				if (!director_movieIDs.containsKey(director)) director_movieIDs.put(director, new HashSet<Integer>());
				director_movieIDs.get(director).add(movieID);
			}
		}
	}
	
	
	public static void save() {
		
		// save hyperedge
		try {
			FileWriter fw_user = new FileWriter(fileOutput_hypergraph);
			HashSet<String> visHyperedges = new HashSet<String>();
			int minCardinality = 2;
			
			String director, hyperedge;
			List<Integer> movieIDs;
			Iterator<String> itr_directors = director_movieIDs.keySet().iterator();
			Iterator<Integer> itr_movieIDs;
			while (itr_directors.hasNext()) {
				director = itr_directors.next();
				itr_movieIDs = director_movieIDs.get(director).iterator();
				
				movieIDs = new ArrayList<Integer>();
				while (itr_movieIDs.hasNext()) {
					movieIDs.add(itr_movieIDs.next());
				}
				Collections.sort(movieIDs);
				
				if (movieIDs.size() < minCardinality) continue;
				
				hyperedge = movieIDs.get(0) + "";
				for (int i = 1; i < movieIDs.size(); i++) {
					hyperedge = hyperedge + "\t" + movieIDs.get(i);
				}
				
				if (!visHyperedges.add(hyperedge)) continue;
				
//				System.out.println("director " + director + ": " + director_movieIDs.get(director).size());
//				System.out.println(movieIDs);
//				System.out.println(hyperedge);
				
				fw_user.write(hyperedge + "\n");
			}
			
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		// save ground truth
		try {
			FileWriter fw_user = new FileWriter(fileOutput_node_groundTruth);
			
			List<Integer> labels;
			Iterator<Integer> itr_labels;
			Iterator<Integer> itr_movieIDs = movies_Labels.keySet().iterator();
			String output;
			int movieID;
			
			while (itr_movieIDs.hasNext()) {
				movieID = itr_movieIDs.next();
				
				labels = new ArrayList<Integer>();
				
				itr_labels = movies_Labels.get(movieID).iterator();
				while (itr_labels.hasNext()) {
					labels.add(itr_labels.next());
				}
				
				Collections.sort(labels);
				
				output = movieID + "";
				for (int i = 0; i < labels.size(); i++) {
					output = output + "\t" + labels.get(i);
				}
				
//				System.out.println("movie " + movieID + ": " + movies_Labels.get(movieID));
//				System.out.println(labels);
//				System.out.println(output);
				
				fw_user.write(output + "\n");
			}
			
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws IOException {
		read();
		save();
	}
}
