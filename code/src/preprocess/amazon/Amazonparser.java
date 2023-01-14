package preprocess.amazon;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import utilities.FilePath_Mon;

public class Amazonparser {
	
	static String dataset = "/raw/Arts_Crafts_and_Sewing";
	static String fileInput = FilePath_Mon.filePathPre + dataset + "/reviews.txt";
	static String fileInput_brand = FilePath_Mon.filePathPre + dataset + "/brands.txt";
	
	static HashMap<String, Integer> asinIDs = new HashMap<String, Integer>();
	static HashMap<String, String> brands = new HashMap<String, String>();
	
	static ArrayList<String> outputs = new ArrayList<String>();
	
	public static void read() throws IOException {
		
		String line, asin, brand;
		String[] strs;
		int asinID = 0;
		try (FileReader reader = new FileReader(fileInput_brand);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				
				asin = strs[0];
				brand = strs[1];
				
				if (!asinIDs.containsKey(asin)) {
					asinIDs.put(asin, asinID++);
					brands.put(asin, brand);
				}
			}
		}
		
		
		System.out.println("number of asins " + asinIDs.size() + " " + brands.size());
		
		int skip = 10;
		
		HashSet<Integer> nodes;
		String output;
		Iterator<Integer> itr;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				
//				if (skip-- < 0) break;
				
				strs = line.split("\t");
				if (strs.length <= 2) continue;
				
//				System.out.println(line);
				
				nodes = new HashSet<Integer>();
				for (int i = 1; i < strs.length; i++) {
					asin = strs[i];
					if (!asinIDs.containsKey(asin)) continue;
					
					asinID = asinIDs.get(asin);
					nodes.add(asinID);
					
//					System.out.println(asin + " " + asinID);
				}
				
				if (nodes.size() < 2) continue;
				
				itr = nodes.iterator();
				output = itr.next() + "";
				while(itr.hasNext()) {
					output = output + "\t" + itr.next();
				}
//				System.out.println(output + "\n");
				outputs.add(output);
			}
		}
	}
	
	public static void save() {
		String fileOutput_reviews = FilePath_Mon.filePathPre + dataset + "/reviewsPro.txt";
		String fileOutput_brands = FilePath_Mon.filePathPre + dataset + "/node_groundTruth.txt";
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_reviews);
			for (String output : outputs) {
				fw_user.write(output + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_brands);
			
			String asin;
			Iterator<String> itr = asinIDs.keySet().iterator();
			while (itr.hasNext()) {
				asin = itr.next();
				fw_user.write(asinIDs.get(asin) + "\t" + brands.get(asin) + "\n");
			}
			fw_user.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String arg[]) throws Exception {
		read();
		save();
	}
}