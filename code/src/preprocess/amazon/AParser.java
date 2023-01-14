package preprocess.amazon;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import utilities.FilePath_Mon;

public class AParser {
	
	static String dic = FilePath_Mon.filePathPre + "/raw/Arts_Crafts_and_Sewing";
	static String fileInput = dic + "/Arts_Crafts_and_Sewing.json";
	static String fileInput_meta = dic + "/meta_Arts_Crafts_and_Sewing.json";
	static String fileOutput = dic + "/reviews.txt";
	static String fileOutput_meta = dic + "/brands.txt";
	
	public static void read() throws IOException {
		
		int stop = 5;
			
		String line, reviewerID, asin, str;
		String[] strs;
		HashMap<String, ArrayList<String>> reviews = new HashMap<String, ArrayList<String>>();
		HashMap<String, String> brands = new HashMap<String, String>();
		
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
//				if (stop-- < 0) break;
//				System.out.println("====");
//				System.out.println(line);
				
				// find authors
				strs = line.split(", ");
				
				reviewerID = "";
				asin = "";
				for (int i = 0; i < strs.length; i++) {
					str = strs[i];
					if (str.contains("\"reviewerID\": ")) {
						str = str.replace("\"reviewerID\": ", "");
						reviewerID = str.replace("\"", "");
						continue;
					}
					
					if (str.contains("\"asin\": ")) {
						str = str.replace("\"asin\": ", "");
						asin = str.replace("\"", "");
						continue;
					}
				}
				
//				System.out.println("reviewerID " + reviewerID + " asin " + asin);
				
				if (reviewerID.equals("") || asin.equals("")) continue;
				
				if (!reviews.containsKey(reviewerID)) reviews.put(reviewerID, new ArrayList<String>());
				reviews.get(reviewerID).add(asin);
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		String brand;
		try (FileReader reader = new FileReader(fileInput_meta);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
//				if (stop-- < 0) break;
//				System.out.println("====");
//				System.out.println(line);
				
				// find authors
				strs = line.split(", ");
				
				brand = "";
				asin = "";
				for (int i = 0; i < strs.length; i++) {
					str = strs[i];
					if (str.contains("\"brand\": ")) {
						str = str.replace("\"brand\": ", "");
						brand = str.replace("\"", "");
						continue;
					}
					
					if (str.contains("\"asin\": ")) {
						str = str.replace("\"asin\": ", "");
						asin = str.replace("\"", "");
						continue;
					}
				}
				
//				System.out.println("brand " + brand + " asin " + asin);
				
				if (brand.equals("") || asin.equals("")) continue;
				
				brands.put(asin, brand);
			}
		}
			
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		int reviewCnt = 0;
		int brandCnt = 0;
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput);
			
			String output;
			ArrayList<String> asins;
			Iterator<String> itr = reviews.keySet().iterator();
			while(itr.hasNext()) {
				reviewerID = itr.next();
				asins = reviews.get(reviewerID);
				
				output = reviewerID + "\t" + asins.get(0);
				for (int i = 1; i < asins.size(); i++) {
					output = output + "\t" + asins.get(i);
				}
				
//				System.out.println("reviewer " + reviewerID + " has " + asins);
//				System.out.println(output);
				
				fw_user.write(output + "\n");
				reviewCnt++;
			}
			
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			FileWriter fw_user = new FileWriter(fileOutput_meta);
			
			Iterator<String> itr = brands.keySet().iterator();
			while(itr.hasNext()) {
				asin = itr.next();
				brand = brands.get(asin);
				
				fw_user.write(asin + "\t" + brand + "\n");
				brandCnt++;
			}
			
			fw_user.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("reviewCnt " + reviewCnt + " brandCnt " + brandCnt);
	}
	
	public static void main(String arg[]) throws Exception {
		read();
	}
}