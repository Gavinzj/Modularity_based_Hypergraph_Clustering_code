package preprocess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import utilities.Constant;
import utilities.FilePath_Mon;

public class newID {
	
	private static String fileInput;
	private static String fileOutput_vID_Idx;
	private static String fileOutput_pro;
	private static String fileOutput_edge;
	private static String fileOutput_inc;
	private static String fileOutput_bipartite;
	
	private static void pro_1() throws IOException {
		
		String[] strs = FilePath_Mon.filePathPre.split("/");
		
		if (!Constant.CONNECTED) {
			fileInput = FilePath_Mon.filePathPre + "/" + strs[strs.length - 1] + ".txt";
			fileOutput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx.txt";
			fileOutput_pro = FilePath_Mon.filePathPre + "/pro.txt";
			fileOutput_edge = FilePath_Mon.filePathPre + "/edge.txt";
			fileOutput_inc = FilePath_Mon.filePathPre + "/inc.txt";
			fileOutput_bipartite = FilePath_Mon.filePathPre + "/bipartite.txt";
		} else {
			fileInput = FilePath_Mon.filePathPre + "/" + strs[strs.length - 1] + "_connect.txt";
			fileOutput_vID_Idx = FilePath_Mon.filePathPre + "/vID_Idx_connect.txt";
			fileOutput_pro = FilePath_Mon.filePathPre + "/pro_connect.txt";
			fileOutput_edge = FilePath_Mon.filePathPre + "/edge_connect.txt";
			fileOutput_inc = FilePath_Mon.filePathPre + "/inc_connect.txt";
			fileOutput_bipartite = FilePath_Mon.filePathPre + "/bipartite_connect.txt";
		}
		
		HashMap<Integer, Integer> Vertices_ID = new HashMap<Integer, Integer>();
		int verticeID = 0;
		int totalDegree = 0;
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_pro);
			
			int skip = 0;
			String line_output, line;
			HashSet<Integer> vs_distinct;
			Iterator<Integer> it;
			int v;
			ArrayList<Integer> vs_sorted;
			try (FileReader reader = new FileReader(fileInput);
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				while ((line = bufferedReader.readLine()) != null) {
					if (skip-- > 0) continue;
					
					strs = line.split("\t");
					
					if (strs.length <= 1) continue;
					
					vs_distinct = new HashSet<Integer>(strs.length);
					for (int i = 0; i < strs.length; i++) {
						vs_distinct.add(Integer.parseInt(strs[i]));
					}
					
					it = vs_distinct.iterator();
					vs_sorted = new ArrayList<Integer>();
			        while (it.hasNext()) {
			        	v = it.next();
			        	
			            if (!Vertices_ID.containsKey(v)) {
							Vertices_ID.put(v, verticeID++);
						}
			            
			            v = Vertices_ID.get(v);
			            vs_sorted.add(v);
			        }
			        
			        totalDegree += vs_sorted.size();
			        Collections.sort(vs_sorted);
			        
			        line_output = vs_sorted.get(0) + "";
			        for (int i = 1; i < vs_sorted.size(); i++) {
			        	line_output = line_output + "\t" + vs_sorted.get(i);
			        }
			        
			        fwCount.write(line_output + "\n");
				}
			}
			
			fwCount.close();
			
			System.out.println("pro_1 vertices size: " + Vertices_ID.size() + " verticeID " + verticeID + " total degree " + totalDegree);
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_vID_Idx);
			for (int vID : Vertices_ID.keySet()) {
				fwCount.write(vID + "\t" + Vertices_ID.get(vID) + "\n");
			}
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		// collect garbage
		Vertices_ID = null;
		Runtime.getRuntime().gc();
	}
	
	private static void pro_2() {
		
		HashSet<String> hyperedges = new HashSet<String>();
		int edgeID = 0;
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput_edge);
			
			String line;
			try (FileReader reader = new FileReader(fileOutput_pro);
					BufferedReader bufferedReader = new BufferedReader((reader))) {
				while ((line = bufferedReader.readLine()) != null) {
					if (!hyperedges.contains(line)) {
						edgeID++;
						hyperedges.add(line);
						fwCount.write(line + "\n");
					}
				}
			}
			
			fwCount.close();
			
			System.out.println("pro_2 edge size: " + hyperedges.size() + " edgeID " + edgeID);
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		// collect garbage
		hyperedges = null;
		Runtime.getRuntime().gc();
	}
	
	private static void pro_3() throws IOException {
		int vertexSize = 0;
		int edgeSize = 0;
		int totalDegree = 0;
		
		Path path = Paths.get(fileOutput_vID_Idx);
		vertexSize = (int) Files.lines(path).count();
		
		List<Integer>[] inc = new ArrayList[vertexSize];
		List<Integer>[] biadjacency = new ArrayList[vertexSize];
		for (int i = 0; i < vertexSize; i++) {
			inc[i] = new ArrayList<Integer>();
			biadjacency[i] = new ArrayList<Integer>();
		}
		
		String line;
		String[] strs;
		int v;
		try (FileReader reader = new FileReader(fileOutput_edge);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			while ((line = bufferedReader.readLine()) != null) {
				strs = line.split("\t");
				
				for (int i = 0; i < strs.length; i++) {
					v = Integer.parseInt(strs[i]);
					inc[v].add(edgeSize);
					biadjacency[v].add((edgeSize + vertexSize));
					totalDegree++;
				}
				
				edgeSize++;
			}
		}
		
		System.out.println("pro_3 vertex size " + vertexSize + " total degree " + totalDegree);
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwAdj = new FileWriter(fileOutput_inc);
			
			String listString;
			for (int i = 0; i < vertexSize; i++) {
				Collections.sort(inc[i]);
				
				listString = inc[i].stream().map(Object::toString).collect(Collectors.joining("\t"));

				fwAdj.write(i + "\t" + listString + "\n");
			}
			
			fwAdj.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter fwAdj = new FileWriter(fileOutput_bipartite);
			
			for (int i = 0; i < vertexSize; i++) {
				Collections.sort(inc[i]);
				
				for (int j = 0; j < inc[i].size(); j++) {
					fwAdj.write(i + "\t" + inc[i].get(j) + "\n");
				}
			}
			
			fwAdj.close();

		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	private static int maxEdge(double purnRatio) throws IOException {
		
		ArrayList<Integer> sizes = new ArrayList<Integer>();
		
		int skip = 0;
		try (FileReader reader = new FileReader(fileInput);
				BufferedReader bufferedReader = new BufferedReader((reader))) {
			String line;
			while ((line = bufferedReader.readLine()) != null) {
				
				if (skip-- > 0) continue;
				
				sizes.add(line.split("\t").length);
			}
		}
		
		Collections.sort(sizes);
		
		int cutOffIdx = (int) (purnRatio * sizes.size());
		int purnEdgeSize = sizes.get(cutOffIdx);
		
		return purnEdgeSize;
	}
	
	public static void pro() throws IOException {
		pro_1();
		pro_2();
		pro_3();
	}
	
	public static void main(String arg[]) throws IOException {
		pro();
	}
}
