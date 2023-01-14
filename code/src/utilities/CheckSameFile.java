package utilities;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Scanner;

public class CheckSameFile {
	public static String file_1 = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect.txt";
	public static String file_2 = FilePath_Mon.filePathPre + "/groundTruth/node_clustering_groundTruth_connect_0.txt";
	
	private static boolean same() throws IOException {
		System.out.println("check " + file_1);
		System.out.println("check " + file_2);
		
		Path path_1 = Paths.get(file_1);
		Scanner scanner_1 = new Scanner(path_1);
		
		Path path_2 = Paths.get(file_2);
		Scanner scanner_2 = new Scanner(path_2);
		
		while (scanner_2.hasNextLine() && scanner_2.hasNextLine()) {
			// process each line
			String line_1 = scanner_1.nextLine();
			String line_2 = scanner_2.nextLine();
			
			if (!line_1.equals(line_2)) {
				System.out.println("1");
				System.out.println(line_1);
				System.out.println(line_2);
				return false;
			}
		}
		
		if (scanner_2.hasNextLine() || scanner_2.hasNextLine()) {
			System.out.println("3");
			return false;
		}
		
		scanner_1.close();
		scanner_2.close();
		
		return true;
	}
	public static void main(String arg[]) throws IOException {
		System.out.println(same());
	}
}
