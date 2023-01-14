package exp.scalability;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import community.pic.PIC_prune12_speed;
import community.pic.PIC_prune1_speed;
import community.pic.PIC_prune2_speed;
import hyperGraph.Hypergraph;
import statistic.Statistic_Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.SingleMergeSortDouble;

public class Scalability {
	
	public static void runningTime_onLambda (boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials, double lambda) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		double[] memories = new double[trials];
		
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune1_speed l = new PIC_prune1_speed(trial, toHigherOrder, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			long startMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.louvain(moveStrategy));
			
			long endMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			long memoryUse = endMem - startMem;
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
			memories[trial] = memoryUse;
		}
		
		double avgRunningTime = 0;
		double avgAllMemoryUse = 0;
		double memory;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			memory = memories[trial];
			System.out.println(runningTime + "," + memory);
			avgRunningTime += runningTime;
			avgAllMemoryUse += memory;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		SSorter.sort(memories);
		double maxAllMemoryUse = memories[memories.length - 1];
		double minAllMemoryUse = memories[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		avgAllMemoryUse -= (maxAllMemoryUse + minAllMemoryUse);
		
		avgRunningTime /= (double) (trials - 2);
		avgAllMemoryUse /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println(maxAllMemoryUse + " " + minAllMemoryUse);
		System.out.println((avgRunningTime / Constant.RUNNING_TIME_UNIT) + " " + (avgAllMemoryUse / Constant.MEMORY_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults/node_runningTimeMem_Louvain_" + moveStrategy + 
				"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_lambda_" + String.format("%.1f",lambda) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.write("AllMemoryUse," + String.format("%.8f", (avgAllMemoryUse / Constant.MEMORY_UNIT)));
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onNodeSize (boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials, double divideRatio) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		double[] memories = new double[trials];
		
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune2_speed l = new PIC_prune2_speed(trial, toHigherOrder, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			long startMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.louvain(moveStrategy));
			
			long endMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			long memoryUse = endMem - startMem;
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
			memories[trial] = memoryUse;
		}
		
		double avgRunningTime = 0;
		double avgAllMemoryUse = 0;
		double memory;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			memory = memories[trial];
			System.out.println(runningTime + "," + memory);
			avgRunningTime += runningTime;
			avgAllMemoryUse += memory;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		SSorter.sort(memories);
		double maxAllMemoryUse = memories[memories.length - 1];
		double minAllMemoryUse = memories[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		avgAllMemoryUse -= (maxAllMemoryUse + minAllMemoryUse);
		
		avgRunningTime /= (double) (trials - 2);
		avgAllMemoryUse /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println(maxAllMemoryUse + " " + minAllMemoryUse);
		System.out.println((avgRunningTime / Constant.RUNNING_TIME_UNIT) + " " + (avgAllMemoryUse / Constant.MEMORY_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize/node_runningTimeMem_Louvain_" + moveStrategy + 
				"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_divide_" + String.format("%.1f",divideRatio) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.write("AllMemoryUse," + String.format("%.8f", (avgAllMemoryUse / Constant.MEMORY_UNIT)));
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
		
	public static void statistic_scalability() {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		Statistic_Hypergraph.basicInfo();
	}
	
	public static void main(String arg[]) throws Exception {
	}
}
