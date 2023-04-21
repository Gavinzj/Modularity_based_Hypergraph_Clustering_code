package exp.scalability;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import community.pic.PIC_noPrune_speed;
import community.pic.PIC_prune12_speed;
import community.pic.PIC_prune1_speed;
import community.pic.PIC_prune2_speed;
import hyperGraph.Hypergraph;
import statistic.Statistic_Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.SingleMergeSortDouble;

public class Scalability {
	
	public static void runningTime_onLambda_PIC_noPrune (String ordering, double ratio, int trials, double lambda) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_noPrune_speed l = new PIC_noPrune_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_lambda_" + String.format("%.1f",lambda) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onLambda_PIC_prune1 (String ordering, double ratio, int trials, double lambda) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune1_speed l = new PIC_prune1_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_lambda_" + String.format("%.1f",lambda) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onLambda_PIC_prune2 (String ordering, double ratio, int trials, double lambda) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune2_speed l = new PIC_prune2_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_lambda_" + String.format("%.1f",lambda) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onLambda_PIC_prune12 (String ordering, double ratio, int trials, double lambda) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune12_speed l = new PIC_prune12_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_lambda_" + String.format("%.1f",lambda) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onNodeSize_PIC_noPrune (String ordering, double ratio, int trials, double divideRatio) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_noPrune_speed l = new PIC_noPrune_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_divide_" + String.format("%.1f",divideRatio) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
		
	public static void runningTime_onNodeSize_PIC_prune1 (String ordering, double ratio, int trials, double divideRatio) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune1_speed l = new PIC_prune1_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_divide_" + String.format("%.1f",divideRatio) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
		
	public static void runningTime_onNodeSize_PIC_prune2 (String ordering, double ratio, int trials, double divideRatio) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune2_speed l = new PIC_prune2_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_divide_" + String.format("%.1f",divideRatio) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void runningTime_onNodeSize_PIC_prune12 (String ordering, double ratio, int trials, double divideRatio) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;
		
		double[] runningTimes = new double[trials];
		
		// run each trial at a time using a single core
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune12_speed l = new PIC_prune12_speed(trial, ordering, ratio, save);
			
			l.initPart1();
			
			Hypergraph.garbbageCollector.gc();
			
			runningTime = l.initPart2();
			runningTime += Double.parseDouble(l.pic());
			
			System.out.println("runningTime " + runningTime);
			runningTimes[trial] = runningTime;
		}
		
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		avgRunningTime -= (maxRunningTime + minRunningTime);
		
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		File theDir = new File(FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize");
		if (!theDir.exists()) theDir.mkdirs();
		
		String fileOutput = FilePath_Mon.filePathPre + "/scalability/runningResults_NodeSize/node_runningTime_pic_move" + 
				"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_divide_" + String.format("%.1f",divideRatio) + ".txt";
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public static void statistic_scalability() {
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		Statistic_Hypergraph.basicInfo();
	}
	
}
