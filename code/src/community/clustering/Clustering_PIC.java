package community.clustering;


import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.CountDownLatch;

import community.pic.PIC_noPrune_quality;
import community.pic.PIC_noPrune_speed;
import community.pic.PIC_prune12_quality;
import community.pic.PIC_prune12_speed;
import community.pic.PIC_prune1_quality;
import community.pic.PIC_prune1_speed;
import community.pic.PIC_prune2_quality;
import community.pic.PIC_prune2_speed;
import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.SingleMergeSortDouble;

public class Clustering_PIC {
	
	Cluster_PIC cluster_PIC;
	Cluster_PIC_noPrune cluster_PIC_noPrune;
	Cluster_PIC_prune1 cluster_PIC_prune1;
	Cluster_PIC_prune2 cluster_PIC_prune2;
	Cluster_PIC_prune12 cluster_PIC_prune12;
	
	CountDownLatch latch;
	SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
	
	public Clustering_PIC() {}
	
	private class Cluster_PIC extends Thread {
		private int trial;
		private String moveStrategy;
		private boolean toHigherOrder;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC(int trial, boolean toHigherOrder, String ordering, String moveStrategy, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.moveStrategy = moveStrategy;
			this.toHigherOrder = toHigherOrder;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_prune12_quality l = new PIC_prune12_quality(trial, toHigherOrder, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.louvain(moveStrategy);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	private class Cluster_PIC_noPrune extends Thread {
		private int trial;
		private String moveStrategy;
		private boolean toHigherOrder;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_noPrune(int trial, boolean toHigherOrder, String ordering, String moveStrategy, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.moveStrategy = moveStrategy;
			this.toHigherOrder = toHigherOrder;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_noPrune_quality l = new PIC_noPrune_quality(trial, toHigherOrder, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.louvain(moveStrategy);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	private class Cluster_PIC_prune1 extends Thread {
		private int trial;
		private String moveStrategy;
		private boolean toHigherOrder;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_prune1(int trial, boolean toHigherOrder, String ordering, String moveStrategy, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.moveStrategy = moveStrategy;
			this.toHigherOrder = toHigherOrder;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_prune1_quality l = new PIC_prune1_quality(trial, toHigherOrder, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.louvain(moveStrategy);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	private class Cluster_PIC_prune2 extends Thread {
		private int trial;
		private String moveStrategy;
		private boolean toHigherOrder;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_prune2(int trial, boolean toHigherOrder, String ordering, String moveStrategy, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.moveStrategy = moveStrategy;
			this.toHigherOrder = toHigherOrder;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_prune2_quality l = new PIC_prune2_quality(trial, toHigherOrder, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.louvain(moveStrategy);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	private class Cluster_PIC_prune12 extends Thread {
		private int trial;
		private String moveStrategy;
		private boolean toHigherOrder;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_prune12(int trial, boolean toHigherOrder, String ordering, String moveStrategy, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.moveStrategy = moveStrategy;
			this.toHigherOrder = toHigherOrder;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_prune12_quality l = new PIC_prune12_quality(trial, toHigherOrder, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.louvain(moveStrategy);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			latch.countDown();
		}
	}
	
	public void doClustering_PIC(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		boolean save = true;
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC = new Cluster_PIC(trial, toHigherOrder, ordering, moveStrategy, ratio, save);
		}
		latch.await();
	}
	
	public void doClustering_PIC_noPrune(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		boolean save = true;
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_noPrune = new Cluster_PIC_noPrune(trial, toHigherOrder, ordering, moveStrategy, ratio, save);
		}
		latch.await();
	}
	
	public void doClustering_PIC_prune1(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		boolean save = true;
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_prune1 = new Cluster_PIC_prune1(trial, toHigherOrder, ordering, moveStrategy, ratio, save);
		}
		latch.await();
	}
	
	public void doClustering_PIC_prune2(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		boolean save = true;
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_prune2 = new Cluster_PIC_prune2(trial, toHigherOrder, ordering, moveStrategy, ratio, save);
		}
		latch.await();
	}
	
	public void doClustering_PIC_prune12(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		boolean save = true;
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_prune12 = new Cluster_PIC_prune12(trial, toHigherOrder, ordering, moveStrategy, ratio, save);
		}
		latch.await();
	}
	
	public void runningTime_PIC(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		double[] runningTimes = new double[trials];
		double[] memories = new double[trials];
		
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune12_speed l = new PIC_prune12_speed(trial, toHigherOrder, ordering, ratio, save);
			
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
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
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
	
	public void runningTime_PIC_noPrune(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		double[] runningTimes = new double[trials];
		double[] memories = new double[trials];
		
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_noPrune_speed l = new PIC_noPrune_speed(trial, toHigherOrder, ordering, ratio, save);
			
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
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
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
	
	public void runningTime_PIC_prune1(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
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
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic1_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic1_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
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
	
	public void runningTime_PIC_prune2(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
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
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic2_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic2_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
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
	
	public void runningTime_PIC_prune12(boolean toHigherOrder, String ordering, String moveStrategy, double ratio, int trials) throws Exception {
		
		Hypergraph.array = new int[Hypergraph.maxArraySize];
		
		boolean save = false;
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		double[] runningTimes = new double[trials];
		double[] memories = new double[trials];
		
		double runningTime;
		for (int trial = 0; trial < trials; trial++) {
			
			PIC_prune12_speed l = new PIC_prune12_speed(trial, toHigherOrder, ordering, ratio, save);
			
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
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic12_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTimeMem_pic12_" + moveStrategy + 
					"_ordered_" + toHigherOrder + "_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
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
	
	
}
