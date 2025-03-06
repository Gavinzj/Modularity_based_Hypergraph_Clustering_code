package community.clustering;


import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.CountDownLatch;

import community.pic.PIC_noPrune_quality;
import community.pic.PIC_noPrune_speed;
import community.pic.PIC_prune12_quality;
import community.pic.PIC_prune12_speed;
import hyperGraph.Hypergraph;
import utilities.Constant;
import utilities.FilePath_Mon;
import utilities.SingleMergeSortDouble;

public class Clustering_PIC {
	
	Cluster_PIC_noPrune cluster_PIC_noPrune;
	Cluster_PIC_prune12 cluster_PIC_prune12;
	
	CountDownLatch latch;
	SingleMergeSortDouble SSorter = new SingleMergeSortDouble();
	
	public Clustering_PIC() {}
	
	private class Cluster_PIC_noPrune extends Thread {
		private int trial;
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_noPrune(int trial, String ordering, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_noPrune_quality l = new PIC_noPrune_quality(trial, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.pic();
				
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
		private String ordering;
		private double ratio;
		private boolean save;

		Cluster_PIC_prune12(int trial, String ordering, double ratio, boolean save) throws InterruptedException {
			this.trial = trial;
			this.ordering = ordering;
			this.ratio = ratio;
			this.save = save;
			start();
		}

		public void run() {
			
			try {
				PIC_prune12_quality l = new PIC_prune12_quality(trial, ordering, ratio, save);
				l.initPart1();
				l.initPart2();
				l.pic();
				
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
	
	public void doClustering_PIC_noPrune(String ordering, double ratio, int trials) throws Exception {
		
		boolean save = true;	// save the clustering results
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		// do clustering using PIC algorithm with no optimization technique
		// run in parallel (as we are not calculating the running time here)
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_noPrune = new Cluster_PIC_noPrune(trial, ordering, ratio, save);
		}
		latch.await();
	}
	
	public void doClustering_PIC_prune12(String ordering, double ratio, int trials) throws Exception {
		
		boolean save = true;	// save the clustering results
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
		// do clustering using PIC algorithm with no optimization technique
		// run in parallel (as we are not calculating the running time here)
		latch = new CountDownLatch(trials);
		for (int trial = 0; trial < trials; trial++) {
			cluster_PIC_prune12 = new Cluster_PIC_prune12(trial, ordering, ratio, save);
		}
		latch.await();
	}
	
	public void runningTime_PIC_noPrune(String ordering, double ratio, int trials) throws Exception {
		
		// A program to calculate the running time of PIC algorithm with no optimization technique
		// Value "trials" is the number of trials that the PIC algorithm will be run
		// When calculate the average running time of the algorithm, we will first remove
		// the largest and smallest time then calculate the average on the remaining 
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;	// the clustering results will not be saved
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
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
		
		// calculate the average running time
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		// sort the running times
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		// remove the largest and the smallest, then calculate the average
		avgRunningTime -= (maxRunningTime + minRunningTime);
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		// save to txt file
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTime_pic_move" + 
					"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTime_pic_move" + 
					"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void runningTime_PIC_prune12(String ordering, double ratio, int trials) throws Exception {
		
		// A program to calculate the running time of PIC algorithm with optimization techniques 1 and 2
		// Value "trials" is the number of trials that the PIC algorithm will be run
		// When calculate the average running time of the algorithm, we will first remove
		// the largest and smallest time then calculate the average on the remaining 
		
		Hypergraph.array = new int[Hypergraph.O_volG];
		
		boolean save = false;	// the clustering results will not be saved
		
		ratio = (double) Math.round(ratio * Constant.PRECISION_ENLARGE) / Constant.PRECISION_ENLARGE;
		
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
		
		// calculate the average running time
		double avgRunningTime = 0;
		for (int trial = 0; trial < runningTimes.length; trial++) {
			runningTime = runningTimes[trial];
			System.out.println(runningTime);
			avgRunningTime += runningTime;
		}
		
		// sort the running times
		SSorter.sort(runningTimes);
		double maxRunningTime = runningTimes[runningTimes.length - 1];
		double minRunningTime = runningTimes[0];
		
		// remove the largest and the smallest, then calculate the average
		avgRunningTime -= (maxRunningTime + minRunningTime);
		avgRunningTime /= (double) (trials - 2);
		
		System.out.println(maxRunningTime + " " + minRunningTime);
		System.out.println("average running time " + (avgRunningTime / Constant.RUNNING_TIME_UNIT));
		
		String fileOutput = "";
		if (!Constant.CONNECTED) {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTime_pic12_move" + 
					"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + ".txt";
		} else {
			fileOutput = FilePath_Mon.filePathPre + "/clustering/pic/node_runningTime_pic12_move" + 
					"_ordered_false_ordering_" + ordering + "_ratio_" + ratio + "_connect.txt";
		}
		
		try {
			FileWriter fwCount = new FileWriter(fileOutput);
			fwCount.write("runningTime," + String.format("%.8f", (avgRunningTime / Constant.RUNNING_TIME_UNIT)) + "\n");
			fwCount.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
		
}
