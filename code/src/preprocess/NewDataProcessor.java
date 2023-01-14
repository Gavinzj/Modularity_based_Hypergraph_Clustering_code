package preprocess;

import java.io.IOException;

import hyperGraph.Hypergraph;
import statistic.Statistic_Hypergraph;
import utilities.Constant;

public class NewDataProcessor {
	
	public static void processMultiSeries() throws IOException, InterruptedException {
		
		Constant.CONNECTED = false;
		newID.pro();
		
		NodeMultiGroundTruthProcess.read();
		NodeMultiGroundTruthProcess.save();
		
		Hypergraph.loadGraph();
		Statistic_Hypergraph.basicInfo();
		Statistic_Hypergraph.saveCardinality();
		
		CheckConnected.saveAll_LargestComponent();
		
		Constant.CONNECTED = true;
		newID.pro();
		
		NodeMultiGroundTruthProcess.read_connect();
		
		Hypergraph.loadGraph();
		Statistic_Hypergraph.basicInfo();
		Statistic_Hypergraph.saveCardinality();
	}
	
	public static void processSingleSeries() throws IOException, InterruptedException {
		Constant.CONNECTED = false;
		newID.pro();
		
		NodeSingleGroundTruthProcess.read();
		NodeSingleGroundTruthProcess.save();
		
		Hypergraph.loadGraph();
		Statistic_Hypergraph.basicInfo();
		Statistic_Hypergraph.saveCardinality();
		
		CheckConnected.saveAll_LargestComponent();
		
		Constant.CONNECTED = true;
		newID.pro();
		
		NodeSingleGroundTruthProcess.read_connect();
		
		Hypergraph.loadGraph();
		Statistic_Hypergraph.basicInfo();
		Statistic_Hypergraph.saveCardinality();
	}
	
	public static void main(String arg[]) throws Exception {
		processMultiSeries();
		
//		processSingleSeries();
	}
}
