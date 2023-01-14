package utilities;

public class Constant {
	/** a very large value */
	public static double large = Math.pow(10, 8);
	
	/** a very small value */
	public static double small = Math.pow(10, -12);
	
	/** the largest value */
	public static double MAX = Double.MAX_VALUE;
	
	public static int UNDERFINE_VALUE = -999999;
	
	// unit of memory measure
	public static double MEMORY_UNIT = 1024*1024;
	
	// unit of running time measure
	public static double RUNNING_TIME_UNIT = 1000;
	
	// initial edge Weight on Graph
	public static double INITIAL_EDGE_WEIGHT = 1;
	
	// initial node weight on edge on Hypergraph
	public static double INITIAL_NODEONEDGE_WEIGHT = 1;
	
	public static boolean CONNECTED = true;
	
	public static double PRECISION_ENLARGE = Math.pow(10, 6);
	
	public static double EPSILON = Math.pow(10, -3);
}
