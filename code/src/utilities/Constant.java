package utilities;

// The global constant variables used in the program. 
public class Constant {
	/** a very large value */
	public static double large = Math.pow(10, 8);
	
	/** a very small value */
	public static double small = Math.pow(10, -12);
	
	/** the largest value */
	public static double MAX = Double.MAX_VALUE;
	
	public static int UNDERFINE_VALUE = -999999;
	
	/** unit of memory measure */
	public static double MEMORY_UNIT = 1024*1024;
	
	/** unit of running time measure */
	public static double RUNNING_TIME_UNIT = 1000;
	
	/** maximum size of an array in Java */
	public static int maxArraySize = Integer.MAX_VALUE-5; // Integer.MAX_VALUE-5;
	
	/** initial edge Weight on dyadic graph */
	public static double INITIAL_EDGE_WEIGHT = 1;
	
	/** initial node weight on an edge on hypergraph */
	public static double INITIAL_NODEONEDGE_WEIGHT = 1;
	
	/** is the hypergraph connected? */
	public static boolean CONNECTED = true;
	
	/** For rounding a number with certain level of precision */
	public static double PRECISION_ENLARGE = Math.pow(10, 6);
	
	/** A variable deciding the termination of the Louvain-like programs */
	public static double EPSILON = Math.pow(10, -3);
}
