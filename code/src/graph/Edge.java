package graph;

public class Edge {
	public int v_1;
	public int v_2;
	public double weight;
	
	public Edge(int v_1, int v_2, double weight) {
		this.v_1 = v_1;
		this.v_2 = v_2;
		this.weight = weight;
	}
	
	public Edge(int v_1, int v_2) {
		this.v_1 = v_1;
		this.v_2 = v_2;
	}
	
	public boolean isEqual(Edge another) {
		// check whether they are the same class object
		if (this.v_1 != another.v_1) return false;
		else {
			if (this.v_2 != another.v_2) return false;
			else return true;
		}
	}
	
	public String toString() {
		String string = this.v_1 + "\t" + this.v_2 + "\t" + weight;
		
		return string;
	}
}
