package graph;

public class adj2EdgeID{
	public int mID;
	public int edgeID;
	
	public adj2EdgeID(int mID, int edgeID) {
		this.mID = mID;
		this.edgeID = edgeID;
	}
	
	public String toString() {
		return mID + "";
	}
	
	@Override
	public boolean equals(Object arg0) {
		adj2EdgeID pair = (adj2EdgeID) arg0;
		
		if (this.mID != pair.mID) return false;
		
		return true;
	}
}
