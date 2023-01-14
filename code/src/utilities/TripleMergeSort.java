package utilities;

public class TripleMergeSort {
	
	public int[] baseArr1;
	public int[] baseArr2;
	public double[] followArr;
	
	int[] base1;
	int[] base2;
	double[] follow;
	
	public TripleMergeSort() {}
	
	// refer to https://github.com/allmycode/sort/blob/master/java/AllSort.java
	public void sort(int[] srbaseArr1, int[] srbaseArr2, double[] srfollowArr) {
		baseArr1 = srbaseArr1;
		baseArr2 = srbaseArr2;
		followArr = srfollowArr;
		
		int start = 0;
		int end = baseArr1.length - 1;
		
		int mid = (start+end)/2;

		mergesortArr(start, mid);  // sort the first part
		mergesortArr(mid+1, end);  // sort the second part

		mergeArr(start, mid, end);
	}
	
	public void mergesortArr(int start, int end) {
		
		if(start < end) {
			int mid = (start+end)/2;

			mergesortArr(start, mid);  // sort the first part
			mergesortArr(mid+1, end);  // sort the second part

			mergeArr(start, mid, end);
		}
	}
	
	private void mergeArr(int start, int mid, int end) {

		int p = start, q = mid+1, k = 0;

		base1 = new int[end - start + 1];
		base2 = new int[end - start + 1];
		follow = new double[end - start + 1];

		for(int i = start; i <= end ; i++ ) {
			if(p > mid) {     //checks if first part comes to an end or not .
				base1[k] = baseArr1[q];
				base2[k] = baseArr2[q];
				follow[k] = followArr[q];
				k++;
				q++;
				
			} else if(q > end) {   //checks if second part comes to an end or not .
				base1[k] = baseArr1[p];
				base2[k] = baseArr2[p];
				follow[k] = followArr[p];
				k++;
				p++;
				
			} else if (baseArr1[p] < baseArr1[q]) {
				base1[k] = baseArr1[p];
				base2[k] = baseArr2[p];
				follow[k] = followArr[p];
				k++;
				p++;
				
			} else if (baseArr1[p] == baseArr1[q]){
				
				if (baseArr2[p] < baseArr2[q]) {
					base1[k] = baseArr1[p];
					base2[k] = baseArr2[p];
					follow[k] = followArr[p];
					k++;
					p++;
				} else {
					base1[k] = baseArr1[q];
					base2[k] = baseArr2[q];
					follow[k] = followArr[q];
					k++;
					q++;
				}
				
			} else {
				base1[k] = baseArr1[q];
				base2[k] = baseArr2[q];
				follow[k] = followArr[q];
				k++;
				q++;
			}
		}

		// updating in the originastartarray ie arr[]
		for(int i = 0; i < k; i++) {
			baseArr1[start] = base1[i];
			baseArr2[start] = base2[i];
			followArr[start] = follow[i];
			start++;
		}

	}
}
