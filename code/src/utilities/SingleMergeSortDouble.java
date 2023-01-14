package utilities;

public class SingleMergeSortDouble {
	
	double[] baseArr; 
	
	double[] base;
	
	public SingleMergeSortDouble() {}
	
	// refer to https://github.com/allmycode/sort/blob/master/java/AllSort.java
	public void sort(double[] srbaseArr) {
		baseArr = srbaseArr;
		
		int start = 0;
		int end = baseArr.length - 1;
		
		int mid = (start+end)/2;

		mergesort(start, mid);  // sort the first part
		mergesort(mid+1, end);  // sort the second part

		merge(start, mid, end);
	}
	
	public void mergesort(int start, int end) {
		
		if(start < end) {
			int mid = (start+end)/2;

			mergesort(start, mid);  // sort the first part
			mergesort(mid+1, end);  // sort the second part

			merge(start, mid, end);
		}
	}
	
	private void merge(int start, int mid, int end) {

		int p = start, q = mid+1, k = 0;

		base = new double[end - start + 1];

		for(int i = start; i <= end ; i++ ) {
			if(p > mid) {     //checks if first part comes to an end or not .
				base[k] = baseArr[q];
				k++;
				q++;
				
			} else if(q > end) {   //checks if second part comes to an end or not .
				base[k] = baseArr[p];
				k++;
				p++;
				
			} else if (baseArr[p] < baseArr[q]) {
				base[k] = baseArr[p];
				k++;
				p++;
				
			} else {
				base[k] = baseArr[q];
				k++;
				q++;
			}
		}

		// updating in the originastartarray ie arr[]
		for(int i = 0; i < k; i++) {
			baseArr[start] = base[i];
			start++;
		}

	}
	
}
