/**
 * 
 */
package utilities;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

/**
 * @author feng zijin
 *
 */
public class Functions {

	public static int binomialCoe(int m) {
		int m_2;

		if (m == 0 || m == 1) {
			m_2 = 0;
		} else {
			m_2 = m * (m - 1) / 2;
		}

		return m_2;
	}

//	public static double factorial(int n) {
//		if (n == 0)
//			return 1;
//
//		return n * factorial(n - 1);
//	}

	public static double factorial(int n) {
		double res = 1;
		for (int i = 2; i <= n; i++)
			res *= i;
		return res;
	}

	public static double combination(int n, int k) {
		double nominator = 1;
		int idx = 0;
		for (int i = n; i >= 2; i--) {
			nominator *= i;
			idx++;
			if (idx >= k) break;
		}
		
		double denominator = 1;
		for (int i = 2; i <= k; i++) {
			denominator *= i;
		}
		
		return nominator / denominator;
	}

	public static String listToString(List<Integer> list) {
		String str = "";

		for (int i = 0; i < list.size() - 1; i++) {
			str += list.get(i) + ",";
		}

		if (list.size() > 0)
			str += list.get(list.size() - 1);

		return str;
	}

	public static String listToString(HashSet<Integer> list) {
		String str = "";

		Iterator<Integer> it = list.iterator();
		while (it.hasNext()) {
			str += it.next() + "\t";
		}

		return str;
	}

	public static String arrToString(String[] list) {
		String str = "";

		for (int i = 0; i < list.length - 1; i++) {
			str += list[i] + "\t";
		}

		if (list.length > 0)
			str += list[list.length - 1];

		return str;
	}

	public static String arrToString(int[] list) {
		String str = "";

		for (int i = 0; i < list.length - 1; i++) {
			str += list[i] + ",";
		}

		if (list.length > 0)
			str += list[list.length - 1];

		return str;
	}

	public static String arrToString(float[] list) {
		String str = "";

		for (int i = 0; i < list.length - 1; i++) {
			str += list[i] + "\t";
		}

		if (list.length > 0)
			str += list[list.length - 1];

		return str;
	}

	public static String arrToString(double[] list) {
		String str = "";

		for (int i = 0; i < list.length - 1; i++) {
			str += list[i] + "\t";
		}

		if (list.length > 0)
			str += list[list.length - 1];

		return str;
	}

	public static int randInt(int min, int max) {
		return ThreadLocalRandom.current().nextInt(min, max + 1);
	}

	public static float max(float[] arr) {
		return arr[arr.length - 1];
	}

	public static float min(float[] arr) {
		return arr[0];
	}
}
