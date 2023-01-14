package statistic;

import java.io.IOException;

public class RandomGenerator {
	
	// paper POWER-LAW DISTRIBUTIONS IN EMPIRICAL DATA https://arxiv.org/pdf/0706.1062.pdf
	public static double[] exp(int n, double lambda, double xmin) {
		double[] predicted = new double[n];
		
		for (int i = 0; i < n; i++) {
			predicted[i] = xmin - ((1 / lambda) * Math.log(1 - Math.random()));
		}
		
		return predicted;
	}
	
	public static void main(String arg[]) throws IOException, InterruptedException {
		double[] predicted = exp(10000, Math.log(2), 1);
		int max = -1;
		
		for (int i = 0; i < predicted.length; i++) {
			int value = (int) Math.ceil(predicted[i]);
			if (max < value) max = value;
		}
		
		int[] cnts = new int[max + 1];
		
		for (int i = 0; i < predicted.length; i++) {
			int value = (int) Math.rint(predicted[i]);
			cnts[value]++;
		}
		
		for (int i = 0; i < cnts.length; i++) {
			System.out.println(i + ";" + cnts[i]);
		}
	}
}
