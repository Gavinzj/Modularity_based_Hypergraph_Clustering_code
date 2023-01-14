package statistic;

public class LossFunctions {
	
	// mean squared error https://en.wikipedia.org/wiki/Mean_squared_error
	public static double loss_MSE(double[] predicted, double[] observed) {
		double sum = 0;
		int idxs = Math.min(predicted.length, observed.length);
		
		for (int i = 0; i < idxs; i++) {
			sum += Math.pow((predicted[i] - observed[i]), 2);
		}
		
		return (sum / idxs);
	}
	
	public static double Kolmogorov_Smirnov (double[] predictedDistribution, double[] observedDistribution) {
		return 0;
	}
}
