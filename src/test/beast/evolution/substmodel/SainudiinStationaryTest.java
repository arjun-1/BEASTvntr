package test.beast.evolution.substmodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.Sainudiin;
import junit.framework.TestCase;


@Description("Test Sainudiin matrix exponentiation")
public class SainudiinStationaryTest extends TestCase {

	public interface Instance {
		Double[] getPi();

		Double getRb();
		Double getIeq();
		Double getG();
		Double getOneOnA1();

		double getDistance();

		double[] getExpectedResult();
	}


	protected Instance test0 = new Instance() {
		@Override
		public Double[] getPi() {
			return new Double[]{0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667,0.066666667};
		}

		@Override
		public Double getRb() {
			return 0.5;
		}
		@Override
		public Double getIeq() {
			return 5.5;
		}
		@Override
		public Double getG() {
			return 0.1;
		}
		@Override
		public Double getOneOnA1() {
			return 1.0;
		}
		
		@Override
		public double getDistance() {
			return 0.5;
		}

		@Override
		public double[] getExpectedResult() {
			return new double[]{
			   0.19065711,0.12877908,0.11435306,0.10634927,0.09821227,0.087937400,0.075374914,0.061382671,0.047266394,0.034306577,0.023420854,0.015018162,0.009036781,0.0050976913,0.0028077555};
		}
	};

	Instance[] all = {test0};

	public void testSainudiinStationary() throws Exception {
		for (Instance test : all) {

			RealParameter f = new RealParameter(test.getPi());

			Frequencies freqs = new Frequencies();
			freqs.initByName("frequencies", f, "estimate", false);


			Sainudiin sainudiin = new Sainudiin();
			sainudiin.setNrOfStates(15);
			sainudiin.setMinRepeat(0);
			sainudiin.initByName("rb", test.getRb().toString(),
				"ieq", test.getIeq().toString(), 
				"g", test.getG().toString(), 
				"oneOnA1", test.getOneOnA1().toString()
				,"frequencies", freqs);

			double distance = test.getDistance();

			double[] mat = new double[15 * 15];
			sainudiin.getTransitionProbabilities(null, distance, 0, 1, mat);

			final double[] result = test.getExpectedResult();

			double[] Eval = sainudiin.eigenDecomposition.getEigenValues();
			double[] Ievc = sainudiin.eigenDecomposition.getInverseEigenVectors();

			double[] stationaryDistribution = sainudiin.findStationaryDistribution(Eval, Ievc);

			for (int k = 0; k < 15; ++k) {
				assertEquals(stationaryDistribution[k], result[k], 1e-8);
				System.out.println(k + " : " + (stationaryDistribution[k] - result[k]));
			}
		}
	}
}