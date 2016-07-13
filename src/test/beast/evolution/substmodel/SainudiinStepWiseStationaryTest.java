package test.beast.evolution.substmodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.SainudiinStepWise;
import beast.evolution.substitutionmodel.Frequencies;
import junit.framework.TestCase;
import beast.evolution.substitutionmodel.EigenSystem;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.DefaultEigenSystem;


@Description("Test SainudiinStepWise stationary distribution")
public class SainudiinStepWiseStationaryTest extends TestCase {

	public interface Instance {
		Double getRb();
		Double getIeq();
		Double getOneOnA1();

		double getDistance();

		double[] getExpectedResult();
	}


	protected Instance test0 = new Instance() {
		@Override
		public Double getRb() {
			return 0.5;
		}
		@Override
		public Double getIeq() {
			return 5.5;
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
				0.15831999,0.12259051,0.11595949,0.11306057,0.10773526,0.09798356,0.083985909,0.067334535,0.050249143,0.034787847,0.022289183,0.013193332,0.0072046799,0.0036258243,0.0016801706
			};
		}
	};

	Instance[] all = {test0};

	public void testSainudiinStepWiseStationary() throws Exception {
		for (Instance test : all) {
			SainudiinStepWise sainudiinstepwise = new SainudiinStepWise();
			sainudiinstepwise.setNrOfStates(15);
			sainudiinstepwise.setMinRepeat(0);
			sainudiinstepwise.initByName("rb", test.getRb().toString(),
				"ieq", test.getIeq().toString(),
				"oneOnA1", test.getOneOnA1().toString());

			double[] mat = new double[15 * 15];
			sainudiinstepwise.getTransitionProbabilities(null, 1.0, 0, 1, mat);

			final double[] result = test.getExpectedResult();

			double[] frequencies = sainudiinstepwise.getFrequencies();
			for (int k = 0; k < 15; ++k) {
				assertEquals(frequencies[k], result[k], 1e-8);
				System.out.println(k + " : " + (frequencies[k] - result[k]));
			}

			double[] Eval = sainudiinstepwise.eigenDecomposition.getEigenValues();
			double[] Ievc = sainudiinstepwise.eigenDecomposition.getInverseEigenVectors();

			double[] stationaryDistribution = sainudiinstepwise.getStationaryDistribution(Eval, Ievc);

			for (int k = 0; k < 15; ++k) {
				assertEquals(stationaryDistribution[k], result[k], 1e-8);
				System.out.println(k + " : " + (stationaryDistribution[k] - result[k]));
			}
			
		}
	}
}