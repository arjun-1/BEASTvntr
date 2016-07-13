/*
* File Sainudiin.java
*
* Copyright (C) 2016 Arjun Dhawan, RIVM <arjun.dhawan@rivm.nl>
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* BEAST is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/

package beast.evolution.substitutionmodel;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.FiniteIntegerData;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.likelihood.ThreadedTreeLikelihood;
import beast.evolution.tree.Node;
import beast.core.util.Log;

@Description(
	"Substitution model of Sainudiin (R. Sainudiin et al., 2004) using\n" +
	"  Wu's modification (C. Wu and A.J. Drummond, 2011) for VNTR evolution.")
@Citation(value = 
	"Raazesh Sainudiin et al. (2004) Microsatellite Mutation Models.\n" +
	"  Genetics 168:383–395\n\n" + 
	"Chieh-Hsi Wu and  Alexei J. Drummond. (2011) Joint Inference of\n" +
	"  Microsatellite Mutation Models, Population History and Genealogies\n" + 
	"  Using Transdimensional Markov Chain Monte Carlo.\n" + 
	"  Genetics 188:151-164",
	DOI= "10.1534/genetics.103.022665", year = 2004, firstAuthorSurname = "sainudiin")

public class Sainudiin extends SubstitutionModel.Base {
	final public Input<RealParameter> rbInput = new Input<>("rb", "force of attraction to equilibrium state i_eq", Validate.REQUIRED);
	final public Input<RealParameter> ieqInput = new Input<>("ieq", "equilibrium state i_eq of mutational bias", Validate.REQUIRED);
	final public Input<RealParameter> gInput = new Input<>("g", "parameter of the geometric distribution of step sizes (1 - g = probability of a mutation being single step)", Validate.REQUIRED);
	final public Input<RealParameter> oneOnA1Input = new Input<>("oneOnA1", "proportionality of mutation rate to repeat length (i - minRepeat)", Validate.REQUIRED);

	protected EigenSystem eigenSystem;
	public EigenDecomposition eigenDecomposition;
	private EigenDecomposition storedEigenDecomposition = null;
	protected double[][] rateMatrix;

	protected boolean updateMatrix = true;
	private boolean storedUpdateMatrix = true;

	int minRepeat;
	double[] rowSum;
	double[] rowSum2;

	public void setNrOfStates(int newNrOfStates) {// required for unit test
		nrOfStates = newNrOfStates;
	}

	public void setMinRepeat(int newMinRepeat) {// required for unit test
		minRepeat = newMinRepeat;
	}

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		updateMatrix = true;

		//navigate the graph of beast objects, to find the alignment and find the nrOfStates and minRepeat.
		for (Object beastObjecti : getOutputs()) {
			if (beastObjecti instanceof SiteModel) {
				SiteModel sitemodel = (SiteModel) beastObjecti;
				for (Object beastObjectj : sitemodel.getOutputs()) {
					if (beastObjectj instanceof ThreadedTreeLikelihood) {
						ThreadedTreeLikelihood likelihood = (ThreadedTreeLikelihood) beastObjectj;
						nrOfStates = likelihood.dataInput.get().getMaxStateCount();
						FiniteIntegerData dataType = (FiniteIntegerData) likelihood.dataInput.get().getDataType();
						minRepeat = dataType.minRepeatInput.get();
						break;
					}
				}
				break;
			}
		}

		// In case initial frequencies of wrong dimension are provided in the beauti template,
		// change them.
		if(nrOfStates != 0 && nrOfStates != frequencies.getFreqs().length) {
			Log.info.println("WARNING: Frequencies has wrong size. Expected " + nrOfStates + 
				", but got " + frequencies.getFreqs().length + 
				". Will change now to correct dimension and " + 
				"assume uniform distribution for initial values.");

			String valuesString = "";
			for (int i = 0; i < nrOfStates; i++) {
				valuesString += 1 / (double) nrOfStates + " ";
			}
			RealParameter freqsRParam = new RealParameter();
			freqsRParam.setID(frequenciesInput.get().frequenciesInput.get().getID());
			freqsRParam.initByName(
							"value", valuesString,
							"lower", 0.0,
							"upper", 1.0,
							"dimension", nrOfStates
			);
			frequenciesInput.get().frequenciesInput.get().assignFrom(freqsRParam);
			frequenciesInput.get().initAndValidate();
		}

		if (nrOfStates != frequencies.getFreqs().length && nrOfStates != 0) {
				throw new IllegalArgumentException("Frequencies has wrong size. Expected " + nrOfStates + ", but got " + frequencies.getFreqs().length + ". Attempted correction failed.");
		}

		rbInput.get().setBounds(Math.max(0.0, rbInput.get().getLower()), rbInput.get().getUpper());
		ieqInput.get().setBounds(ieqInput.get().getLower(), ieqInput.get().getUpper());
		gInput.get().setBounds(Math.max(0.0, gInput.get().getLower()), Math.min(1.0, gInput.get().getUpper()));
		oneOnA1Input.get().setBounds(Math.max(0.0, oneOnA1Input.get().getLower()), oneOnA1Input.get().getUpper());
		
		eigenSystem = new DefaultEigenSystem(nrOfStates);
		rateMatrix = new double[nrOfStates][nrOfStates];
	}

	// copied from GeneralSubstitutionModel.java
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		double distance = (startTime - endTime) * rate;

		int i, j, k;
		double temp;

		// this must be synchronized to avoid being called simultaneously by
		// two different likelihood threads - AJD
		synchronized (this) {
		 if (updateMatrix) {
				setupRateMatrix();
				eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
				updateMatrix = false;
		 }
		}

		// is the following really necessary?
		// implemented a pool of iexp matrices to support multiple threads
		// without creating a new matrix each call. - AJD
		// a quick timing experiment shows no difference - RRB
		double[] iexp = new double[nrOfStates * nrOfStates];
		// Eigen vectors
		double[] Evec = eigenDecomposition.getEigenVectors();
		// inverse Eigen vectors
		double[] Ievc = eigenDecomposition.getInverseEigenVectors();
		// Eigen values
		double[] Eval = eigenDecomposition.getEigenValues();

		double stationaryDistribution[] = getStationaryDistribution(Eval, Ievc);
		double normalization = 0.0;

		for (i = 0; i < nrOfStates; i++) {
			normalization += stationaryDistribution[i] * rowSum2[i];
		}

		distance /= normalization;

		for (i = 0; i < nrOfStates; i++) {
			temp = Math.exp(distance * Eval[i]);
			for (j = 0; j < nrOfStates; j++) {
				iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
			}
		}

		int u = 0;
		for (i = 0; i < nrOfStates; i++) {
			for (j = 0; j < nrOfStates; j++) {
				temp = 0.0;
				for (k = 0; k < nrOfStates; k++) {
						temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
				}
				matrix[u] = Math.abs(temp);
				u++;
			}
		}	
	} // getTransitionProbabilities

	/*
	* access to (copy of) rate matrix *
	*/
	protected double[][] getRateMatrix() {
		return rateMatrix.clone();
	}

	public double[] getStationaryDistribution(double[] Eval, double[] Ievc) {
		//find smallest eigenvalue
		int index = 0;
		double smallest = Math.abs(Eval[0]);
		for (int i = 0; i < nrOfStates; i++) {
			if(Math.abs(Eval[i]) < smallest) {
				index = i;
				smallest = Math.abs(Eval[i]);
			}
		}
		
		//make sure the smallest eigenvalue is zero
		// assertEquals(Math.abs(Eval[index]), 0, 1e-12);
		// System.out.println("Eigenvalue" + " : " + Eval[index]);

		//normalize the eigenvector
		double sum = 0.0;
		double[] statDist = new double[nrOfStates];
		for (int k = 0; k < nrOfStates; k++) {
			sum += Ievc[index * nrOfStates + k];
		}

		for (int k = 0; k < nrOfStates; k++) {
			statDist[k] = Ievc[index * nrOfStates + k] / sum;
		}
		return statDist;
	}

	protected void setupRateMatrix() {
		// Note that in setting up the rate matrix, we always assume minRepeat=0, since
		// the data is already corrected for minRepeat in FiniteIntegerData
		final double rb = rbInput.get().getValue();
		final double ieq = ieqInput.get().getValue() - minRepeat;//translation to provide correct output in the logs
		final double g = gInput.get().getValue();
		final double oneOnA1 = oneOnA1Input.get().getValue();

		double b0 = rb * 1 / Math.sqrt(1 + 1 / (ieq * ieq));
		double b1 = rb * -1 / (Math.sqrt(ieq * ieq + 1));

		double alpha, beta, gamma = 1;
		rowSum = new double[nrOfStates];
		rowSum2 = new double[nrOfStates];
	
		for (int i = 0; i < nrOfStates; i++) {
			rowSum[i] = 0.0;
			rowSum2[i] = 0.0;

			/*
			/* The following is equivalent to:
			/* 1.0 + oneOnA1 * (i - 0)
			*/
			alpha = oneOnA1 + (i - 0);
			beta = 1 / (1 + Math.exp(-(b0 + b1 * (i - 0))));

			for (int j = 0; j < nrOfStates; j++) {
				if (j == i + 1) {
					gamma = (1 - g) / (1 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
					if(Double.isNaN(gamma)) { // if g = 1.0
						gamma = 1 / (double) (nrOfStates - 1 - i);
					}
					rateMatrix[i][j] = alpha * beta * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);

				} else if (j > i + 1) {
					gamma = (1 - g) / (1 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1);
					if(Double.isNaN(gamma)) { // if g = 1.0
						gamma = 1 / (double) (nrOfStates - 1 - i);
					}
					rateMatrix[i][j] = alpha * beta * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				} else if (j == i - 1) {
					gamma = (1 - g) / (1 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
					if(Double.isNaN(gamma)) { // if g = 1.0
						gamma = 1.0 / (double) i;
					}
					rateMatrix[i][j] = alpha * (1 - beta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				} else if (j < i - 1) {
					gamma = (1 - g) / (1 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1);
					if(Double.isNaN(gamma)) { // if g = 1.0
						gamma = 1.0 / (double) i;
					}
					rateMatrix[i][j] = alpha * (1 - beta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				}
			}
			rateMatrix[i][i] = -rowSum[i];
		}
	}

	@Override  // copied from GeneralSubstitutionModel.java
	protected boolean requiresRecalculation() {
		// we only get here if something is dirty
		updateMatrix = true;
		return true;
	}

	@Override // copied from GeneralSubstitutionModel.java
	public void store() {
		storedUpdateMatrix = updateMatrix;
		if( eigenDecomposition != null ) {
			storedEigenDecomposition = eigenDecomposition.copy();
		}
		super.store();
	}

	@Override // copied from GeneralSubstitutionModel.java
	public void restore() {
		updateMatrix = storedUpdateMatrix;
		if( storedEigenDecomposition != null ) {
			EigenDecomposition tmp = storedEigenDecomposition;
			storedEigenDecomposition = eigenDecomposition;
			eigenDecomposition = tmp;
		}
		super.restore();
	}

	@Override // copied from GeneralSubstitutionModel.java
	public EigenDecomposition getEigenDecomposition(Node node) {
		synchronized (this) {
			if (updateMatrix) {
				setupRateMatrix();
				eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
				updateMatrix = false;
			}
		}
		return eigenDecomposition;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof FiniteIntegerData;
	}
}
