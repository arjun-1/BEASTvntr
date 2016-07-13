/*
* File SainudiinFreqsComputed.java
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
	"  Wu's modification (C. Wu and A.J. Drummond, 2011) for VNTR evolution,\n" +
	"  where the frequencies of the root node are computed from the model.")
@Citation(value = 
	"Raazesh Sainudiin et al. (2004) Microsatellite Mutation Models.\n" +
	"  Genetics 168:383–395\n\n" + 
	"Chieh-Hsi Wu and  Alexei J. Drummond. (2011) Joint Inference of\n" +
	"  Microsatellite Mutation Models, Population History and Genealogies\n" + 
	"  Using Transdimensional Markov Chain Monte Carlo.\n" + 
	"  Genetics 188:151-164",
	DOI= "10.1534/genetics.103.022665", year = 2004, firstAuthorSurname = "sainudiin")

public class SainudiinFreqsComputed extends Sainudiin {
	public SainudiinFreqsComputed() {
		// this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
		frequenciesInput.setRule(Validate.OPTIONAL);
	}
	double[] frequencies;

	@Override
	public void initAndValidate() {
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

		rbInput.get().setBounds(Math.max(0.0, rbInput.get().getLower()), rbInput.get().getUpper());
		ieqInput.get().setBounds(ieqInput.get().getLower(), ieqInput.get().getUpper());
		gInput.get().setBounds(Math.max(0.0, gInput.get().getLower()), Math.min(1.0, gInput.get().getUpper()));
		oneOnA1Input.get().setBounds(Math.max(0.0, oneOnA1Input.get().getLower()), oneOnA1Input.get().getUpper());
		
		eigenSystem = new DefaultEigenSystem(nrOfStates);
		rateMatrix = new double[nrOfStates][nrOfStates];

		if (frequenciesInput.get() != null) {
				throw new RuntimeException("Frequencies must not be specified in SainudiinStepWise. The Frquencies are calculated from the other parameters.");
		}
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
		frequencies = stationaryDistribution;
		double normalization = 0.0;

		for (i = 0; i < nrOfStates; i++) {
			//normalization += statFromEigen[i] * rowSum[i];
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

	@Override
	public double[] getFrequencies() {
		if(frequencies == null) { // This is the case only during initialization
			frequencies = new double [nrOfStates];
			for(int i = 0; i < nrOfStates; i++) {
				frequencies[i] = 1.0 / (double) nrOfStates;
			}
		}
		return frequencies;
	}
}