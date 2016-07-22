/*
* File SainudiinVanilla.java
*
* Copyright (C) 2016 Arjun Dhawan, RIVM <arjun.dhawan@rivm.nl>
*
* This file is part of BEASTvntr.
*
* BEASTvntr is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BEASTvntr is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with BEASTvntr.  If not, see <http://www.gnu.org/licenses/>.
*/

package beast.evolution.substitutionmodel;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.parameter.IntegerParameter;
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
	"  Genetics 168:383–395", year = 2004, firstAuthorSurname = "sainudiin")
@Citation(value =
	"Chieh-Hsi Wu and  Alexei J. Drummond. (2011) Joint Inference of\n" +
	"  Microsatellite Mutation Models, Population History and Genealogies\n" + 
	"  Using Transdimensional Markov Chain Monte Carlo.\n" + 
	"  Genetics 188:151-164",
	DOI= "10.1534/genetics.103.022665", year = 2011, firstAuthorSurname = "wu")

public class SainudiinVanilla extends Sainudiin {
	final public Input<RealParameter> b0Input = new Input<>("b0", "constant bias parameter of mutational bias beta", Validate.REQUIRED);
	final public Input<RealParameter> b1Input = new Input<>("b1", "linear bias parameter of mutational bias beta", Validate.REQUIRED);
	final public Input<RealParameter> a1Input = new Input<>("a1", "the proportionality of mutation rate to repeat length (i - minimum repeat)", Validate.REQUIRED);

	public SainudiinVanilla() {
		biasMagnitudeInput.setRule(Validate.OPTIONAL);
		focalPointInput.setRule(Validate.OPTIONAL);
		oneOnA1Input.setRule(Validate.OPTIONAL);
	}
	
	@Override
	protected void setParameterBounds() {
		b0Input.get().setBounds(b0Input.get().getLower(), b0Input.get().getUpper());
		b1Input.get().setBounds(b1Input.get().getLower(), b1Input.get().getUpper());
		gInput.get().setBounds(Math.max(0.0, gInput.get().getLower()), Math.min(1.0, gInput.get().getUpper()));
		a1Input.get().setBounds(Math.max(0.0, a1Input.get().getLower()), a1Input.get().getUpper());
	}
	
	protected void setupRateMatrix() {
		// Since the data is already corrected for minRepeat in FiniteIntegerData,
		// we always assume minRepeat is 0 in the substitution model. Except for
		// parameters focalPoint, which are not from FiniteIntegerData.
		final double g = gInput.get().getValue();
		final double a1 = a1Input.get().getValue();

		final double b0 = b0Input.get().getValue();
		final double b1 = b1Input.get().getValue();

		rowSum = new double[nrOfStates];
		rowSum2 = new double[nrOfStates];
	
		for (int i = 0; i < nrOfStates; i++) {
			rowSum[i] = 0.0;
			rowSum2[i] = 0.0;

			// Note that 1.0 + oneOnA1 * (i - 0) and oneOnA1 + (i - 0) are equivalent.
			double alpha = 1.0 + a1 * (i - 0);
			double oneOnbeta = (1.0 + Math.exp(-(b0 + b1 * (i - 0))));

			for (int j = 0; j < nrOfStates; j++) {
				if (j == i + 1) {
					double gamma = (1 - g) * (Math.pow(g, (int) Math.abs(i - j) - 1) / (1 - Math.pow(g, nrOfStates - 1 - i)));
					// If g = 1.0, we assume the limiting case for gamma
					if(Double.isNaN(gamma)) { 
						gamma = 1 / (double) (nrOfStates - 1 - i);
					}
					rateMatrix[i][j] = (alpha / oneOnbeta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				} else if (j > i + 1) {
					double gamma = (1 - g) * (Math.pow(g, (int) Math.abs(i - j) - 1) / (1 - Math.pow(g, nrOfStates - 1 - i)));
					if(Double.isNaN(gamma)) {
						gamma = 1 / (double) (nrOfStates - 1 - i);
					}
					rateMatrix[i][j] = (alpha / oneOnbeta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				} else if (j == i - 1) {
					double gamma = (1 - g) * (Math.pow(g, (int) Math.abs(i - j) - 1) / (1 - Math.pow(g, i - 0)));
					if(Double.isNaN(gamma)) {
						gamma = 1.0 / (double) i;
					}
					rateMatrix[i][j] = (alpha - alpha / oneOnbeta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				} else if (j < i - 1) {
					double gamma = (1 - g) * (Math.pow(g, (int) Math.abs(i - j) - 1) / (1 - Math.pow(g, i - 0))); 
					if(Double.isNaN(gamma)) {
						gamma = 1.0 / (double) i;
					}
					rateMatrix[i][j] = (alpha - alpha / oneOnbeta) * gamma;
					rowSum[i] += rateMatrix[i][j];
					rowSum2[i] += rateMatrix[i][j] * Math.abs(i - j);
				}
			}
			rateMatrix[i][i] = -rowSum[i];
		}
	}
}
