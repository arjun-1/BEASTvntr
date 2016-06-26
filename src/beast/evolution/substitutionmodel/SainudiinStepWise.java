/*
* File SainudiinStepWise.java
*
* Copyright (C) 2016 Arjun Dhawan arjun.dhawan@rivm.nl, RIVM
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

@Description("Substitution model of Sainudiin (R. Sainudiin et al., 2004) using Wu's modification (C. Wu and A.J. Drummond, 2011) for VNTR evolution.")
@Citation(value = "Raazesh Sainudiin et al. (2004) Microsatellite Mutation Models. Genetics 168:83–395", DOI= "10.1534/genetics.103.022665", year = 2004, firstAuthorSurname = "sainudiin")
public class SainudiinStepWise extends Sainudiin {
  public SainudiinStepWise() {
    // this is added to avoid a parsing error inherited from superclass because frequencies and g are not provided.
    frequenciesInput.setRule(Validate.OPTIONAL);
    gInput.setRule(Validate.OPTIONAL);
  }
  double[] frequencies;

  @Override
	public void initAndValidate() {
		updateMatrix = true;

    for (Object beastObjecti : getOutputs()) {
      if (beastObjecti instanceof SiteModel) {
        SiteModel sitemodel = (SiteModel) beastObjecti;
        for (Object beastObjectj : sitemodel.getOutputs()) {
          if (beastObjectj instanceof ThreadedTreeLikelihood) {
            ThreadedTreeLikelihood likelihood = (ThreadedTreeLikelihood) beastObjectj;
            nrOfStates = likelihood.dataInput.get().getMaxStateCount();
            break;
          }
        }
        break;
      }
    }

    rbInput.get().setBounds(Math.max(0.0, rbInput.get().getLower()), rbInput.get().getUpper());
    ieqInput.get().setBounds(Math.max(0.0, ieqInput.get().getLower()), Math.min(ieqInput.get().getUpper(), (float) nrOfStates - 1.0));
    a1Input.get().setBounds(Math.max(0.0, a1Input.get().getLower()), a1Input.get().getUpper());
				
		eigenSystem = new DefaultEigenSystem(nrOfStates);
		rateMatrix = new double[nrOfStates][nrOfStates];

    if (frequenciesInput.get() != null) {
        throw new RuntimeException("Frequencies must not be specified in SainudiinStepWise. They are calculated internally.");
    }
    frequencies = new double[nrOfStates];
	}

  @Override
	protected void setupRateMatrix() {
    // Note that in setting up the rate matrix, we always assume iMin=0, since
    // the data is already corrected for iMin when parsing in 
    // BeautiAlignmentProvider.java
    // This rate matrix is almost the same as that of Sainudiin.java. 
    // Here however, we calculate the freqs from the stationary distribution.
    final double rb = rbInput.get().getValue();
    final double ieq = ieqInput.get().getValue();
    final double g = 0.0; //stepwise    
    final double a1 = a1Input.get().getValue();

    double b0;
    double b1;
    if (ieq == 0) {
      b0 = 0;
      b1 = -1 * rb;
    } else {
      b0 = rb * 1.0 / Math.sqrt(1.0 + 1.0 / (ieq * ieq));
      b1 = rb * -1.0 / (Math.sqrt(ieq * ieq + 1.0));
    }

		double alpha = 1.0, beta = 1.0, gamma = 1.0, rowSum;
    // extra variables needed for calculation of frequencies:
    double alphaold, betaold, sumfreq = 0.0, birthj, deathjplus1;
  
    for (int i = 0; i < nrOfStates; i++) {
      rowSum = 0.0;

      alphaold = alpha;
      betaold = beta;
      alpha = 1.0 + a1 * (i - 0);
      beta = 1.0 / (1.0 + Math.exp(-(b0 + b1 * (i - 0))));

      for (int j = 0; j < nrOfStates; j++) {
        if (j == i + 1) {
          if(g == 0.0) {
            gamma = 1.0;
          } else if (g == 1) {
            gamma = 1.0 / (double) (nrOfStates - 1 - i);
          } else {
            gamma = (1.0 - g) / (1.0 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
          }
          rateMatrix[i][j] = alpha * beta * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j > i + 1) {
          if(g == 0) {
            gamma = 0.0;
          } else if (g == 1.0) {
            gamma = 1 / (double) (nrOfStates - 1 - i);
          } else {
            gamma = (1.0 - g) / (1.0 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1);
          }
          rateMatrix[i][j] = alpha * beta * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j == i - 1) {
          if(g == 0.0) {
            gamma = 1.0;
          } else if (g == 1.0) {
            gamma = 1.0 / (double) i;//
          } else {
            gamma = (1.0 - g) / (1.0 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
          }
          rateMatrix[i][j] = alpha * (1 - beta) * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j < i - 1) {
          if(g == 0.0) {
            gamma = 0.0;
          } else if (g == 1.0) {
            gamma = 1.0 / (double) i;
          } else {
            gamma = (1.0 - g) / (1.0 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1);
          }
          rateMatrix[i][j] = alpha * (1.0 - beta) * gamma;
          rowSum += rateMatrix[i][j];
        }       
      }
      rateMatrix[i][i] = -rowSum;

      // Calc frequencies from stationary distribution, which is a special case
      // of the birth death chain
      birthj      = alphaold * betaold;
      deathjplus1 = alpha * (1.0 - beta);

      if(i == 0) 
        frequencies[i] = 1.0;
      else
        frequencies[i] = frequencies[i - 1] * birthj / deathjplus1;
      sumfreq += frequencies[i];
    }
    for(int i = 0; i < nrOfStates; i++) {
      frequencies[i] /= sumfreq;
    }
	}

  @Override
  public double[] getFrequencies() {
    return frequencies;
  }

  @Override
  public boolean canHandleDataType(DataType dataType) {
    return dataType instanceof FiniteIntegerData;
  }
}
