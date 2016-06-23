/*
* File Sainudiin.java
*
* Copyright (C) 2016 Arjun Dhawan arjun.dhawan@rivm.nl, RIVM
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
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;

@Description("Substitution model of Sainudiin (R. Sainudiin et al., 2004) using Wu's modification (C. Wu and A.J. Drummond, 2011) for VNTR evolution.")
@Citation(value = "Raazesh Sainudiin et al. (2004) Microsatellite Mutation Models. Genetics 168:83–395", DOI= "10.1534/genetics.103.022665", year = 2004, firstAuthorSurname = "sainudiin")

public class Sainudiin extends SubstitutionModel.Base {
  final public Input<RealParameter> rbInput = new Input<>("rb", "force of attraction to i_eq", Validate.REQUIRED);
	final public Input<RealParameter> ieqInput = new Input<>("ieq", "equilibrium state i_eq - i_min", Validate.REQUIRED);
  final public Input<RealParameter> gInput = new Input<>("g", "parameter of the geometric distribution of step sizes (1 - g = probability of a mutation being single step)", Validate.REQUIRED);
  final public Input<RealParameter> a1Input = new Input<>("a1", "proportionality of mutation rate to repeat length i - i_min", Validate.REQUIRED);
  final public Input<RealParameter> nrOfStatesInput = new Input<>("nrOfStates", "number of states", Validate.REQUIRED);

	protected EigenSystem eigenSystem;
  protected EigenDecomposition eigenDecomposition;
  private EigenDecomposition storedEigenDecomposition = null;
  protected double[][] rateMatrix;

  protected boolean updateMatrix = true;
  private boolean storedUpdateMatrix = true;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		updateMatrix = true;


    frequencies = frequenciesInput.get();
    nrOfStates = nrOfStatesInput.get().getValue().intValue();

    rbInput.get().setBounds(Math.max(0.0, rbInput.get().getLower()), rbInput.get().getUpper());
    ieqInput.get().setBounds(Math.max(0.0, ieqInput.get().getLower()), Math.min(ieqInput.get().getUpper(), nrOfStates - 1.0));
    gInput.get().setBounds(Math.max(0.0, gInput.get().getLower()), Math.min(1.0, gInput.get().getUpper()));
    a1Input.get().setBounds(Math.max(0.0, a1Input.get().getLower()), a1Input.get().getUpper());
		
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

	/**
     * access to (copy of) rate matrix *
     */
  protected double[][] getRateMatrix() {
      return rateMatrix.clone();
  }

	protected void setupRateMatrix() {
    // Note that in setting up the rate matrix, we always assume iMin=0, since
    // the data is already corrected for iMin when parsing in 
    // BeautiAlignmentProvider.java
    final double rb = rbInput.get().getValue();
    final double ieq = ieqInput.get().getValue();
    final double g = gInput.get().getValue();
    final double a1 = a1Input.get().getValue();

    double b0;
    double b1;
    if (ieq == 0) { // make sure we don't divide by zero
      b0 = 0;
      b1 = -1 * rb;
    } else {
      b0 = rb * 1 / Math.sqrt(1 + 1 / (ieq * ieq));
      b1 = rb * -1 / (Math.sqrt(ieq * ieq + 1));
    }

		double alpha, beta, gamma = 1, rowSum;
  
    for (int i = 0; i < nrOfStates; i++) {
      rowSum = 0;

      alpha = 1 + a1 * (i - 0);
      beta = 1 / (1 + Math.exp(-(b0 + b1 * (i - 0))));

      for (int j = 0; j < nrOfStates; j++) {
        if (j == i + 1) {
          if(g == 0) {
            gamma = 1.0;
          } else if (g == 1) {
            gamma = 1 / (double) (nrOfStates - 1 - i);
          } else {
            gamma = (1 - g) / (1 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
          }
          rateMatrix[i][j] = alpha * beta * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j > i + 1) {
          if(g == 0) {
            gamma = 0.0;
          } else if (g == 1) {
            gamma = 1 / (double) (nrOfStates - 1 - i);
          } else {
            gamma = (1 - g) / (1 - Math.pow(g, nrOfStates - 1 - i)) * Math.pow(g, (int) Math.abs(i - j) - 1);
          }
          rateMatrix[i][j] = alpha * beta * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j == i - 1) {
          if(g == 0) {
            gamma = 1.0;
          } else if (g == 1) {
            gamma = 1 / (double) i;//
          } else {
            gamma = (1 - g) / (1 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1); 
          }
          rateMatrix[i][j] = alpha * (1 - beta) * gamma;
          rowSum += rateMatrix[i][j];
        } else if (j < i - 1) {
          if(g == 0) {
            gamma = 0.0;
          } else if (g == 1) {
            gamma = 1 / (double) i;
          } else {
            gamma = (1 - g) / (1 - Math.pow(g, i - 0)) * Math.pow(g, (int) Math.abs(i - j) - 1);
          }
          rateMatrix[i][j] = alpha * (1 - beta) * gamma;
          rowSum += rateMatrix[i][j];
        }       
      }
      rateMatrix[i][i] = -rowSum;
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
    if (dataType instanceof IntegerData) {
      dataType.setStateCount(nrOfStates);
      return true;
    } else {
      return false;
    }
		//return dataType instanceof IntegerData;
	}
}
