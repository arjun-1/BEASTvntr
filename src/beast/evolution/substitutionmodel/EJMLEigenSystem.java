/*
* File EJMLEigenSystem.java
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

import java.util.Arrays;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.ops.CommonOps;
import org.ejml.ops.EigenOps;

public class EJMLEigenSystem implements EigenSystem {
  private final int stateCount;

  public EJMLEigenSystem(int stateCount) {
    this.stateCount = stateCount;
  }

  @Override
  public EigenDecomposition decomposeMatrix(double[][] qMatrix) {
    DenseMatrix64F rateMatrix = new DenseMatrix64F(qMatrix);
    org.ejml.interfaces.decomposition.EigenDecomposition<DenseMatrix64F> eig = DecompositionFactory.eig(stateCount, true);
    eig.decompose(rateMatrix);
    DenseMatrix64F eigenVectors = EigenOps.createMatrixV(eig);
    //DenseMatrix64F eigenValues = EigenOps.createMatrixD(eig);
    DenseMatrix64F inverseEigenVectors = new DenseMatrix64F(stateCount, stateCount);
    CommonOps.invert(eigenVectors, inverseEigenVectors);

    double[] flatEvec = eigenVectors.getData();
    double[] flatIevc = inverseEigenVectors.getData();

    double[] Eval = new double[stateCount];
    for(int i=0; i<stateCount; i++) {
      Eval[i] = eig.getEigenvalue(i).real;
    }
    return new EigenDecomposition(flatEvec, flatIevc, Eval);
  }
}

