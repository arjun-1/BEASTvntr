/*
* File SainudiinComputedFrequencies.java
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
import beast.core.Input.Validate;

@Description(
  "Substitution model of Sainudiin (R. Sainudiin et al., 2004) using\n" +
  "  Wu's modification (C. Wu and A.J. Drummond, 2011) for VNTR evolution,\n" +
  "  where the frequencies of the root node are computed from the model.")
@Citation(value = 
  "Raazesh Sainudiin et al. (2004) Microsatellite Mutation Models.\n" +
  "  Genetics 168:383-395", year = 2004, firstAuthorSurname = "sainudiin")
@Citation(value =
  "Chieh-Hsi Wu and  Alexei J. Drummond. (2011) Joint Inference of\n" +
  "  Microsatellite Mutation Models, Population History and Genealogies\n" + 
  "  Using Transdimensional Markov Chain Monte Carlo.\n" + 
  "  Genetics 188:151-164",
  DOI= "10.1534/genetics.103.022665", year = 2011, firstAuthorSurname = "wu")

public class SainudiinComputedFrequencies extends Sainudiin {
  public SainudiinComputedFrequencies() {
    // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
    frequenciesInput.setRule(Validate.OPTIONAL);
  }

  @Override
  public void initAndValidate() {
    updateMatrix = true;
    setStateBoundsFromAlignment();
    setParameterBounds();
    
    //eigenSystem = new DefaultEigenSystem(nrOfStates);
    eigenSystem = new EJMLEigenSystem(nrOfStates);
    rateMatrix = new double[nrOfStates][nrOfStates];

    if (frequenciesInput.get() != null) {
      throw new RuntimeException("Frequencies must not be specified in SainudiinComputedFrequencies. The frequencies are calculated from the other parameters.");
    }
  }

  @Override
  public double[] getFrequencies() {
    // During initialization, the stationary distribution is not yet set.
    // If this is the case, set it, and assume all frequencies equal.
    if(stationaryDistribution == null) {
      stationaryDistribution = new double[nrOfStates];
      for(int i = 0; i < nrOfStates; i++) {
        stationaryDistribution[i] = 1.0 / (double) nrOfStates;
      }
    }
    return stationaryDistribution;
  }
}