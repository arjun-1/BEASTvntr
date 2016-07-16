/*
* File CSVImporter.java
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

package beast.app.fileimporters;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.IOException;

import javax.swing.JOptionPane;
import javax.swing.*;



import beast.core.BEASTInterface;
import beast.app.beauti.AlignmentImporter;

import beast.evolution.datatype.FiniteIntegerData;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;


public class CSVImporter implements AlignmentImporter {
	
	public String [] getFileExtensions() {
		return new String[]{"csv"};
	}

	enum ParseOption {
		REPEATS, REPEATS_INHOMOGEN, NUCLEOTIDES
	}

	private int[][] sequenceData;
	private String[][] sequenceDataNucleotides;
	private String[] taxaNames;
	private int nrOfTaxa;
	private int nrOfLoci;

	private void parseFileAsRepeats(File file) throws IOException {
		BufferedReader fin = new BufferedReader(new FileReader(file));
    List<String> sequenceList = new ArrayList<String>();

    while (fin.ready()) {
    	String line = fin.readLine();
    	if ( line.trim().length() == 0 ) {
					continue;  // Skip blank lines
			}
			sequenceList.add(line);
		}
		fin.close();

		nrOfTaxa = sequenceList.size();
		nrOfLoci = sequenceList.get(0).split(",").length - 1;
		sequenceData = new int[nrOfTaxa][nrOfLoci];
		taxaNames = new String[nrOfTaxa];

		int taxonIndex = 0;
		for(String sequenceString : sequenceList) {
			String[] splitSequenceString = sequenceString.split(",");
    	taxaNames[taxonIndex] = splitSequenceString[0];

    	if (taxaNames[taxonIndex] == null || taxaNames[taxonIndex].trim().length() == 0) {
				// check if not null and not only white space
				throw new IllegalArgumentException("Expected taxon defined on first line");
			}
			for (int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
				int parsedAllel = Integer.parseInt(splitSequenceString[locusIndex + 1]);
				sequenceData[taxonIndex][locusIndex] = parsedAllel;
			}
			taxonIndex += 1;
		}
	}

	private void parseFileAsNucleotides(File file) throws IOException {
		BufferedReader fin = new BufferedReader(new FileReader(file));
    List<String> sequenceList = new ArrayList<String>();

    while (fin.ready()) {
    	String line = fin.readLine();
    	if ( line.trim().length() == 0 ) {
					continue;  // Skip blank lines
			}
			sequenceList.add(line);
		}
		fin.close();

		nrOfTaxa = sequenceList.size();
		nrOfLoci = sequenceList.get(0).split(",").length - 1;
		sequenceDataNucleotides = new String[nrOfTaxa][nrOfLoci];
		taxaNames = new String[nrOfTaxa];

		int taxonIndex = 0;
		for(String sequenceString : sequenceList) {
			String[] splitSequenceString = sequenceString.split(",");
    	taxaNames[taxonIndex] = splitSequenceString[0];

    	if (taxaNames[taxonIndex] == null || taxaNames[taxonIndex].trim().length() == 0) {
				// check if not null and not only white space
				throw new IllegalArgumentException("Expected taxon defined on first line");
			}
			for (int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
				String parsedAllel = splitSequenceString[locusIndex + 1];
				sequenceDataNucleotides[taxonIndex][locusIndex] = parsedAllel;
			}
			taxonIndex += 1;
		}
	}

	public List<BEASTInterface> loadFile(File file) {    
		Alignment alignment;
		FiniteIntegerData finiteIntegerData = new FiniteIntegerData();;
		ArrayList<BEASTInterface> alignmentList = new ArrayList<BEASTInterface>();;

		int minRepeat = 1, maxRepeat = 15;
		ParseOption choice = ParseOption.REPEATS;

		String[] parseOptions = {"Repeats", "Repeats (Inhomogen)", "Nucleotides"};
		String	selectedOption = (String) JOptionPane.showInputDialog(
				null,
				"Choose Parsing Option (default: Repeats)",
				"Choose Parsing Option",
				JOptionPane.PLAIN_MESSAGE,
				null,
				parseOptions,
				parseOptions[0]
			);
			switch (selectedOption) {
				case "Repeats":
					choice = ParseOption.REPEATS;
					break;
				case "Repeats (Inhomogen)":
					choice = ParseOption.REPEATS_INHOMOGEN;
					break;
				case "Nucleotides":
					choice = ParseOption.NUCLEOTIDES;
					break;
				default:
					throw new IllegalArgumentException("No parsing option specified");
			}
			
			switch (choice) {
				case REPEATS:
				case REPEATS_INHOMOGEN:
		      JTextField minRepeatField = new JTextField(5);
		      JTextField maxRepeatField = new JTextField(5);

		      JPanel myPanel = new JPanel();
		      myPanel.setLayout(new BoxLayout(myPanel, BoxLayout.PAGE_AXIS));

		      myPanel.add(new JLabel("<html>Choose the Minimum and Maximum repeat<br>(Maximum repeat - Minimum repeat >= 0 must hold)</html>"));
		      myPanel.add(new JLabel("Minimum repeat ( >= 0 ):"));
		      myPanel.add(minRepeatField);
		      myPanel.add(new JLabel("Maximum repeat:"));
		      myPanel.add(maxRepeatField);

		      int result = JOptionPane.showConfirmDialog(null, myPanel, 
		               "Please Enter Minimum and Maximum Repeat", JOptionPane.OK_CANCEL_OPTION);
		      if (result == JOptionPane.OK_OPTION) {
		      	minRepeat = Integer.parseInt(minRepeatField.getText());
		        maxRepeat = Integer.parseInt(maxRepeatField.getText());
					} else {
						throw new IllegalArgumentException("Minimum and maximum repeat were not specified");
					}
					break;
				case NUCLEOTIDES:
					break;
		}


		String ID = file.getName();
		ID = ID.substring(0, ID.lastIndexOf('.')).replaceAll("\\..*", "");

		try {
			// Put the sequences into alignments
			switch(choice) {
				case REPEATS:
				case REPEATS_INHOMOGEN:
					parseFileAsRepeats(file);
					finiteIntegerData.setInputValue("minRepeat", minRepeat);
					finiteIntegerData.setInputValue("maxRepeat", maxRepeat);
					finiteIntegerData.initAndValidate();
				break;
				case NUCLEOTIDES:
					parseFileAsNucleotides(file);
			}
			switch(choice){
				case REPEATS:
					// case of homogeneous
					alignment = new Alignment();
					for(int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
						StringBuilder sb = new StringBuilder();
						for(int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
							sb.append(sequenceData[taxonIndex][locusIndex]);
			    		sb.append(",");
						}
						Sequence sequence = new Sequence();
						sequence.init(maxRepeat - minRepeat + 1, taxaNames[taxonIndex], sb.toString());
						sequence.setID("seq_" + taxaNames[taxonIndex]);
						alignment.sequenceInput.setValue(sequence, alignment);
					}
					alignment.setID(ID);
					alignment.setInputValue("userDataType", finiteIntegerData);
					alignment.initAndValidate();
					alignmentList.add(alignment);
					break;
				case REPEATS_INHOMOGEN:
					for(int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
						alignment = new Alignment();
						for(int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
							Sequence sequence = new Sequence();
							sequence.init(maxRepeat - minRepeat + 1, taxaNames[taxonIndex], Integer.toString(sequenceData[taxonIndex][locusIndex]) + ",");
							sequence.setID("seq_" + taxaNames[taxonIndex]);
							alignment.sequenceInput.setValue(sequence, alignment);
						}
						alignment.setID("VNTR" + String.format("%02d", locusIndex + 1));
						alignment.setInputValue("userDataType", finiteIntegerData);
						alignment.initAndValidate();
						alignmentList.add(alignment);
					}
					break;
				case NUCLEOTIDES:
					alignment = new Alignment();
					for(int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
						StringBuilder sb = new StringBuilder();
						for(int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
							sb.append(sequenceDataNucleotides[taxonIndex][locusIndex]);
			    		sb.append(",");
						}
						Sequence sequence = new Sequence();
						sequence.init(4, taxaNames[taxonIndex], sb.toString());
						sequence.setID("seq_" + taxaNames[taxonIndex]);
						alignment.sequenceInput.setValue(sequence, alignment);
					}
					alignment.setID(ID);
					alignment.dataTypeInput.setValue("nucleotide", alignment);
					alignment.initAndValidate();
					alignmentList.add(alignment);
					break;
			}
			return alignmentList;
		} catch (Exception e) {
			e.printStackTrace();
			JOptionPane.showMessageDialog(null, "Loading of " + file.getName() + " failed: " + e.getMessage());
			return null;
		}
	}
}
