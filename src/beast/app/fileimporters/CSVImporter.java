/*
* File CSVImporter.java
*
* Copyright (C) 2016 Arjun Dhawan, RIVM <arjun.dhawan@rivm.nl>
* Copyright (C) 2017 Arjun Dhawan
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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

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

  private int nrOfLoci;
  private int nrOfTaxa;

  public enum ParseOption {
    REPEATS_HOMOGEN("Repeats (single partition)"), REPEATS_INHOMOGEN("Repeats (multiple partitions)"),
      NUCLEOTIDES("Nucleotides");
    private final String display;

    ParseOption(String display) {
      this.display = display;
    }
    @Override
    public String toString() {
      return display;
    }
  }

  /**
   * Returns a 2d array of Strings, containing values read from a CSV file.
   */
  private String[][] parseFile(File file, boolean hasHeader) throws IOException {
    BufferedReader fin = new BufferedReader(new FileReader(file));
    List<String[]> taxaData = new ArrayList<>();
    nrOfTaxa = 0;
    int lineNumber = 0;
    while (fin.ready()) {
      lineNumber++;
      String line = fin.readLine();

      if (hasHeader && lineNumber == 1) {
        // Skip header
        continue;
      }
      if (line.trim().length() == 0) {
        // Skip blank lines
        continue;
      }

      String[] taxonData = line.split(",");
      if(nrOfTaxa > 0 && nrOfLoci != taxonData.length - 1) {
        throw new IllegalArgumentException("Wrong number of loci at line " + lineNumber + ". Expected: " + nrOfLoci + ", got: " + (taxonData.length - 1));
      }
      nrOfLoci = taxonData.length - 1;
      taxaData.add(taxonData);
      nrOfTaxa++;
    }
    fin.close();

    if (taxaData.size() == 0)
      throw new IllegalArgumentException("No taxa specified in " + file.getName());

    return taxaData.toArray(new String[][] {});
  }

  /**
   * Returns a list of Alignments
   */
  public List<BEASTInterface> loadFile(File file) {
    // Ask the user for parsing options
    int minRepeat = -1, maxRepeat = -1;
    JPanel myPanel = new JPanel();
    myPanel.add(new JLabel("Parse as:"));
    JComboBox<ParseOption> comboBox = new JComboBox<>(ParseOption.values());
    myPanel.add(comboBox);
    JCheckBox hasHeader = new JCheckBox("Skip header");
    myPanel.add(hasHeader);

    int result = JOptionPane.showConfirmDialog(null, myPanel,
             "Choose parsing method", JOptionPane.OK_CANCEL_OPTION);
    switch (result) {
      case JOptionPane.OK_OPTION:
        break;
      case JOptionPane.OK_CANCEL_OPTION:
      default:
        // User clicked cancel
        return new ArrayList<>();
    }
    ParseOption parseOption = (ParseOption) comboBox.getSelectedItem();
    switch (parseOption) {
      case REPEATS_HOMOGEN:
      case REPEATS_INHOMOGEN:
        JTextField minRepeatField = new JTextField(5);
        JTextField maxRepeatField = new JTextField(5);
        myPanel = new JPanel();
        myPanel.setLayout(new BoxLayout(myPanel, BoxLayout.PAGE_AXIS));
        myPanel.add(new JLabel("<html>Choose the Minimum and Maximum repeat<br>(Maximum repeat - Minimum repeat >= 0 must hold)</html>"));
        myPanel.add(new JLabel("Minimum repeat ( >= 0 ):"));
        myPanel.add(minRepeatField);
        myPanel.add(new JLabel("Maximum repeat:"));
        myPanel.add(maxRepeatField);

        result = JOptionPane.showConfirmDialog(null, myPanel,
                 "Please Enter Minimum and Maximum Repeat", JOptionPane.OK_CANCEL_OPTION);
        while (result == JOptionPane.OK_OPTION) {
          try {
            minRepeat = Integer.parseInt(minRepeatField.getText());
            maxRepeat = Integer.parseInt(maxRepeatField.getText());
            break;
          } catch (NumberFormatException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(null, "Not an integer: " + e.getMessage());
            result = JOptionPane.showConfirmDialog(null, myPanel,
                    "Please Enter Minimum and Maximum Repeat", JOptionPane.OK_CANCEL_OPTION);
          }
        }
        if (result != JOptionPane.OK_OPTION) {
          // User clicked cancel
          return new ArrayList<>();
        }
        break;
      case NUCLEOTIDES:
        break;
    }
    // Parse the CSV file, and put the results into alignment(s)
    String ID = file.getName();
    ID = ID.substring(0, ID.lastIndexOf('.')).replaceAll("\\..*", "");
    List<BEASTInterface> alignmentList = new ArrayList<>();
    try {
      String[][] taxaData = parseFile(file, hasHeader.isSelected());
      Alignment alignment;
      List<Sequence> sequenceList = new ArrayList<>();
      FiniteIntegerData finiteIntegerData = new FiniteIntegerData();

      switch(parseOption){
        case REPEATS_HOMOGEN:
          finiteIntegerData.setInputValue("minRepeat", minRepeat);
          finiteIntegerData.setInputValue("maxRepeat", maxRepeat);
          finiteIntegerData.initAndValidate();
          for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
            String taxonName = taxaData[taxonIndex][0];
            String dataString = String.join(",", Arrays.copyOfRange(taxaData[taxonIndex], 1, nrOfLoci + 1));
            Sequence sequence = new Sequence();
            sequence.totalCountInput.setValue(maxRepeat - minRepeat + 1, sequence);
            sequence.taxonInput.setValue(taxonName, sequence);
            sequence.dataInput.setValue(dataString, sequence);
            sequence.setID("seq_" + taxonName);
            sequenceList.add(sequence);
          }
          alignment = new Alignment();
          alignment.sequenceInput.setValue(sequenceList, alignment);
          alignment.setID(ID);
          alignment.userDataTypeInput.setValue(finiteIntegerData, alignment);
          alignment.initAndValidate();
          alignmentList.add(alignment);
          break;
        case REPEATS_INHOMOGEN:
          finiteIntegerData.setInputValue("minRepeat", minRepeat);
          finiteIntegerData.setInputValue("maxRepeat", maxRepeat);
          finiteIntegerData.initAndValidate();
          for (int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
            sequenceList = new ArrayList<>();
            for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
              String taxonName = taxaData[taxonIndex][0];
              Sequence sequence = new Sequence();
              sequence.totalCountInput.setValue(maxRepeat - minRepeat + 1, sequence);
              sequence.taxonInput.setValue(taxonName, sequence);
              sequence.dataInput.setValue(taxaData[taxonIndex][locusIndex + 1], sequence);
              sequence.setID("seq_" + taxonName + "_VNTR" + String.format("%02d", locusIndex + 1));
              sequenceList.add(sequence);
            }
            alignment = new Alignment();
            alignment.sequenceInput.setValue(sequenceList, alignment);
            alignment.userDataTypeInput.setValue(finiteIntegerData, alignment);
            alignment.setID(ID + String.format("%02d", locusIndex + 1));
            alignment.initAndValidate();
            alignmentList.add(alignment);
          }
          break;
        case NUCLEOTIDES:
          for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
            String taxonName = taxaData[taxonIndex][0];
            StringBuilder sb = new StringBuilder();
            for (int locusIndex = 0; locusIndex < nrOfLoci; locusIndex++) {
              sb.append(String.valueOf(taxaData[taxonIndex][locusIndex + 1]));
            }
            Sequence sequence = new Sequence();
            sequence.totalCountInput.setValue(4, sequence);
            sequence.taxonInput.setValue(taxonName, sequence);
            sequence.dataInput.setValue(sb.toString(), sequence);
            sequence.setID("seq_" + taxonName);
            sequenceList.add(sequence);
          }
          alignment = new Alignment();
          alignment.sequenceInput.setValue(sequenceList, alignment);
          alignment.dataTypeInput.setValue("nucleotide", alignment);
          alignment.setID(ID);
          alignment.initAndValidate();
          alignmentList.add(alignment);
          break;
      }
      return alignmentList;
    } catch (Exception e) {
      e.printStackTrace();
      JOptionPane.showMessageDialog(null, "Loading of " + file.getName() + " failed: " + e.getMessage());
      return alignmentList;
    }
  }
}
