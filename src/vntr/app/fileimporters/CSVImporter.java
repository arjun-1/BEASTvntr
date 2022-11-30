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

package vntr.app.fileimporters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.ArrayList;

//import javax.swing.*;

import beast.base.core.BEASTInterface;
import beastfx.app.beauti.ThemeProvider;
import beastfx.app.inputeditor.AlignmentImporter;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import vntr.evolution.datatype.FiniteIntegerData;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;

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
    	  fin.close();
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
    HBox myPanel0 = FXUtils.newHBox();
    Label label = new Label("Parse as:");
    myPanel0.getChildren().add(label);
    label.setPadding(new Insets(5));
    ComboBox<ParseOption> comboBox = new ComboBox<>();
    for (ParseOption po : ParseOption.values()) {
    	comboBox.getItems().add(po);
    }
    comboBox.setMinWidth(200);
    comboBox.getSelectionModel().select(0);
    myPanel0.getChildren().add(comboBox);
    CheckBox hasHeader = new CheckBox("Skip header");
    hasHeader.setMinWidth(200);
    hasHeader.setPadding(new Insets(5));
    myPanel0.getChildren().add(hasHeader);

	Dialog dlg = new Dialog();
	dlg.getDialogPane().setContent(myPanel0);
	dlg.getDialogPane().getButtonTypes().addAll(Alert.OK_CANCEL_OPTION);
	dlg.setResizable(true);
	Scene node = dlg.getDialogPane().getScene();
	dlg.setX(node.getX() + node.getWidth()/2);
	dlg.setY(node.getY() + node.getHeight()/2);
	ThemeProvider.loadStyleSheet(dlg.getDialogPane().getScene());
	Optional result = dlg.showAndWait();

	if (result.toString().toLowerCase().contains("cancel")) { 
        // User clicked cancel
        return new ArrayList<>();
    }
	
    ParseOption parseOption = (ParseOption) comboBox.getSelectionModel().getSelectedItem();
    switch (parseOption) {
      case REPEATS_HOMOGEN:
      case REPEATS_INHOMOGEN:
        TextField minRepeatField = new TextField();
        minRepeatField.setPrefColumnCount(5);
        TextField maxRepeatField = new TextField();
        maxRepeatField.setPrefColumnCount(5);
        VBox myPanel = FXUtils.newVBox();
        myPanel.getChildren().add(new Label("Choose the Minimum and Maximum repeat\n(Maximum repeat - Minimum repeat >= 0 must hold)"));
        myPanel.getChildren().add(new Label("Minimum repeat ( >= 0 ):"));
        myPanel.getChildren().add(minRepeatField);
        myPanel.getChildren().add(new Label("Maximum repeat:"));
        myPanel.getChildren().add(maxRepeatField);


    	Dialog dlg2 = new Dialog();
    	dlg2.getDialogPane().setContent(myPanel);
    	dlg2.getDialogPane().getButtonTypes().addAll(Alert.OK_CANCEL_OPTION);
    	dlg2.setHeaderText("Please Enter Minimum and Maximum Repeat");
    	dlg2.setResizable(true);
    	Scene node2 = dlg2.getDialogPane().getScene();
    	dlg2.setX(node2.getX() + node.getWidth()/2);
    	dlg2.setY(node2.getY() + node.getHeight()/2);
    	ThemeProvider.loadStyleSheet(dlg.getDialogPane().getScene());
    	Optional result2 = dlg2.showAndWait();

        while (result2.toString().toLowerCase().contains("ok")) {
          try {
            minRepeat = Integer.parseInt(minRepeatField.getText());
            maxRepeat = Integer.parseInt(maxRepeatField.getText());
            break;
          } catch (NumberFormatException e) {
            e.printStackTrace();
            Alert.showMessageDialog(null, "Not an integer: " + e.getMessage());
            result2 =  dlg.showAndWait();
          }
        }
        if (!result2.toString().toLowerCase().contains("ok")) {
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
      Alert.showMessageDialog(null, "Loading of " + file.getName() + " failed: " + e.getMessage());
      return alignmentList;
    }
  }
}
