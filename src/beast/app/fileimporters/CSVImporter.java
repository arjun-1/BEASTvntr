package beast.app.fileimporters;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.List;
import java.util.Arrays;

import javax.swing.JOptionPane;

import beast.core.BEASTInterface;
import beast.app.beauti.AlignmentImporter;

import beast.evolution.datatype.FiniteIntegerData;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;

public class CSVImporter implements AlignmentImporter {

	public String [] getFileExtensions() {
		return new String[]{"csv"};
	}

	public List<BEASTInterface> loadFile(File file) {
  	try {
			//Set the number of states to initialize the sequence with.
			int min_i = 0;
			int max_i = 14;
			FiniteIntegerData type = new FiniteIntegerData();
      type.setInputValue("min_i", min_i);
      type.setInputValue("max_i", max_i);
    	type.initAndValidate();
      
			// start reading the csv file
			BufferedReader fin = new BufferedReader(new FileReader(file));
			
			// the possible datatypes we can parse from the csv
			String[] choicesDataType = {"finiteinteger", "nucleotide"};
			String datatype = (String) JOptionPane.showInputDialog(
				null,
				"Choose DataType (default: integer)",
				"Choose DataType",
				JOptionPane.PLAIN_MESSAGE,
				null,
				choicesDataType,
				choicesDataType[0]
				);

			// if the user clicks 'cancel', we go for integer
			if (datatype == null) {
				datatype = "finiteinteger";
				// alert the user of our assumption
				JOptionPane.showMessageDialog(null, "DataType is assumed integer");
			}

      Alignment m_alignment = new Alignment();

      while (fin.ready()) {
      	String line = fin.readLine();
      	if ( line.trim().length() == 0 ) {
						continue;  // Skip blank lines
				}
      	String[] row = line.split(",");
      	String currentTaxon = row[0];

				if (currentTaxon == null || currentTaxon.trim().length() == 0) {
					// check if not null and not only white space
					fin.close();
					throw new RuntimeException("Expected taxon defined on first line");
				}
      	StringBuilder sb = new StringBuilder();
      	int rowLength = row.length;
      	
      	for (int i = 1; i < row.length; i++) {
      		if (datatype == "finiteinteger") {
      				int parsedAllel = Integer.parseInt(row[i]);
      				if (parsedAllel - min_i < 0) {
								throw new RuntimeException("Encountered i < min_i");
							} else if (parsedAllel > max_i) {
								throw new RuntimeException("Encountered i > max_i");
							}
      				sb.append(parsedAllel);
      				sb.append(",");
      		} else if (datatype == "nucleotide") {
      			sb.append(row[i]);
      		}
      	}
      	String sequenceData = sb.toString();
      	Sequence sequence = new Sequence();
      	
      	if (datatype == "finiteinteger") {
      		sequence.init(max_i - min_i + 1, currentTaxon, sequenceData);
      	} else if (datatype == "nucleotide") {
      		sequence.init(4, currentTaxon, sequenceData);
      	}
      	sequence.setID("seq_" + currentTaxon);
      	m_alignment.sequenceInput.setValue(sequence, m_alignment);
      }
      fin.close();
      
      String ID = file.getName();
      ID = ID.substring(0, ID.lastIndexOf('.')).replaceAll("\\..*", "");
      m_alignment.setID(ID);

      if (datatype == "finiteinteger") {
      	m_alignment.setInputValue("userDataType", type);
      } else if (datatype == "nucleotide") {
      	m_alignment.dataTypeInput.setValue(datatype, m_alignment);
      }

      m_alignment.initAndValidate();
      return Arrays.asList(m_alignment);
  		
  	} catch (Exception e) {
			e.printStackTrace();
			JOptionPane.showMessageDialog(null, "Loading of " + file.getName() + " failed: " + e.getMessage());
			Alignment dummy = new Alignment();
			return Arrays.asList(dummy);
  	}
	}
}
