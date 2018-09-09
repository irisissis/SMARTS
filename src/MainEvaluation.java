package implementation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import edu.emory.mathcs.backport.java.util.Collections;

public class MainEvaluation {
	/*
	 * Parameter (1): path to base directory to evaluate
	 * Parameter (2): suffix for identification of evaluation result files (optional)
	 */
	public static void main(String... args) throws IOException {
		String directory = args[0];
		String suffix = args.length == 2 ? "-" + args[1] : "";
		Map<String, String> solutions = getSolutions(directory + File.separator + "metadata.csv");
		
		PrintWriter writer = new PrintWriter(directory + File.separator + "evaluation" + File.separator + "evaluation" + suffix + ".csv", "UTF-8");
		
		for(String challenge : solutions.keySet()) {
			String finderFile = directory + File.separator + "results" + File.separator + challenge + "-FragmentFinder" + suffix + ".csv";
			String smartsFile = directory + File.separator + "results" + File.separator + challenge + "-FragmentFinderSmarts" + suffix + ".csv";
			String withoutAromaticityFile = directory + File.separator + "results" + File.separator + challenge + "-FragmentFinderWithoutAromaticity" + suffix + ".csv";
			
			List<String> finderResults = countResultsInFile(finderFile);
			List<String> smartsResults = countResultsInFile(smartsFile);
			List<String> withoutAromaticityResults = countResultsInFile(withoutAromaticityFile);
			
			writer.print(challenge);
			writer.println(";FragmentFinder;" + (finderResults.contains(solutions.get(challenge)) ? "yes" : "no") + ";" + finderResults.size());
			writer.println(";FragmentFinderSmarts;" + (smartsResults.contains(solutions.get(challenge)) ? "yes" : "no") + ";" + smartsResults.size());
			writer.println(";FragmentFinderWithoutAromaticity;" + (withoutAromaticityResults.contains(solutions.get(challenge)) ? "yes" : "no") + ";" + withoutAromaticityResults.size());
			
			System.out.println(challenge + ": " + solutions.get(challenge));
			System.out.println(finderResults);
		}
		writer.close();
		System.out.println("FINISHED");
	}
	
	private static Map<String, String> getSolutions(String metaDataCsv) throws IOException {
		Map<String, String> solutions = new HashMap<String, String>();
		
		FileReader metaDataReader = new FileReader(metaDataCsv);
		BufferedReader metaDataRowReader = new BufferedReader(metaDataReader);
		
		String row = metaDataRowReader.readLine();
		while(row != null) {
			String[] columns = row.split("[|]");
			String challenge = columns[1];
			String smiles = columns[6];
			solutions.put(challenge, smiles);
			row = metaDataRowReader.readLine();
		}
		
		metaDataReader.close();
		metaDataRowReader.close();
		
		return solutions;
	}
	
	private static List<String> countResultsInFile(String filePath) throws IOException {
		List<String> resultSmiles = Lists.newArrayList();
		FileReader metaDataReader = new FileReader(filePath);
		BufferedReader metaDataRowReader = new BufferedReader(metaDataReader);
		String row = metaDataRowReader.readLine();
		while(row != null) {
			if(row.equals("")) {
				row = metaDataRowReader.readLine();
				continue;
			}
			if(!row.contains(";")) {
				row = metaDataRowReader.readLine();
				resultSmiles.add(row);
			}
			row = metaDataRowReader.readLine();
		}
		
		metaDataReader.close();
		metaDataRowReader.close();
		
		return resultSmiles;
	}
}
