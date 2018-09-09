package implementation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;

public class MainFragmentsForChallenges {
	/*
	 * Parameter (1): path to directory with challenge files
	 * Parameter (2): path to metadata.csv
	 * Parameter (3): path to directory for challenge fragment output files
	 * Parameter (4): path to file where all fragments have to be written
	 */
	public static void main(String[] args) throws Exception {
		String challengeFolder = args[0];
		String metaDataCsv = args[1];
		PrintWriter writer = new PrintWriter(args[3], "UTF-8");
		
		FileReader metaDataReader = new FileReader(metaDataCsv);
		BufferedReader metaDataRowReader = new BufferedReader(metaDataReader);
		
		String row = metaDataRowReader.readLine();
		long time = 0;
		
		while(row != null) {
			String[] columns = row.split("[|]");
			PrintWriter testWriter = new PrintWriter(args[2] + "\\" + columns[1] + ".csv", "UTF-8");
			String challengeFile = columns[1] + ".txt";
			double exactMass = Double.valueOf(columns[2]);
			String smiles = columns[7];

			List<Double> massesList = Lists.newArrayList();
			List<Double> intsList = Lists.newArrayList();
			
			System.out.println("Reading Challenge file: " + challengeFile);
			FileReader challengeReader = new FileReader(challengeFolder + challengeFile);
			BufferedReader challengeRowReader = new BufferedReader(challengeReader);
			String challengeRow = challengeRowReader.readLine();
			while(challengeRow != null) {
				String[] values = challengeRow.split("\\t");
				massesList.add(Double.valueOf(values[0]));
				intsList.add(Double.valueOf(values[1]));
				
				challengeRow = challengeRowReader.readLine();
			}
			
			Double[] masses = new Double[massesList.size()];
			massesList.toArray(masses);
			Double[] ints = new Double[intsList.size()];
			intsList.toArray(ints);

			long timeStart = System.currentTimeMillis();
			List<FragmentSmilesWithMass> fragmentSmiles = FragmentSingleMoleculeIris.getFragment(smiles, masses, ints, exactMass, writer);
			long timeEnd = System.currentTimeMillis();
			time += (timeEnd - timeStart);
			
			for(FragmentSmilesWithMass fragment : fragmentSmiles) {
				String fragmentSmarts = replaceHydrogens(fragment.getFragmentSmiles());
				testWriter.println(fragmentSmarts + ";" + fragment.getMass());
			}
			
			
			challengeReader.close();
			testWriter.close();
			row = metaDataRowReader.readLine();
		}
		writer.println("\nTime elapsed: " + time + "ms");
		
		System.out.println("FINISHED");
		
		writer.close();
		metaDataReader.close();
		metaDataRowReader.close();
	}
	
	private static String replaceHydrogens(String smarts) {
		String regex = "\\[[a-z]H\\d\\]";
		Pattern pattern = Pattern.compile(regex);
		Matcher matcher = pattern.matcher(smarts);
		while(matcher.find()) {
			String match = matcher.group();
			char atomSymbol = match.charAt(1);
			smarts = matcher.replaceFirst(String.valueOf(atomSymbol));
			matcher = pattern.matcher(smarts);
		}
		return smarts;
	}
}
