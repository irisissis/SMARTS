package implementation;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;

public class FragmentFinderMZDatabase {
	/*
	 * Parameter (1): path to base directory
	 */
	public static void run(String directory) throws Exception {
		List<String> challenges = getChallenges(directory + File.separator + "metadata.csv");
		
		ListMultimap<String, String> fragmentSmarts = ArrayListMultimap.create();
		for(String challenge : challenges) {
			String fragmentsFile = directory + File.separator + "fragments" + File.separator + challenge + ".csv";
			fragmentSmarts(fragmentSmarts, fragmentsFile);
		}
		
		long timeFragmentFinder = 0;
		long timeFragmentFinderSmarts = 0;
		long timeFragmentFinderWithoutAromaticity = 0;
		for(String challenge : challenges) {
			System.out.println("Challenge: " + challenge);
			String moleculeFile = directory + File.separator + "molecules" + File.separator + challenge + ".csv";
			Map<String, String> moleculeSmiles = moleculeSmiles(moleculeFile);
			String mzFile = directory + File.separator + "fragments" + File.separator + challenge + ".csv";
			List<String> mzValues = mzValues(mzFile);
			
			Map<String, Collection<String>> fragmentSmartsMap = fragmentSmarts.asMap();
			List<Collection<String>> fragments = Lists.newArrayList();
			for(String mzValue : mzValues) {
				fragments.add(fragmentSmartsMap.get(mzValue));
			}
			
			long timeStart = System.currentTimeMillis();
			fragmentFinder(moleculeSmiles, fragments, directory, challenge);
			long timeEnd = System.currentTimeMillis();
			timeFragmentFinder += (timeEnd - timeStart);
			System.out.println("FragmentFinder finished");

			timeStart = System.currentTimeMillis();
			fragmentFinderSmarts(moleculeSmiles, fragments, directory, challenge);
			timeEnd = System.currentTimeMillis();
			timeFragmentFinderSmarts += (timeEnd - timeStart);
			System.out.println("FragmentFinderSmarts finished");
			
			timeStart = System.currentTimeMillis();
			fragmentFinderWithoutAromaticity(moleculeSmiles, fragments, directory, challenge);
			timeEnd = System.currentTimeMillis();
			timeFragmentFinderWithoutAromaticity += (timeEnd - timeStart);
			System.out.println("FragmentFinderWithoutAromaticity finished");
		}
		
		PrintWriter timesWriter = new PrintWriter(directory + File.separator + "results-mit-mz" + File.separator + "times.txt", "UTF-8");
		timesWriter.println("FragmentFinder                   " + timeFragmentFinder + "ms");
		timesWriter.println("FragmentFinderSmarts             " + timeFragmentFinderSmarts + "ms");
		timesWriter.println("FragmentFinderWithoutAromaticity " + timeFragmentFinderWithoutAromaticity + "ms");
		timesWriter.close();		
		System.out.println("FINISHED");
	}
	
	private static List<String> getChallenges(String csvFilePath) throws IOException {
		FileReader metaDataReader = new FileReader(csvFilePath);
		BufferedReader metaDataRowReader = new BufferedReader(metaDataReader);
		
		List<String> challenges = Lists.newArrayList();
		
		String row = metaDataRowReader.readLine();
		while(row != null) {
			String[] columns = row.split("[|]");
			String smiles = columns[1];
			challenges.add(smiles);
			row = metaDataRowReader.readLine();
		}
		
		metaDataReader.close();
		metaDataRowReader.close();
		return challenges;
	}
	
	private static List<String> mzValues(String inputFile) throws IOException {
		FileReader mzReader = new FileReader(inputFile);
		BufferedReader mzRowReader = new BufferedReader(mzReader);
		
		List<String> mzValues = Lists.newArrayList();
		
		String row = mzRowReader.readLine();
		while(row != null) {
			String[] columns = row.split(";");
			String value = columns[1];
			mzValues.add(value);
			row = mzRowReader.readLine();
		}
		
		mzReader.close();
		mzRowReader.close();
		return mzValues;
	}

	private static Map<String, String> moleculeSmiles(String csvFilePath) throws IOException {
		FileReader challengeReader = new FileReader(csvFilePath);
		BufferedReader challengeRowReader = new BufferedReader(challengeReader);
		
		Map<String, String> molecules = new HashMap<String, String>();
		
		String row = challengeRowReader.readLine();
		row = challengeRowReader.readLine();
		while(row != null) {
			String[] columns = row.split("\"");
			String smiles = columns[5];
			String name = columns[1];
			molecules.put(smiles, name);
			row = challengeRowReader.readLine();
		}
		
		challengeReader.close();
		challengeRowReader.close();
		return molecules;
	}
	
	private static void fragmentSmarts(ListMultimap<String, String> fragmentSmarts, String inputFilePath) throws IOException {
		FileReader inputFileReader = new FileReader(inputFilePath);
		BufferedReader inputFileRowReader = new BufferedReader(inputFileReader);
		
		String row = inputFileRowReader.readLine();
		while(row != null) {
			String[] columns = row.split(";");
			String fragment = columns[0];
			String mzValue = columns[1];
			fragmentSmarts.put(mzValue, fragment);
			row = inputFileRowReader.readLine();
		}
		
		inputFileReader.close();
		inputFileRowReader.close();
	}
	
	private static void fragmentFinder(Map<String, String> moleculeSmiles, List<Collection<String>> fragmentSmarts, String directory, String challenge) throws FileNotFoundException, UnsupportedEncodingException, CDKException {
		PrintWriter writer = new PrintWriter(directory + File.separator + "results" + File.separator + challenge + "-FragmentFinder.csv", "UTF-8");
		for(String molecule : moleculeSmiles.keySet()) {
			molecule = molecule.replace("\"", "");
			boolean allFound = true;
			ListMultimap<String, Integer> multimap = ArrayListMultimap.create();
			for(Collection<String> fragment : fragmentSmarts) {
				int found = 0;
				for(String frag : fragment) {
					found += FragmentFinder.findFragments(molecule, frag);
					multimap.put(frag, found);
				}
				if(found == 0) {
					allFound = false;
					break;
				}
			}
			if(allFound) {
				writer.println(molecule);
				writer.println(moleculeSmiles.get(molecule));
				for(String fragment : multimap.keySet()) {
					List<Integer> counts = multimap.get(fragment);
					for(Integer count : counts) {
						writer.println(";" + fragment + ";" + count);
					}
				}
				writer.println();
			}
		}
		writer.close();
	}
	
	private static void fragmentFinderSmarts(Map<String, String> moleculeSmiles, List<Collection<String>> fragmentSmarts, String directory, String challenge) throws Exception {
		PrintWriter writer = new PrintWriter(directory + File.separator + "results" + File.separator + challenge + "-FragmentFinderSmarts.csv", "UTF-8");
		for(String molecule : moleculeSmiles.keySet()) {
			boolean allFound = true;
			ListMultimap<String, Integer> multimap = ArrayListMultimap.create();
			for(Collection<String> fragment : fragmentSmarts) {
				int found = 0;
				for(String frag : fragment) {
					found += FragmentFinderSmarts.findFragments(molecule, frag);
					multimap.put(frag, found);
				}
				if(found == 0) {
					allFound = false;
					break;
				}
			}
			if(allFound) {
				writer.println(molecule);
				writer.println(moleculeSmiles.get(molecule));
				for(String fragment : multimap.keySet()) {
					List<Integer> counts = multimap.get(fragment);
					for(Integer count : counts) {
						writer.println(";" + fragment + ";" + count);
					}
				}
				writer.println();
			}
		}
		
		writer.close();
	}
	
	private static void fragmentFinderWithoutAromaticity(Map<String, String> moleculeSmiles, List<Collection<String>> fragmentSmarts, String directory, String challenge) throws FileNotFoundException, UnsupportedEncodingException, CDKException {
		PrintWriter writer = new PrintWriter(directory + File.separator + "results" + File.separator + challenge + "-FragmentFinderWithoutAromaticity.csv", "UTF-8");
		for(String molecule : moleculeSmiles.keySet()) {
			boolean allFound = true;
			ListMultimap<String, Integer> multimap = ArrayListMultimap.create();
			for(Collection<String> fragment : fragmentSmarts) {
				int found = 0;
				for(String frag : fragment) {
					found += FragmentFinderWithoutAromaticity.findFragments(molecule, frag);
					multimap.put(frag, found);
				}
				if(found == 0) {
					allFound = false;
					break;
				}
			}
			if(allFound) {
				writer.println(molecule);
				writer.println(moleculeSmiles.get(molecule));
				for(String fragment : multimap.keySet()) {
					List<Integer> counts = multimap.get(fragment);
					for(Integer count : counts) {
						writer.println(";" + fragment + ";" + count);
					}
				}
				writer.println();
			}
		}
		writer.close();
	}
}
