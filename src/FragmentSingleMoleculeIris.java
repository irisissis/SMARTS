package implementation;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;
import de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException;
import de.ipbhalle.metfraglib.fragment.DefaultBitArrayFragment;
import de.ipbhalle.metfraglib.interfaces.IFragment;
import de.ipbhalle.metfraglib.interfaces.IMatch;
import de.ipbhalle.metfraglib.list.MatchList;
import de.ipbhalle.metfraglib.list.SortedScoredCandidateList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.process.CombinedMetFragProcess;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;

public class FragmentSingleMoleculeIris {
	public static List<FragmentSmilesWithMass> getFragment(String smiles, Double[] masses, Double[] ints, double exactMass, PrintWriter writer) throws Exception {
		
		/**
		 * 1st: set all parameters and inputs
		 */
		// set metfrag parameters
		double mzppm = 0;
		double mzabs = 0.5;
		boolean posCharge = true;
		int mode = 1;
		int treeDepth = 2;
		
		// init molecule as atomcontainer
		IAtomContainer molecule = null;
		try {
			molecule = MoleculeFunctions.parseSmiles(smiles);
		} catch (InvalidSmilesException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			System.exit(1);
		}
		
		// set candidate list as needed by metfrag
		IAtomContainer[] moleculeAsArray = {molecule};

		// generate peak list string as needed by metfrag
		String peaksString = "";
		if(masses.length > 0) peaksString += masses[0] + " " + ints[0];
		for(int i = 1; i < masses.length; i++) {
			peaksString += "\n" + masses[i] + " " + ints[i];
		}
		
		// set scoring parameters
		String[] score_names = {VariableNames.METFRAG_FRAGMENTER_SCORE_NAME};
		Double[] score_weights = {1.0};
		
		// init metfrag settings object
		MetFragGlobalSettings settings = new MetFragGlobalSettings();
		settings.set(VariableNames.MOLECULES_IN_MEMORY, moleculeAsArray);
		settings.set(VariableNames.PEAK_LIST_STRING_NAME, peaksString);
		settings.set(VariableNames.METFRAG_DATABASE_TYPE_NAME, "LocalInMemoryDatabase");
		settings.set(VariableNames.METFRAG_PEAK_LIST_READER_NAME, "de.ipbhalle.metfraglib.peaklistreader.FilteredStringTandemMassPeakListReader");
		settings.set(VariableNames.METFRAG_SCORE_TYPES_NAME, score_names);
		settings.set(VariableNames.METFRAG_SCORE_WEIGHTS_NAME, score_weights);
		settings.set(VariableNames.RELATIVE_MASS_DEVIATION_NAME, mzppm);
		settings.set(VariableNames.ABSOLUTE_MASS_DEVIATION_NAME, mzabs);
		settings.set(VariableNames.IS_POSITIVE_ION_MODE_NAME, posCharge);
		settings.set(VariableNames.PRECURSOR_ION_MODE_NAME, mode);
		settings.set(VariableNames.PRECURSOR_NEUTRAL_MASS_NAME, exactMass);
		settings.set(VariableNames.MAXIMUM_TREE_DEPTH_NAME, (byte)treeDepth);
		settings.set(VariableNames.USE_SMILES_NAME, true);
		
		/**
		 * 2nd: run the MetFrag fragmentation 
		 */
		// init metfrag process
		CombinedMetFragProcess mp = new CombinedMetFragProcess(settings);

		// retrieve candidates in this case from the "InMemoryDatabase" as set previously
		try {
			mp.retrieveCompounds();
		} catch (Exception e2) {
			System.err.println("Error retrieving candidates");
		}
		
		// run metfrag
		try {
			mp.run();
		} catch (Exception e) {
			System.err.println("Error running MetFrag process");
		}
		
		/**
		 * 3rd: get the results and print the output
		 */
		
		// get all scored candidates (here only benzene of course)
		SortedScoredCandidateList scoredCandidateList = (SortedScoredCandidateList) mp.getCandidateList();
		
		// get the match list (each match includes the matched fragments)
		// peak - fragments
		MatchList assignedFragmentList = scoredCandidateList.getElement(0).getMatchList();
		// init array for all fragments
		IAtomContainer[] assignedFragments = new IAtomContainer[assignedFragmentList.getNumberElements()];
		
		List<FragmentSmilesWithMass> fragmentSmiles = Lists.newArrayList();
		// iterate over all matches and extract the fingerprints and fragments
		for(int i = 0; i < assignedFragmentList.getNumberElements(); i++) {
			// get the match
			IMatch match = assignedFragmentList.getElement(i);
			// only get the best fragment of the match
			IFragment fragment = match.getBestMatchedFragment();
			// init the atomcontainer
			try {
				scoredCandidateList.getElement(0).setUseSmiles(true);
				scoredCandidateList.getElement(0).initialisePrecursorCandidate();
			} catch (AtomTypeNotKnownFromInputListException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			IAtomContainer con = scoredCandidateList.getElement(0).getPrecursorMolecule().getStructureAsIAtomContainer();
			
			IAtomContainer currentFragmentAC = fragment.getStructureAsIAtomContainer(scoredCandidateList.getElement(0).getPrecursorMolecule());
			assignedFragments[i] = currentFragmentAC;

			String fp = ((DefaultBitArrayFragment)fragment).getAtomsFastBitArray().toString();

			String metFragFragment = fragment.getSmiles(scoredCandidateList.getElement(0).getPrecursorMolecule());
			
			FragmentSmilesWithMass fragmentSmilesWithMass = new FragmentSmilesWithMass();
			fragmentSmilesWithMass.setFragmentSmiles(FragmentBuilderIris.printFragmentToFile(con, fp, smiles, metFragFragment, writer, true));
			fragmentSmilesWithMass.setMetFragFragmentSmiles(metFragFragment);
			fragmentSmilesWithMass.setMass(masses[i]);
			fragmentSmiles.add(fragmentSmilesWithMass);
		}
		for(FragmentSmilesWithMass fragment : fragmentSmiles) {
			writer.println(fragment.getFragmentSmiles());
			System.out.println(fragment.getFragmentSmiles() + " " + fragment.getMetFragFragmentSmiles());
		}
		return fragmentSmiles;
	}
}
