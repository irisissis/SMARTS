package implementation;

import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;

import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;

/*
 * method to find fragments in molecule over their SMILES/SMARTS
 * input molecule can only handle as SMILES
 * input fragment must be generated with the FragmentBuilder
 * create an output file output_fragment_finder.csv with contains the tested molecules, the tested fragments and the detected mappings 
 */
public class FragmentFinderWithoutAromaticity {
	public static int findFragments(String moleculeSmiles, String fragmentSmiles) throws CDKException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		sp.kekulise(false);
		
		IAtomContainer molecule = MoleculeFunctions.parseSmilesImplicitHydrogen(moleculeSmiles);
		
	    Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		
		SmilesParser fragmentparser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		//kelkulise is false, because the fragment SMILES have an user defined form 
		fragmentparser.kekulise(false);

		//convert the fragment SMILES in IAtomContainer and set the bonds order, because is missing of the SMILES form
		IAtomContainer fragment = fragmentparser.parseSmiles(fragmentSmiles);
		
		if(!hasIons(molecule) && hasIons(fragment)) {
			return 0;
		}
		
		//test is fragment SMILES subgraph of the molecule SMILES
		Isomorphism comparison = new Isomorphism(Algorithm.SubStructure, true);
		// set molecules (fragment, molecule), remove hydrogens: true, clean and configure molecule: false
		comparison.init(fragment, molecule, true, false);
		// set chemical filter all false
		comparison.setChemFilters(false, false, false);
		
		//List with all detected mappings
		List<Map<Integer, Integer>> mappings = comparison.getAllMapping();
		//no mappings detected
		if(mappings == null) {
			return 0;
		}
		
		return mappings.size();
	}
	
	private static boolean hasIons(IAtomContainer molecule) {
		for(IAtom atom : molecule.atoms()) {
			if(atom.getFormalCharge() != 0) {
				return true;
			}
		}
		return false;
	}
}
