package implementation;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;

import com.google.common.collect.Lists;

import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;


/*
 * method to find fragments in molecule over their SMILES/SMARTS
 * input molecule can only handle as SMILES
 * input fragment must be generated with the FragmentBuilder
 * create an output file output_fragment_finder.csv with contains the tested molecules, the tested fragments and the detected mappings 
 */
public class FragmentFinder {
	public static int findFragments(String moleculeSmiles, String fragmentSmiles) throws CDKException {
		//System.out.println(moleculeSmiles);
//		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//		sp.kekulise(false);
//		IAtomContainer molecule = sp.parseSmiles(moleculeSmiles);
		IAtomContainer molecule = null;
		try {
			molecule = MoleculeFunctions.parseSmilesImplicitHydrogen(moleculeSmiles);
		} catch (Exception e) {
			System.out.println(moleculeSmiles);
			return 0;
		}

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

		//delete all mappings where the aromaticity of the bonds and atoms are different
		// eliminate mappings from aromatic fragments in aliphatic structures 
		List<Integer> indexesToDelete = Lists.newArrayList();
		
		for(int i = 0; i < mappings.size(); i++) {
			if(mappings.get(i).size() != fragment.getAtomCount()) {
				indexesToDelete.add(i);
			}
		}
		Collections.reverse(indexesToDelete);
		for(Integer index : indexesToDelete) {
			mappings.remove(index.intValue());
		}
		
		indexesToDelete = Lists.newArrayList();
		
		//foreach mapping, the indices of atoms are used, because the mapping list make new IAtomContainer
//		System.out.println(moleculeSmiles + " " + fragmentSmiles + " " + mappings.size());
		
		for(int i = 0; i < mappings.size(); i++) {
			Map<Integer, Integer> mapping = mappings.get(i);
			boolean deleteMapping = false;
			Set<Integer> fragmentAtoms = mapping.keySet();
			for(Integer fragmentAtom1 : fragmentAtoms) {
				Integer moleculeAtom1 = mapping.get(fragmentAtom1);
				if(!atomsAreCorrect(fragment.getAtom(fragmentAtom1), molecule.getAtom(moleculeAtom1))) {
//					System.out.println(fragment.getAtom(fragmentAtom1).isAromatic() + " " + molecule.getAtom(moleculeAtom1).isAromatic());
					deleteMapping = true;
					break;
				}
				
				//System.out.println(fragmentAtom1.hashCode());
				//get bond order to test, that the right mapping is detected, because algorithm only checked the element symbols
				//List which atom is connected to which atom in fragment
				List<IAtom> connectedFragmentAtoms = fragment.getConnectedAtomsList(fragment.getAtom(fragmentAtom1));
				//foreach atom in connectedFragmentAtoms search bond order and check is the same as in fragment 
				for(IAtom fragmentAtom2 : connectedFragmentAtoms) {
					IBond fragmentBond = fragment.getBond(fragment.getAtom(fragmentAtom1), fragmentAtom2);
					Integer moleculeAtom2 = mapping.get(fragment.indexOf(fragmentAtom2));
					IBond moleculeBond = molecule.getBond(molecule.getAtom(moleculeAtom1), molecule.getAtom(moleculeAtom2));
					if(!fragmentBond.getOrder().equals(moleculeBond.getOrder()) && !fragmentBond.getOrder().equals(Order.UNSET)) {
//						System.out.println(fragmentBond.getOrder() + " " + moleculeBond.getOrder());
						deleteMapping = true;
						break;
					}
				}
				if(deleteMapping) {
					break;
				}
			}
			//delete false mappings
			if(deleteMapping) {
				indexesToDelete.add(i);
			}
		}
		
		//delete false detected mappings
		Collections.reverse(indexesToDelete);
		for(Integer index : indexesToDelete) {
			mappings.remove(index.intValue());
		}
		
		indexesToDelete = Lists.newArrayList();
		List<List<Integer>> moleculeIndexLists = Lists.newArrayList(Lists.newArrayList());
		
		//test to eliminate the same mappings, because the algorithm has no reading direction
		for(int i = 0; i < mappings.size(); i++) {
			Map<Integer, Integer> mapping = mappings.get(i);
			List<Integer> values = Lists.newArrayList(mapping.values());
			Collections.sort(values);
			if(!moleculeIndexLists.contains(values)) {
				moleculeIndexLists.add(values);
			} else {
				indexesToDelete.add(i);
			}
		}
		
		//System.out.println(mappings);
		//delete false detected mappings
		Collections.reverse(indexesToDelete);
		for(Integer index : indexesToDelete) {
			mappings.remove(index.intValue());
		}
		
		//System.out.println(mappings);
		
//		System.out.println("Matching molecule indices:");
//		for(Map<Integer, Integer> mapping : mappings) {
//			System.out.print("- " + mapping.values());
//			for(Integer index : mapping.values()) {
//				System.out.print(" " + molecule.getAtom(index).getSymbol());
//			}
//			System.out.println();
//		}
		//return the detected mapping and how often for one molecule
//		System.out.println(mappings.size());
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
	
	private static boolean atomsAreCorrect(IAtom fragmentAtom, IAtom moleculeAtom) {
		boolean isAromatic = fragmentAtom.isAromatic() == moleculeAtom.isAromatic();
		boolean chargesAreEqual = fragmentAtom.getFormalCharge() == moleculeAtom.getFormalCharge();
		boolean hydrogenCountsAreEqual = true;
//		if(fragmentAtom.getSymbol() != "C" && (fragmentAtom.getFormalCharge() != 0 || moleculeAtom.getFormalCharge() != 0)) {
//			hydrogenCountsAreEqual = fragmentAtom.getImplicitHydrogenCount() == moleculeAtom.getImplicitHydrogenCount();
//		}
		return isAromatic && chargesAreEqual && hydrogenCountsAreEqual;
	}
}
