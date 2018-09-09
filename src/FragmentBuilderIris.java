package implementation;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import com.google.common.collect.Lists;


/*
 * function to convert a fingerprint in a fragment SMILES/SMARTS
 * needed fingerprint and related IAtomContainer as input
 * return a fragment SMILES/SMARTS and fingerprint
 * 
 * 
 */
public class FragmentBuilderIris {
	/*
	 * parameter:
	 * - molecule: IAtomContainer input
	 * - fingerprint: fingerprint for fragment from input
	 * - writer: file writer for output, null for printing on terminal
	 * - printTrace: boolean flag - true: print whole program stack trace, false: print nothing
	 */
	public static String printFragmentToFile(IAtomContainer molecule, String fingerprint, String origSmiles, String metFragFragment, PrintWriter writer, boolean printTrace) throws CDKException {
		AtomContainerManipulator manipulator = new AtomContainerManipulator();
		manipulator.suppressHydrogens(molecule);
		
		//test if fingerprint length equal to number of atoms from molecule
		// if to short or to long than print ERROR and stop
		if(fingerprint.length() != molecule.getAtomCount()) {
			print(writer, "ERROR: fingerprint has invalid length. Number of atoms: " + molecule.getAtomCount() + " fingerprint length: " + fingerprint.length(), printTrace);
			return null;
		}
    
		//test that fingerprint contains at least one atom
		//if not print ERROR and stop
		boolean fragmentExist = true;
		for(int i = 0; i < fingerprint.length(); i++) {
    	if(fingerprint.charAt(i) == '1') {
    		fragmentExist = false;
    		break;
    	}
    }
    if(fragmentExist) {
    	print(writer, "ERROR: no fragment found.", printTrace);
    	return null;
    }

    //test that molecule contains a ring
    RingSearch ringSearch = new RingSearch(molecule);
    int[] cyclic = ringSearch.cyclic();
    
    // aronmatize molecule
    Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
	aromaticity.apply(molecule);

    //IAtomContainer for fragment
    AtomContainer fragment = new AtomContainer();
    
    // foreach atom in molecule
    //test is atom part from fragment and if true: copy atom to fragment 
    for(int i = 0; i < molecule.getAtomCount(); i++) {
    	if(fingerprint.toCharArray()[i] == '1') {
    		fragment.addAtom(molecule.getAtom(i));
    	}
    }
    
    // foreach bond in molecule 
    // test that both atoms connected to  the bond are part of the fragment, if true: add bond to fragment   
    for(IBond bond : molecule.bonds()) {
    	if(fragmentContainsAtom(fragment, bond.getAtom(0)) && fragmentContainsAtom(fragment, bond.getAtom(1))) {
    		fragment.addBond(bond);
    	}
    }

    //set SmilesGenerator null, because test if linear or ring molecule
    //if linear: create non-canonical aromatic SMILES
    //if ring: create canonical aromatic SMILES
    SmilesGenerator fragmentgenerator = null;
    if(cyclic.length != 0) {
    	fragmentgenerator = new SmilesGenerator(SmiFlavor.Absolute | SmiFlavor.UseAromaticSymbols);
    } else {
    	fragmentgenerator = new SmilesGenerator(SmiFlavor.Generic | SmiFlavor.UseAromaticSymbols);
    }
    
    //create fragment SMILES
    String fragmentSmile = fragmentgenerator.create(fragment);//.replace("[", "").replace("]", "");
     
    //if there more than one fragment, get Warning and print all fragments
    if(fragmentSmile.contains(".")) {
    	print(writer, "WARNING: too many fragments found for fingerprint: " + fingerprint, printTrace);
    }

    return fragmentSmile;
	}
	
	// test: if atom part of fragment 
	private static boolean fragmentContainsAtom(AtomContainer fragment, IAtom atom) {
		for(IAtom fragmentAtom : fragment.atoms()) {
			if(fragmentAtom.equals(atom)) {
				return true;
			}
		}
		return false;
	}
	
	//method if writer handed over or write in terminal
	private static void print(PrintWriter writer, String message, boolean printTrace) {
		if(writer == null) {
			if(printTrace) {
				System.out.println(message);
			}
			return;
		}
		writer.println(message);
	}
}
