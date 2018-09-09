package implementation;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class FragmentFinderSmarts {
	public static int findFragments(String smiles, String smarts) throws Exception {
		SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		IAtomContainer molecule = smilesParser.parseSmiles(smiles);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		SMARTSQueryTool smartsQueryTool = new SMARTSQueryTool(replaceHydrogens(smarts), DefaultChemObjectBuilder.getInstance());
		smartsQueryTool.setAromaticity(new Aromaticity(ElectronDonation.cdk(), Cycles.all()));
		smartsQueryTool.matches(molecule);
		return smartsQueryTool.countMatches();
	}
	
	private static String replaceHydrogens(String smarts) {
		String regex = "\\[[a-zA-Z]H((\\d\\+{0,1}\\])|(\\+{0,1}\\]))";
		Pattern pattern = Pattern.compile(regex);
		Matcher matcher = pattern.matcher(smarts);
		while(matcher.find()) {
			String match = matcher.group();
			String atomSymbol = String.valueOf(match.charAt(1));
			smarts = match.contains("+") ? matcher.replaceFirst('[' + atomSymbol + "+]") : matcher.replaceFirst(atomSymbol);
			matcher = pattern.matcher(smarts);
		}
		return smarts;
	}
}
