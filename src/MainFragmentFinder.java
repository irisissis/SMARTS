package implementation;

public class MainFragmentFinder {
	/*
	 * Parameter (1): path to base directory
	 * Parameter (2): "1" - use FragmentFinderDatabase
	 *                "2" - use FragmentFinderMZDatabase
	 */
	public static void main(String... args) throws Exception {
		String baseDirectory = args[0];
		String modus = args[1];
		switch(modus) {
		case "1": FragmentFinderDatabase.run(baseDirectory); break;
		case "2": FragmentFinderMZDatabase.run(baseDirectory); break;
		}
	}
}
