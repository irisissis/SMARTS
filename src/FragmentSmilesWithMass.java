package implementation;

/*
 * Utility class for saving for a fragmentSmile its original SMILE from MetFrag and its mass 
 */
public class FragmentSmilesWithMass {
	private String fragmentSmiles;
	private String metFragFragmentSmiles;
	private Double mass;
	
	FragmentSmilesWithMass() {
		
	}
	
	FragmentSmilesWithMass(String fragmentSmiles, Double mass) {
		this.fragmentSmiles = fragmentSmiles;
		this.mass = mass;
	}
	
	public String getFragmentSmiles() {
		return fragmentSmiles;
	}
	
	public void setFragmentSmiles(String fragmentSmiles) {
		this.fragmentSmiles = fragmentSmiles;
	}
	
	public String getMetFragFragmentSmiles() {
		return metFragFragmentSmiles;
	}
	
	public void setMetFragFragmentSmiles(String metFragFragmentSmiles) {
		this.metFragFragmentSmiles = metFragFragmentSmiles;
	}
	
	public Double getMass() {
		return mass;
	}
	
	public void setMass(Double mass) {
		this.mass = mass;
	}
}
