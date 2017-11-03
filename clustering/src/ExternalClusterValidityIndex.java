/**
 *
 */
public abstract class ExternalClusterValidityIndex extends ClusterValidityIndex {

	/**
	 * 
	 * 
	 * @param clustering
	 *            The cluster id of object i is stored in the clustering array
	 *            at position i-1.
	 * @param goldStandard
	 *            The class id of object i is stored in the goldStandard array
	 *            at position i-1.
	 * @return A validity for the given clustering.
	 */
	public abstract double getValidity(int[] clustering, int[] goldStandard);
}
