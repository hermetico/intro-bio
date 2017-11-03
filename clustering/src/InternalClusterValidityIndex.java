/**
 *
 */
public abstract class InternalClusterValidityIndex extends ClusterValidityIndex {

	/**
	 * 
	 * 
	 * @param clustering
	 *            The cluster id of object i is stored in the clustering array
	 *            at position i-1.
	 * @param similarities
	 *            The similarity between object i and j are stored at position
	 *            [i-1,j-1] and [j-1,i-1] (since we assume symmetric
	 *            simmilarities.
	 * @return A validity for the given clustering.
	 */
	public abstract double getValidity(int[] clustering, double[][] similarities);
}
