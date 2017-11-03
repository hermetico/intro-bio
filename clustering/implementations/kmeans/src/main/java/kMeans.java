import java.util.ArrayList;

/**
 * This is a simple implementation of kmeans without nstart parameter.
 */
public class kMeans {

    protected Utils utils;
	protected int k;
	protected double[][] dataset;
	protected double[][] centroids;
    protected int[] assignations;
	protected int iterations;


	/**
	 * @param iterations
	 *            The number of total iterations to perform
	 * 
	 */
	public kMeans(int iterations) {
		this.utils = new Utils();
		this.iterations = iterations;
	}

	/**
	 * This method clusters the given dataset.
	 * 
	 * @param dataset
	 *            Row [i][] of the dataset contains measurements (e.g. gene
	 *            expression levels) of object (e.g. patient) i+1.
	 * @return The clustering of the dataset. The length of this array has to be
	 *         equal to the number of rows of the dataset. The cluster id of
	 *         object i is stored at position i-1 of this array.
	 */
	public int[] clusterData(double[][] dataset, int k) {
		this.k = k;
		this.dataset = dataset;
		this.centroids = new double[k][dataset[0].length];
		this.assignations = new int[dataset.length];

		seedRandomCentroids(dataset);
		int iterationsDone = 0;
		while (iterationsDone < this.iterations) {
			assignObjectsToClosestCentroids();
            //showCentroids();
			updateCentroids();
			//showCentroids();
            iterationsDone++;

		}
		return getClusterAssignment();
	}

	public void seedRandomCentroids_old(double[][] dataset) {
	    double[] mins = utils.mins(dataset);
	    double[] maxs = utils.maxs(dataset);

	    for(int cent = 0; cent < centroids.length; cent++)
        {
            int patient = (int)utils.random(0, dataset.length );
            for( int dim = 0; dim < mins.length; dim++){
                centroids[cent][dim] = utils.random(mins[dim], maxs[dim]);
            }
        }

	}
    public void seedRandomCentroids(double[][] dataset) {
        for (int cent = 0; cent < centroids.length; cent++) {
            int patient = (int) utils.random(0, dataset.length);
            for (int dim = 0; dim < centroids[0].length; dim++) {
                centroids[cent][dim] = dataset[patient][dim];
            }
        }
    }


    public void assignObjectsToClosestCentroids() {
	    for(int patient = 0; patient < dataset.length; patient++){
	        double[] distances = utils.distances(dataset[patient], centroids);
            this.assignations[patient] = utils.argMin(distances);
        }
	}

	public void updateCentroids() {
	    for(int cent = 0; cent < centroids.length; cent++){ // each centroid
            for(int dim = 0; dim < centroids[0].length; dim++){ // each centroid dimension
                ArrayList<Double> values = new ArrayList<>();
                // gather the dimension values of the patients associated to the centroid
                // in order to compute the average position for each dimension
                for(int patient = 0; patient < assignations.length; patient++){
                    if(assignations[patient] == cent){ // patients associated to that centroid
                        values.add( dataset[patient][dim] );
                    }
                }

                centroids[cent][dim] = utils.mean(values);
            }
        }
	}

	public void showCentroids(){
	    for(int i = 0; i < centroids.length; i++){
	        for(int j = 0; j < centroids[0].length; j++){
	            System.out.print(centroids[i][j] + " ");
            }
	        System.out.println();
        }
    }

	public int[] getClusterAssignment() {
		return assignations;
	}
}
