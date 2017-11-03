import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import java.lang.*;
public class Utils {

    /**
     * Returns the minimum value per dimension in a matrix
     * @param data
     * @return
     */
    public double[] mins(double[][] data){
        double mins[] = new double[data[0].length];
        for(int dim = 0; dim < data[0].length; dim++){
            mins[dim] = data[0][dim];
            for (int patient = 0; patient < data.length; patient++){
                if(data[patient][dim] < mins[dim]){
                    mins[dim] = data[patient][dim];
                }
            }
        }
        return mins;
    }

    /**
     * Returns the maximum value per dimension in a matrix
     * @param data
     * @return
     */
    public double[] maxs(double[][] data){
        double maxs[] = new double[data[0].length];
        for(int dim = 0; dim < data[0].length; dim++){
            maxs[dim] = data[0][dim];
            for (int patient = 0; patient < data.length; patient++){
                if(data[patient][dim] > maxs[dim]){
                    maxs[dim] = data[patient][dim];
                }
            }
        }
        return maxs;
    }

    public double random(double min, double max) {
        return ThreadLocalRandom.current().nextDouble(min, max);

    }

    public double[] distances(double[] data, double[][] centers){
        double[] distances = new double[centers.length];
        for(int i = 0; i < centers.length; i++){
            distances[i] = distance(data, centers[i]);
        }
        return distances;
    }

    public double distance(double[] data, double[] center){
        double distance = 0;
        for(int i = 0; i < data.length; i++)
        {
            distance += Math.pow( data[i] - center[i], 2);
        }
        return Math.sqrt(distance);
    }

    public int argMin(double[] data){
        double bestValue = data[0];
        int bestIndex = 0;
        for(int i = 0; i < data.length; i++){
            if(data[i] < bestValue)
            {
                bestValue = data[i];
                bestIndex = i;
            }

        }
        return bestIndex;
    }

    public double mean(ArrayList<Double> data){
        double avg = 0;
        for (int i = 0; i < data.size(); i++){
            avg += data.get(i);
        }
        avg /= data.size();
        return avg;
    }
}
