import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;

public class Main {

    @Option(name="-f", aliases = {"--file"}, usage="The input File", required = true)
    private String inputFileName;

    @Option(name="-o", aliases = {"--output"}, usage="The output file name", required = true)
    private String outputFileName;

    @Option(name="-k", usage="The number of k centroids (default: 3)")
    private int centroids = 3;

    @Option(name="-i", aliases = {"--iterations"}, usage="The numberof iterations (default: 20)")
    private int iterations = 20;

    private Map<String, Integer> ids_to_index;
    private Map<Integer, String> ids_from_index;
    private double[][] dataset;

    public static void main(String args[]){
        new Main(args);
    }

    public Main(String args[]){
        parseInputParameters(args);
        parseInputData(readInputFile());
        kMeans kmeans = new kMeans(iterations);
        int[] assignations = kmeans.clusterData(dataset.clone(), centroids);
        writeOutput(assignations);
    }

    public void parseInputParameters(String[] args){
        CmdLineParser parser = new CmdLineParser(this);
        parser.setUsageWidth(80); // width of the error display area

        try {
            // parse the arguments.
            parser.parseArgument(args);

          } catch( CmdLineException e ) {
            // if there's a problem in the command line,
            // you'll get this exception. this will report
            // an error message.
            System.err.println(e.getMessage());
            System.err.println("java intro-bio-kmeans [options...] arguments...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println();

            // print option sample. This is useful some time
            System.exit(1);

            return;
        }

    }

    public ArrayList<ArrayList<Double>> readInputFile(){
        ArrayList<ArrayList<Double>> data = new ArrayList<>();
        ids_to_index = new HashMap<>();
        ids_from_index = new HashMap<>();
        int index = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(this.inputFileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                String pieces[] = line.split("\t");
                String id = pieces[0];

                ArrayList<Double> values = new ArrayList<>(pieces.length - 1);
                for(int i = 1; i < pieces.length; i++){
                    values.add(Double.parseDouble(pieces[i]));
                }

                ids_to_index.put(pieces[0], index);
                ids_from_index.put(index, pieces[0]);

                data.add(values);
                index++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        return data;
    }

    public void parseInputData(ArrayList<ArrayList<Double>> data) {
        this.dataset = new double[data.size()][data.get(0).size()];
        int index = 0;
        for(ArrayList<Double> row : data){
            double[] newRow = new double[row.size()];
            for(int j = 0; j < row.size(); j++){
                newRow[j] = (double) row.get(j);
            }
            dataset[index] = newRow;
            index++;
        }
    }

    public void writeOutput(int[] assignations)
    {
        PrintWriter writer = null;
        String filename  = this.centroids + "_" + this.outputFileName;
        try {
            writer = new PrintWriter(filename);
            for(int centroid = 0; centroid < centroids; centroid++ )
            {
                writer.print(centroid);
                for(int id = 0; id < assignations.length; id++)
                {
                    if( assignations[id] == centroid) {
                        writer.print("\t" + ids_from_index.get(id));
                    }
                }
                writer.print("\n");
            }
            writer.close();
            System.out.println("Output stored at " + filename);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }
}
