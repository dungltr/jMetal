package org.uma.jmetal.util.experiment.util;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.experiment.Experiment;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.fileoutput.writeCSV;
//.nsgav.utilsnsgav.writeCSV;
import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Class defining tasks for the execution of algorithms in parallel.
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class ExperimentAlgorithm<S extends Solution<?>, Result>  {
  private Algorithm<Result> algorithm;
  private String algorithmTag;
  private String problemTag;

  /**
   * Constructor
   */
  public ExperimentAlgorithm(
          Algorithm<Result> algorithm,
          String algorithmTag,
          String problemTag) {
    this.algorithm = algorithm;
    this.algorithmTag = algorithmTag;
    this.problemTag = problemTag;
  }

  public ExperimentAlgorithm(
          Algorithm<Result> algorithm,
          String problemTag) {
    this(algorithm, algorithm.getName(), problemTag) ;
  }

  public void runAlgorithm(int id, Experiment<?, ?> experimentData) {
    String outputDirectoryName = experimentData.getExperimentBaseDirectory()
            + "/data/"
            + algorithmTag
            + "/"
            + problemTag;

    File outputDirectory = new File(outputDirectoryName);
    if (!outputDirectory.exists()) {
      boolean result = new File(outputDirectoryName).mkdirs();
      if (result) {
        JMetalLogger.logger.info("Creating " + outputDirectoryName);
      } else {
        JMetalLogger.logger.severe("Creating " + outputDirectoryName + " failed");
      }
    }

    String funFile = outputDirectoryName + "/FUN" + id + ".tsv";
    String varFile = outputDirectoryName + "/VAR" + id + ".tsv";
    JMetalLogger.logger.info(
            " Running algorithm: " + algorithmTag +
                    ", problem: " + problemTag +
                    ", run: " + id +
                    ", funFile: " + funFile);


    algorithm.run();
    ////////Dung edit//////
    AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
            .execute() ;
    long computingTime = algorithmRunner.getComputingTime() ;
    double[] array = new double[2];
    array[0] = (double)id;
    array[1] = (double)computingTime;
    JMetalLogger.logger.info("Total execution time of "+algorithm.getName()+": " + computingTime + "ms");
    try {
      writeCSV.addArray2Csv(ReadFile.readhome("HOME")+algorithm.getName()+"_"+ problemTag +".csv",array);
      writeCSV.addNumberCsv(ReadFile.readhome("HOME")+algorithm.getName()+"_"+ problemTag+"_computeTime" +".csv",(double)computingTime);
    } catch (IOException e) {
      e.printStackTrace();
    }
    ///////////////////////////////////////////////////////////////////////////////
    Result population = algorithm.getResult();

    new SolutionListOutput((List<S>) population)
            .setSeparator("\t")
            .setVarFileOutputContext(new DefaultFileOutputContext(varFile))
            .setFunFileOutputContext(new DefaultFileOutputContext(funFile))
            .print();
  }

  public Algorithm<Result> getAlgorithm() {
    return algorithm;
  }

  public String getAlgorithmTag() {
    return algorithmTag;
  }

  public String getProblemTag() {
    return problemTag;
  }
}
