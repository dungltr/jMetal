package org.uma.jmetal.experiment;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.NSGAIIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgav.NSGAVBuilder;
import org.uma.jmetal.algorithm.multiobjective.smpso.SMPSOBuilder;
import org.uma.jmetal.algorithm.multiobjective.spea2.SPEA2Builder;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.problem.multiobjective.zdt.*;
import org.uma.jmetal.problem.multiobjective.dtlz.*;
import org.uma.jmetal.problem.multiobjective.UF.*;
import org.uma.jmetal.qualityindicator.impl.*;
import org.uma.jmetal.qualityindicator.impl.hypervolume.PISAHypervolume;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;
import org.uma.jmetal.util.archive.impl.CrowdingDistanceArchive;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.experiment.Experiment;
import org.uma.jmetal.util.experiment.ExperimentBuilder;
import org.uma.jmetal.util.experiment.component.*;
import org.uma.jmetal.util.experiment.util.ExperimentAlgorithm;
import org.uma.jmetal.util.experiment.util.ExperimentProblem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Example of experimental study based on solving the ZDT problems with the algorithms NSGAII,
 * SPEA2, and SMPSO
 *
 * This experiment assumes that the reference Pareto front are known, so the names of files
 * containing them and the directory where they are located must be specified.
 *
 * Six quality indicators are used for performance assessment.
 *
 * The steps to carry out the experiment are: 1. Configure the experiment 2. Execute the algorithms
 * 3. Compute que quality indicators 4. Generate Latex tables reporting means and medians 5.
 * Generate R scripts to produce latex tables with the result of applying the Wilcoxon Rank Sum Test
 * 6. Generate Latex tables with the ranking obtained by applying the Friedman test 7. Generate R
 * scripts to obtain boxplots
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */

public class DTLZStudy {
    private static final int INDEPENDENT_RUNS = 25;

    public static void main(String[] args) throws IOException {
        /*if (args.length != 1) {
            throw new JMetalException("Missing argument: experimentBaseDirectory");
        }*/
        String experimentBaseDirectory = ReadFile.readhome("HOME")+"/jMetalData";//args[0];

        List<ExperimentProblem<DoubleSolution>> problemList = new ArrayList<>();
        problemList.add(new ExperimentProblem<>(new DTLZ1()));
        problemList.add(new ExperimentProblem<>(new DTLZ2()));
        problemList.add(new ExperimentProblem<>(new DTLZ3()));
        /*problemList.add(new ExperimentProblem<>(new DTLZ4()));
        problemList.add(new ExperimentProblem<>(new DTLZ5()));
        problemList.add(new ExperimentProblem<>(new DTLZ6()));
        problemList.add(new ExperimentProblem<>(new DTLZ7()));
        problemList.add(new ExperimentProblem<>(new UF1()));
        problemList.add(new ExperimentProblem<>(new UF2()));
        problemList.add(new ExperimentProblem<>(new UF3()));
        /*problemList.add(new ExperimentProblem<>(new UF4()));
        problemList.add(new ExperimentProblem<>(new UF5()));
        problemList.add(new ExperimentProblem<>(new UF6()));
        problemList.add(new ExperimentProblem<>(new UF7()));
        problemList.add(new ExperimentProblem<>(new UF8()));
        problemList.add(new ExperimentProblem<>(new UF9()));
        problemList.add(new ExperimentProblem<>(new UF10()));

        problemList.add(new ExperimentProblem<>(new ZDT1()));
        problemList.add(new ExperimentProblem<>(new ZDT2()));
        problemList.add(new ExperimentProblem<>(new ZDT3()));
        */
        List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> algorithmList =
                configureAlgorithmList(problemList);

        List<String> referenceFrontFileNames =
                Arrays.asList("DTLZ1.8D.pf", "DTLZ2.8D.pf", "DTLZ3.8D.pf");//,"ZDT1.pf","ZDT2.pf","ZDT3.pf"
                        //"UF1.pf", "UF2.pf", "UF3.pf"// "DTLZ4.3D.pf", "DTLZ5.3D.pf","DTLZ6.3D.pf", "DTLZ7.3D.pf");
                //,"UF1.pf", "UF2.pf", "UF3.pf");//, "UF4.pf", "UF5.pf", "UF6.pf", "UF7.pf", "UF8.pf", "UF9.pf", "UF10.pf");

        Experiment<DoubleSolution, List<DoubleSolution>> experiment =
                new ExperimentBuilder<DoubleSolution, List<DoubleSolution>>("NSGAII and NSGAV with DTLZ 8D")
                        .setAlgorithmList(algorithmList)
                        .setProblemList(problemList)
                        .setReferenceFrontDirectory("/pareto_fronts")
                        .setReferenceFrontFileNames(referenceFrontFileNames)
                        .setExperimentBaseDirectory(experimentBaseDirectory)
                        .setOutputParetoFrontFileName("FUN")
                        .setOutputParetoSetFileName("VAR")
                        .setIndicatorList(Arrays.asList(
                                //new Epsilon<DoubleSolution>(),
                                //new Spread<DoubleSolution>(),
                                new GenerationalDistance<DoubleSolution>(),
                                new PISAHypervolume<DoubleSolution>()))//,
                                //new InvertedGenerationalDistance<DoubleSolution>(),
                                //new InvertedGenerationalDistancePlus<DoubleSolution>()))
                        .setIndependentRuns(INDEPENDENT_RUNS)
                        .setNumberOfCores(8)
                        .build();

        new ExecuteAlgorithms<>(experiment).run();
        new ComputeQualityIndicators<>(experiment).run();
        new GenerateLatexTablesWithStatistics(experiment).run();
        new GenerateWilcoxonTestTablesWithR<>(experiment).run();
        new GenerateFriedmanTestTables<>(experiment).run();
        new GenerateBoxplotsWithR<>(experiment).setRows(3).setColumns(3).setDisplayNotch().run();
    }

    /**
     * The algorithm list is composed of pairs {@link Algorithm} + {@link Problem} which form part of
     * a {@link ExperimentAlgorithm}, which is a decorator for class {@link Algorithm}.
     */
    static List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> configureAlgorithmList(
            List<ExperimentProblem<DoubleSolution>> problemList) {
        List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> algorithms = new ArrayList<>();
        /*
        for (int i = 0; i < problemList.size(); i++) {
            double mutationProbability = 1.0 / problemList.get(i).getProblem().getNumberOfVariables();
            double mutationDistributionIndex = 20.0;
            Algorithm<List<DoubleSolution>> algorithm = new SMPSOBuilder((DoubleProblem) problemList.get(i).getProblem(),
                    new CrowdingDistanceArchive<DoubleSolution>(100))
                    .setMutation(new PolynomialMutation(mutationProbability, mutationDistributionIndex))
                    .setMaxIterations(250)
                    .setSwarmSize(100)
                    .setSolutionListEvaluator(new SequentialSolutionListEvaluator<DoubleSolution>())
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm, problemList.get(i).getTag()));
        }
        */
        for (int i = 0; i < problemList.size(); i++) {
            Algorithm<List<DoubleSolution>> algorithm = new NSGAIIBuilder<DoubleSolution>(
                    problemList.get(i).getProblem(),
                    new SBXCrossover(1.0, 20.0),
                    new PolynomialMutation(1.0 / problemList.get(i).getProblem().getNumberOfVariables(), 20.0))
                    .setMaxEvaluations(10000)
                    .setPopulationSize(100)
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm, problemList.get(i).getTag()));

        }
        for (int i = 0; i < problemList.size(); i++) {
            Problem<DoubleSolution> problem;
            CrossoverOperator<DoubleSolution> crossover;
            MutationOperator<DoubleSolution> mutation;
            SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
            problem = problemList.get(i).getProblem();

            double crossoverProbability = 1.0;
            double crossoverDistributionIndex = 20.0 ;
            crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;

            double mutationProbability = 1.0 / problem.getNumberOfVariables() ;
            double mutationDistributionIndex = 20.0 ;
            mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex) ;

            selection = new BinaryTournamentSelection<DoubleSolution>();
            Algorithm<List<DoubleSolution>> algorithm = new NSGAIIIBuilder<DoubleSolution>(
                    problemList.get(i).getProblem())
                    .setCrossoverOperator(crossover)
                    .setMutationOperator(mutation)
                    .setSelectionOperator(selection)
                    .setMaxIterations(10000)
                    .setPopulationSize(100)
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm,problemList.get(i).getTag()));
        }
        for (int i = 0; i < problemList.size(); i++) {
            Algorithm<List<DoubleSolution>> algorithm = new NSGAVBuilder<DoubleSolution>(
                    problemList.get(i).getProblem(),
                    new SBXCrossover(1.0, 20.0),
                    new PolynomialMutation(1.0 / problemList.get(i).getProblem().getNumberOfVariables(), 20.0))
                    .setMaxEvaluations(10000)
                    .setPopulationSize(100)
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm, "NSGAIV", problemList.get(i).getTag()));

        }

        /*
        for (int i = 0; i < problemList.size(); i++) {
            Algorithm<List<DoubleSolution>> algorithm = new SPEA2Builder<DoubleSolution>(
                    problemList.get(i).getProblem(),
                    new SBXCrossover(1.0, 10.0),
                    new PolynomialMutation(1.0 / problemList.get(i).getProblem().getNumberOfVariables(), 20.0))
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm, problemList.get(i).getTag()));
        }
        */
        return algorithms;
    }
}