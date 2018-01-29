package org.uma.jmetal.experiment;
import org.uma.jmetal.algorithm.multiobjective.nsgav.utilsnsgav.*;
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
    private static final int INDEPENDENT_RUNS = 10;

    public static void main(String[] args) throws IOException {
        /*if (args.length != 1) {
            throw new JMetalException("Missing argument: experimentBaseDirectory");
        }*/
        int variables = 10;
        int objecitves = 4;
        String experimentBaseDirectory = ReadFile.readhome("HOME_jMetal");//args[0];

        List<ExperimentProblem<DoubleSolution>> problemList = new ArrayList<>();
        //problemList.add(new ExperimentProblem<>(new DTLZ1(variables,objecitves-2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ2(variables,objecitves-2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ3(variables,objecitves-2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ4(variables,objecitves-2)));
        problemList.add(new ExperimentProblem<>(new DTLZ1()));
        problemList.add(new ExperimentProblem<>(new DTLZ2()));
        problemList.add(new ExperimentProblem<>(new DTLZ3()));
        problemList.add(new ExperimentProblem<>(new DTLZ4()));
        //problemList.add(new ExperimentProblem<>(new DTLZ1(variables,objecitves+2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ2(variables,objecitves+2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ3(variables,objecitves+2)));
        //problemList.add(new ExperimentProblem<>(new DTLZ4(variables,objecitves+2)));
        /*problemList.add(new ExperimentProblem<>(new DTLZ5()));
        problemList.add(new ExperimentProblem<>(new DTLZ6()));
        problemList.add(new ExperimentProblem<>(new DTLZ7()));
        */
        //problemList.add(new ExperimentProblem<>(new UF1()));
        //problemList.add(new ExperimentProblem<>(new UF2()));
        //problemList.add(new ExperimentProblem<>(new UF3()));
        //problemList.add(new ExperimentProblem<>(new UF4()));
        //problemList.add(new ExperimentProblem<>(new UF5()));
        //problemList.add(new ExperimentProblem<>(new UF6()));
        //problemList.add(new ExperimentProblem<>(new UF7()));
        //problemList.add(new ExperimentProblem<>(new UF8()));
        //problemList.add(new ExperimentProblem<>(new UF9()));
        //problemList.add(new ExperimentProblem<>(new UF10()));

        //problemList.add(new ExperimentProblem<>(new ZDT1()));
        //problemList.add(new ExperimentProblem<>(new ZDT2()));
        //problemList.add(new ExperimentProblem<>(new ZDT3()));
        //problemList.add(new ExperimentProblem<>(new ZDT4()));

        List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> algorithmList =
                configureAlgorithmList(problemList);
        /*
        List<String> referenceFrontFileNames =
                Arrays.asList("DTLZ1.3D.pf", "DTLZ2.3D.pf", "DTLZ3.3D.pf","DTLZ4.3D.pf"
                        //"UF1.pf", "UF2.pf", "UF3.pf","UF4.pf"// "UF5.pf"//, "UF6.pf"//, "UF7.pf", "UF8.pf", "UF9.pf", "UF10.pf"
                        //"DTLZ1.2D.pf", "DTLZ2.2D.pf", "DTLZ3.2D.pf","DTLZ4.2D.pf"//, "DTLZ5.2D.pf"//, "DTLZ3.2D.pf"
                        //"ZDT1.pf","ZDT2.pf","ZDT3.pf","ZDT4.pf"
                        //"UF1.pf", "UF2.pf", "UF3.pf"
                        // "DTLZ4.3D.pf", "DTLZ5.3D.pf","DTLZ6.3D.pf"
                        // , "DTLZ7.3D.pf");

                );//);
        */
        String experimentName = "3AlgorithmsDTLZ3Objective_Pop100Max10000";
        String homeFile = ReadFile.readhome("HOME_jMetal")+"/"+experimentName;

        Experiment<DoubleSolution, List<DoubleSolution>> experiment =
                new ExperimentBuilder<DoubleSolution, List<DoubleSolution>>(experimentName)
                        .setAlgorithmList(algorithmList)
                        .setProblemList(problemList)
                        .setReferenceFrontDirectory("/pareto_fronts")
                        //.setReferenceFrontFileNames(referenceFrontFileNames)
                        .setExperimentBaseDirectory(experimentBaseDirectory)
                        .setOutputParetoFrontFileName("FUN")
                        .setOutputParetoSetFileName("VAR")
                        .setIndicatorList(Arrays.asList(
                                //new Epsilon<DoubleSolution>(),
                                //new Spread<DoubleSolution>(),
                                //new GenerationalDistance<DoubleSolution>()))//,
                                //new PISAHypervolume<DoubleSolution>(),
                                new InvertedGenerationalDistance<DoubleSolution>()))//,
                                //new InvertedGenerationalDistancePlus<DoubleSolution>()))//
                        .setIndependentRuns(INDEPENDENT_RUNS)
                        .setNumberOfCores(8)
                        .build();

        new ExecuteAlgorithms<>(experiment).run();
        new ComputeQualityIndicators<>(experiment).run();
        new GenerateLatexTablesWithStatistics(experiment).run();
        new GenerateWilcoxonTestTablesWithR<>(experiment).run();
        new GenerateFriedmanTestTables<>(experiment).run();
        new GenerateBoxplotsWithR<>(experiment).setRows(3).setColumns(3).setDisplayNotch().run();

        String Caption = "Inverted Generational Distance";
        //Caption = "";

        List<String> Problems = new ArrayList<>();
        for (int i = 0; i< problemList.size(); i++){
            Problems.add(problemList.get(i).getTag());
            //System.out.println(problemList.get(i).getTag());
        }

        String [] algorithms = new String [algorithmList.size()];
        for (int i = 0; i< algorithmList.size(); i++){
            algorithms[i] = algorithmList.get(i).getAlgorithmTag();
            //System.out.println(algorithmList.get(i).getAlgorithmTag());
        }
        GeneratorLatexTable.GeneratorComputeTimeToLatex(homeFile, Caption, Problems, algorithms);
    }

    /**
     * The algorithm list is composed of pairs {@link Algorithm} + {@link Problem} which form part of
     * a {@link ExperimentAlgorithm}, which is a decorator for class {@link Algorithm}.
     */
    static List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> configureAlgorithmList(
            List<ExperimentProblem<DoubleSolution>> problemList) {
        List<ExperimentAlgorithm<DoubleSolution, List<DoubleSolution>>> algorithms = new ArrayList<>();
        double crossoverProbability = 1.0;
        double crossoverDistributionIndex = 5.0 ;

        double mutationProbability = 1.0;
        double mutationDistributionIndex = 10.0 ;
        int MaxEvaluations = 10000;
        int PopulationSize = 100;
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
                    new SBXCrossover(crossoverProbability, crossoverDistributionIndex),
                    new PolynomialMutation(mutationProbability / problemList.get(i).getProblem().getNumberOfVariables(), mutationDistributionIndex))
                    .setMaxEvaluations(MaxEvaluations)
                    .setPopulationSize(PopulationSize)
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm, problemList.get(i).getTag()));

        }

        for (int i = 0; i < problemList.size(); i++) {
            Problem<DoubleSolution> problem;
            CrossoverOperator<DoubleSolution> crossover;
            MutationOperator<DoubleSolution> mutation;
            SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
            problem = problemList.get(i).getProblem();


            crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;


            mutation = new PolynomialMutation(mutationProbability / problemList.get(i).getProblem().getNumberOfVariables(), mutationDistributionIndex) ;

            selection = new BinaryTournamentSelection<DoubleSolution>();
            Algorithm<List<DoubleSolution>> algorithm = new NSGAIIIBuilder<DoubleSolution>(
                    problemList.get(i).getProblem())
                    .setCrossoverOperator(crossover)
                    .setMutationOperator(mutation)
                    .setSelectionOperator(selection)
                    .setMaxIterations(MaxEvaluations)
                    .setPopulationSize(PopulationSize)
                    .build();
            algorithms.add(new ExperimentAlgorithm<>(algorithm,problemList.get(i).getTag()));
        }

        for (int i = 0; i < problemList.size(); i++) {
            Algorithm<List<DoubleSolution>> algorithm = new NSGAVBuilder<DoubleSolution>(
                    problemList.get(i).getProblem(),
                    new SBXCrossover(crossoverProbability, crossoverDistributionIndex),
                    new PolynomialMutation(mutationProbability / problemList.get(i).getProblem().getNumberOfVariables(), mutationDistributionIndex))
                    .setMaxEvaluations(MaxEvaluations)
                    .setPopulationSize(PopulationSize)
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
