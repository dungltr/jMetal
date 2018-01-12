package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgav.NSGAVBuilder;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.problem.multiobjective.dtlz.*;
import org.uma.jmetal.problem.multiobjective.zdt.ZDT2;
import org.uma.jmetal.runner.AbstractAlgorithmRunner;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;

import java.io.FileNotFoundException;
import java.util.List;

/**
 * Class to configure and run the NSGA-II algorithm
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class NSGAVRunner extends AbstractAlgorithmRunner {
    /**
     * @param args Command line arguments.
     * @throws JMetalException
     * @throws FileNotFoundException
     * Invoking command:
    java org.uma.jmetal.runner.multiobjective.NSGAIIRunner problemName [referenceFront]
     */
    public static void main(String[] args) throws JMetalException, FileNotFoundException {
        String[] Args = args;
        test(Args,"NSGAV");
        test(Args,"NSGAII");

    }
    public static void test(String[] args, String algorithmName)throws JMetalException, FileNotFoundException{
        Problem<DoubleSolution> problem;
        Algorithm<List<DoubleSolution>> algorithm;
        CrossoverOperator<DoubleSolution> crossover;
        MutationOperator<DoubleSolution> mutation;
        SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
        String referenceParetoFront = "" ;

        String problemName ;
        if (args.length == 1) {
            problemName = args[0];
        } else if (args.length == 2) {
            problemName = args[0] ;
            referenceParetoFront = args[1] ;
        } else {
            //problemName = "org.uma.jmetal.problem.multiobjective.zdt.ZDT1";
            problemName = "org.uma.jmetal.problem.multiobjective.dtlz.dtlz1";
            //referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/ZDT1.pf" ;
            referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/DTLZ1.3D.pf" ;
        }
        problem = new DTLZ1();
        //problem = new ZDT2(); //ProblemUtils.<DoubleSolution> loadProblem(problemName);

        double crossoverProbability = 0.9 ;
        double crossoverDistributionIndex = 20.0 ;
        crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;

        double mutationProbability = 1.0 / problem.getNumberOfVariables() ;
        double mutationDistributionIndex = 20.0 ;
        mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex) ;

        selection = new BinaryTournamentSelection<DoubleSolution>(
                new RankingAndCrowdingDistanceComparator<DoubleSolution>());
        if (algorithmName.toLowerCase().contains("nsgav")){
            algorithm = new NSGAVBuilder<DoubleSolution>(problem, crossover, mutation)
                    .setSelectionOperator(selection)
                    .setMaxEvaluations(25000)
                    .setPopulationSize(300)
                    .build() ;

            List<DoubleSolution> population = algorithm.getResult() ;
            printFinalSolutionSet(population);
            if (!referenceParetoFront.equals("")) {
                printQualityIndicators(population, referenceParetoFront) ;
            }
        }
        if (algorithmName.toLowerCase().contains("nsgaii")){
            algorithm = new NSGAIIBuilder<DoubleSolution>(problem, crossover, mutation)
                    .setSelectionOperator(selection)
                    .setMaxEvaluations(25000)
                    .setPopulationSize(300)
                    .build() ;

            List<DoubleSolution> population = algorithm.getResult() ;
            printFinalSolutionSet(population);
            if (!referenceParetoFront.equals("")) {
                printQualityIndicators(population, referenceParetoFront) ;
            }
        }


    }
}
