package org.uma.jmetal.experiment;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;

import java.io.IOException;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class DungExperiments {
    private NSGAIIBuilder<DoubleSolution> builder;
    private Problem<DoubleSolution> problem;
    private static final double EPSILON = 0.000000000000001;
    private static final int NUMBER_OF_VARIABLES_OF_THE_MOCKED_PROBLEM = 20;
    private CrossoverOperator<DoubleSolution> crossover;
    private MutationOperator<DoubleSolution> mutation;

    @Before public void startup() {
        problem = mock(Problem.class);
        when(problem.getNumberOfVariables()).thenReturn(NUMBER_OF_VARIABLES_OF_THE_MOCKED_PROBLEM);

        double crossoverProbability = 0.9 ;
        double crossoverDistributionIndex = 20.0 ;
        crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;

        double mutationProbability = 1.0 / problem.getNumberOfVariables() ;
        double mutationDistributionIndex = 20.0 ;
        mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex) ;

        builder = new NSGAIIBuilder<DoubleSolution>(problem, crossover, mutation);
    }

    @After public void cleanup() {
        problem = null;
        builder = null;
    }
    @Test public static void main(String[] args) throws IOException{
    int variables = 10;
    int objectives = 10;
        for (int i = 3; i<objectives;i++){
            ConstraintProblemsStudy.ProblemsStudyRun(variables,i);
        }
    }
}
