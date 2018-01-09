package org.uma.jmetal.algorithm.multiobjective.nsgav;

import org.uma.jmetal.algorithm.multiobjective.nsgav.NSGAV;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.AlgorithmBuilder;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.util.List;


/** Builder class */
public class NSGAVBuilder<S extends Solution<?>> implements AlgorithmBuilder<NSGAV<S>>{
  // no access modifier means access from classes within the same package
  private Problem<S> problem ;
  private int maxIterations ;
  private int populationSize ;
  private CrossoverOperator<S> crossoverOperator ;
  private MutationOperator<S> mutationOperator ;
  private SelectionOperator<List<S>, S> selectionOperator ;

  private SolutionListEvaluator<S> evaluator ;

  /** Builder constructor */
  public NSGAVBuilder(Problem<S> problem) {
    this.problem = problem ;
    maxIterations = 250 ;
    populationSize = 100 ;
    evaluator = new SequentialSolutionListEvaluator<S>() ;
  }

  public NSGAVBuilder<S> setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations ;

    return this ;
  }

  public NSGAVBuilder<S> setPopulationSize(int populationSize) {
    this.populationSize = populationSize ;

    return this ;
  }

  public NSGAVBuilder<S> setCrossoverOperator(CrossoverOperator<S> crossoverOperator) {
    this.crossoverOperator = crossoverOperator ;

    return this ;
  }

  public NSGAVBuilder<S> setMutationOperator(MutationOperator<S> mutationOperator) {
    this.mutationOperator = mutationOperator ;

    return this ;
  }

  public NSGAVBuilder<S> setSelectionOperator(SelectionOperator<List<S>, S> selectionOperator) {
    this.selectionOperator = selectionOperator ;

    return this ;
  }

  public NSGAVBuilder<S> setSolutionListEvaluator(SolutionListEvaluator<S> evaluator) {
    this.evaluator = evaluator ;

    return this ;
  }

  public SolutionListEvaluator<S> getEvaluator() {
    return evaluator;
  }

  public Problem<S> getProblem() {
    return problem;
  }

  public int getMaxIterations() {
    return maxIterations;
  }

  public int getPopulationSize() {
    return populationSize;
  }

  public CrossoverOperator<S> getCrossoverOperator() {
    return crossoverOperator;
  }

  public MutationOperator<S> getMutationOperator() {
    return mutationOperator;
  }

  public SelectionOperator<List<S>, S> getSelectionOperator() {
    return selectionOperator;
  }

  public NSGAV<S> build() {
    return new NSGAV<>(this) ;
  }
}
