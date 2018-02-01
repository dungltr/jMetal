package org.uma.jmetal.algorithm.multiobjective.nsgav;

import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.AlgorithmBuilder;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.util.List;

//import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIMeasures;
//import org.uma.jmetal.algorithm.multiobjective.nsgaii.SteadyStateNSGAII;
//import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder.NSGAIIVariant;
//import org.uma.jmetal.algorithm.multiobjective.nsgav.NSGAV;


/** Builder class */
public class NSGAVBuilder<S extends Solution<?>> implements AlgorithmBuilder<NSGAV<S>>{
	public enum NSGAVVariant {NSGAV}

	  /**
	   * NSGAIIBuilder class
	   */
	  private final Problem<S> problem;
	  private int maxEvaluations;
	  private int populationSize;
	  private int gridPoint;
	  private CrossoverOperator<S>  crossoverOperator;
	  private MutationOperator<S> mutationOperator;
	  private SelectionOperator<List<S>, S> selectionOperator;
	  private SolutionListEvaluator<S> evaluator;

	  private NSGAVVariant variant;

  /** Builder constructor */
public NSGAVBuilder(Problem<S> problem, CrossoverOperator<S> crossoverOperator,
		      MutationOperator<S> mutationOperator) {
		    this.problem = problem;
		    maxEvaluations = 25000;
		    populationSize = 100;
		    gridPoint = 2;
		    this.crossoverOperator = crossoverOperator ;
		    this.mutationOperator = mutationOperator ;
		    selectionOperator = new BinaryTournamentSelection<S>(new RankingAndCrowdingDistanceComparator<S>()) ;
		    evaluator = new SequentialSolutionListEvaluator<S>();

		    this.variant = NSGAVVariant.NSGAV ;
		  }

  public NSGAVBuilder<S> setMaxEvaluations(int maxEvaluations) {
    if (maxEvaluations < 0) {
      throw new JMetalException("maxEvaluations is negative: " + maxEvaluations);
    }
    this.maxEvaluations = maxEvaluations;

    return this;
  }
    public NSGAVBuilder<S> setGridPoint(int gridPoint) {
        if (gridPoint < 0) {
            throw new JMetalException("gridPoint is negative: " + gridPoint);
        }
        this.gridPoint = gridPoint;
        return this;
    }
  public NSGAVBuilder<S> setPopulationSize(int populationSize) {
    if (populationSize < 0) {
      throw new JMetalException("Population size is negative: " + populationSize);
    }

    this.populationSize = populationSize;

    return this;
  }

  public NSGAVBuilder<S> setSelectionOperator(SelectionOperator<List<S>, S> selectionOperator) {
    if (selectionOperator == null) {
      throw new JMetalException("selectionOperator is null");
    }
    this.selectionOperator = selectionOperator;

    return this;
  }

  public NSGAVBuilder<S> setSolutionListEvaluator(SolutionListEvaluator<S> evaluator) {
    if (evaluator == null) {
      throw new JMetalException("evaluator is null");
    }
    this.evaluator = evaluator;

    return this;
  }


  public NSGAVBuilder<S> setVariant(NSGAVVariant variant) {
    this.variant = variant;

    return this;
  }

  public NSGAV<S> build() {
    NSGAV<S> algorithm = null ;
    if (variant.equals(NSGAVVariant.NSGAV)) {
      algorithm = new NSGAV<S>(problem, maxEvaluations, populationSize, crossoverOperator,
          mutationOperator, selectionOperator, evaluator, gridPoint);
    } 
    return algorithm ;
  }

  /* Getters */
  public Problem<S> getProblem() {
    return problem;
  }

  public int getMaxIterations() {
    return maxEvaluations;
  }

  public int getGridPoint() {
        return gridPoint;
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

  public SolutionListEvaluator<S> getSolutionListEvaluator() {
    return evaluator;
  }
}