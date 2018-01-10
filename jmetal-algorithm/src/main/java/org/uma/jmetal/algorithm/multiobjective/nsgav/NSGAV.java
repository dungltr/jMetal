package org.uma.jmetal.algorithm.multiobjective.nsgav;

//import org.moeaframework.core.Population;
import org.uma.jmetal.algorithm.impl.AbstractGeneticAlgorithm;
import org.uma.jmetal.algorithm.multiobjective.nsgav.NSGAVBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgav.utilsPopulation.utilsPopulation;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.EnvironmentalSelection;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.ReferencePoint;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DominanceRanking;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

/**
 * Created by ajnebro on 30/10/14.
 * Modified by Juanjo on 13/11/14
 *
 * This implementation is based on the code of Tsung-Che Chiang
 * http://web.ntnu.edu.tw/~tcchiang/publications/nsga3cpp/nsga3cpp.htm
 */
@SuppressWarnings("serial")
public class NSGAV<S extends Solution<?>> extends AbstractGeneticAlgorithm<S, List<S>> {
  protected final int maxEvaluations;

  protected final SolutionListEvaluator<S> evaluator;

  protected int evaluations;

  /**
   * Constructor
   */
  public NSGAV(Problem<S> problem, int maxEvaluations, int populationSize,
      CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, SolutionListEvaluator<S> evaluator) {
    super(problem);
    this.maxEvaluations = maxEvaluations;
    setMaxPopulationSize(populationSize); ;

    this.crossoverOperator = crossoverOperator;
    this.mutationOperator = mutationOperator;
    this.selectionOperator = selectionOperator;

    this.evaluator = evaluator;
  }

    @Override
    protected void initProgress() {
        evaluations = getMaxPopulationSize();
      }

      @Override protected void updateProgress() {
        evaluations += getMaxPopulationSize() ;
      }

      @Override protected boolean isStoppingConditionReached() {
        return evaluations >= maxEvaluations;
      }

      @Override protected List<S> evaluatePopulation(List<S> population) {
        population = evaluator.evaluate(population, getProblem());

        return population;
      }

    @Override
    protected List<S> selection(List<S> population) {
        List<S> matingPopulation = new ArrayList<>(population.size()) ;
        for (int i = 0; i < getMaxPopulationSize(); i++) {
            S solution = selectionOperator.execute(population);
            matingPopulation.add(solution) ;
        }

        return matingPopulation;
    }

    @Override
    protected List<S> reproduction(List<S> population) {
        List<S> offspringPopulation = new ArrayList<>(getMaxPopulationSize());
        for (int i = 0; i < getMaxPopulationSize(); i+=2) {
            List<S> parents = new ArrayList<>(2);
            parents.add(population.get(i));
            parents.add(population.get(Math.min(i + 1, getMaxPopulationSize()-1)));

            List<S> offspring = crossoverOperator.execute(parents);

            mutationOperator.execute(offspring.get(0));
            mutationOperator.execute(offspring.get(1));

            offspringPopulation.add(offspring.get(0));
            offspringPopulation.add(offspring.get(1));
        }
        return offspringPopulation ;
    }
    @Override
    protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
		System.out.println("The begining of replacement function******************************************************");
        List<S> jointPopulation = new ArrayList<>();
        jointPopulation.addAll(population) ;
        jointPopulation.addAll(offspringPopulation) ;
        Ranking<S> ranking = computeRanking(jointPopulation);
        //List<Solution> pop = crowdingDistanceSelection(ranking);
        List<S> pop = new ArrayList<>();
        List<List<S>> fronts = new ArrayList<>();
        int rankingIndex = 0;
        int candidateSolutions = 0;
        /*------------------------------------------------------------------------------------------*/
        while (candidateSolutions <= getMaxPopulationSize()) {
        	//System.out.println("\n The ranks of Fronts is:="+fronts.size());
			//System.out.println("\n The size of pop of Fronts is:="+candidateSolutions);
            fronts.add(ranking.getSubfront(rankingIndex));// update fronts
            candidateSolutions += ranking.getSubfront(rankingIndex).size();// update number of fronts add to population
            if ((pop.size() + ranking.getSubfront(rankingIndex).size()) <= getMaxPopulationSize())
                addRankedSolutionsToPopulation(ranking, rankingIndex, pop);
            rankingIndex++;
        }
		//System.out.println("The ranks of Fronts after while () is:="+fronts.size());
		System.out.println("The size of candidates after while is:="+candidateSolutions);
		System.out.println("The size of pop after while is:="+pop.size());
        // A copy of the reference list should be used as parameter of the environmental selection
        //EnvironmentalSelection<S> selection =
        //        new EnvironmentalSelection<>(fronts,getMaxPopulationSize(),getReferencePointsCopy(),
        //                getProblem().getNumberOfObjectives());
        /* When candidateSolution > MaxPopulation */
        //pop = selection.execute(pop);
        /*------------------------------------------------------------------------------------------*/
        if(pop.size()<getMaxPopulationSize()){
			List<S> currentFront = new ArrayList<>();
			List<S> previousFront = new ArrayList<>();
			currentFront = ranking.getSubfront(rankingIndex-1);
			previousFront = ranking.getSubfront(rankingIndex-2);
			//System.out.println("The size of previousFront is:="+previousFront.size());
			//System.out.println("The size of currentFront is:="+currentFront.size());
			//pop.clear();
			//for(int i=0; i<rankingIndex-1;i++) {
			//	pop.addAll(ranking.getSubfront(i));
			//}
			//pop.removeAll(currentFront);
			//System.out.println("The size of population before filting is:="+pop.size());
			int addMore = getMaxPopulationSize() - pop.size();//candidateSolutions + ranking.getSubfront(rankingIndex).size();
			//System.out.println("The size of results is required:="+addMore);
			List<S> resultFilter =
					filter(previousFront, currentFront, addMore);
			//System.out.println("The size of resultFilter is:="+resultFilter.size());
			for (int i=0;i<resultFilter.size();i++) {
				pop.add(resultFilter.get(i));
			}

			/*for (int i=0;i<currentFront.size();i++) {
				pop.add(currentFront.get(i));
			}
			*/
		}else {
        	if (pop.size()==getMaxPopulationSize())
				System.out.println("The size of pop is Max and no need to run truncate");
        	else System.out.println("The size of pop bigger than Max and should be checked error");
		}
		System.out.println("The size of population after filting is:="+pop.size());
		System.out.println("The end of replacement function *********************************************************");
        return pop;
    }
    protected List<double []> updateDeltas(List<S> previousFront, double epsilon){
		List<double []> Deltas = new ArrayList<double[]>();
		for (int i = 0; i< previousFront.size(); i++) {
			double [] temp = new double [previousFront.get(i).getNumberOfObjectives()];
			for (int j = 0; j< previousFront.get(i).getNumberOfObjectives(); j++) {
				temp[j] = previousFront.get(i).getObjective(j)*epsilon;
			}
			Deltas.add(temp);
			//System.out.println("\n This is the solution:");
			//NSGAIV.utilsPopulation.printSolution(solution);
			//System.out.println("\n This is the delta");
			//utilsPopulation.printArray(temp);
		}
		return Deltas;
	}
    protected S findMaxSolution (List<S> resultFilter){
		double[] Distance = new double[resultFilter.size()];
		for (int i=0; i<resultFilter.size();i++) {
			Distance[i] = 0;
			for (int j = 0; j < resultFilter.get(i).getNumberOfObjectives(); j++) {
				Distance[i] += resultFilter.get(i).getObjective(j)*resultFilter.get(i).getObjective(j);
			}
		}
		double max = 0;
		int index = 0;
		for (int i=0; i<Distance.length;i++){
			max = Math.max(max, Distance[i]);
			if (max==Distance[i]) index = i;
		}		
		return resultFilter.get(index);
	}
    protected S updateObjecitves(S solution, double[] epsilon){
    	S tempSolution = solution;
		double [] temp = new double[solution.getNumberOfObjectives()];
		for (int i = 0; i< solution.getNumberOfObjectives(); i++){
			temp[i] = epsilon[i];
		}
		for (int i = 0; i< solution.getNumberOfObjectives(); i++){
			temp[i] = temp[i] + solution.getObjective(i);
			tempSolution.setObjective(i, solution.getObjective(i) + epsilon[i]);
		}
		return tempSolution;
	}
    protected List<S> filter (List<S> previousFront, List<S> currentFront, int newSize) {//,Comparator<? super Solution> comparator) {
    		List<S> resultFilter = new ArrayList<>();	
		double epsilon = 0.001;
		int k=0;
		int size=0;
		List<S> temp = new ArrayList<>();
		for (S solution: previousFront) {
			temp.add(solution);
		}
		//System.out.println("\nThis is the Store");
		List<double []> Store = updateDeltas (temp, 1);
		//System.out.println("\nThis is the Deltas");
		List<double []> Deltas = updateDeltas (temp, epsilon);
		System.out.println("The newSize is"+newSize);
		resultFilter.clear();
		while (resultFilter.size()<newSize) {
			resultFilter.clear();
			//System.out.println("\nThis is the previousFront");
			//utilsPopulation UtilsPopulation = new utilsPopulation();
			//UtilsPopulation.printPopulation(previousFront);
			//System.out.println("\nThis is the tempFront");
			//UtilsPopulation.printPopulation(temp);
			k++;
			int[][] dominanceChecks = new int[currentFront.size()][previousFront.size()];
			for (int i = 0; i < currentFront.size(); i++) {
				S si = currentFront.get(i);				
				for (int j = 0; j < previousFront.size(); j++) {						 
						S sj = updateObjecitves(temp.get(j),Deltas.get(j));//The previous Front value will be changed after this
						//temp.set(j, sj);
						//sj.setObjectives(updateObjecitves(temp.get(j). getObjectives(),Deltas.get(j)));
						dominanceChecks[i][j] = compare(si, sj);//compareSolutionAproximate(si, sj, k, epsilon);//.compare(si, sj);
				}
			}

			for (int i = 0; i < currentFront.size(); i++) {
				List<Integer> dominates = new ArrayList<Integer>();
				int dominatedCount = 0;	
				int dominatesCount = 0;	
				for (int j = 0; j < previousFront.size(); j++) {
						if (dominanceChecks[i][j] < 0) {
							dominates.add(j);
							dominatesCount += 1;
						} else {
							dominatedCount += 1;
						}
				}			
				if (dominatedCount == 0) {
					S solution = currentFront.get(i);
					resultFilter.add(solution);
				}		
			}
			if (resultFilter.size()>=newSize) {
				if (resultFilter.size()==newSize){
					//System.out.println("The size of results interrupted at k = "+k+"and Size:="+resultFilter.size()+"and newSize is:="+newSize);
					for (int i = 0; i<previousFront.size();i++){
						for (int j = 0; j<previousFront.get(i).getNumberOfObjectives();j++) {
							temp.get(i).setObjective(j, Store.get(i)[j]);//;setObjectives(Store.get(m));
						}
					}
					//currentFront.clear();
					//currentFront.addAll(resultFilter);
//					return resultFilter;
				}else{
					//System.out.println("**************************Need to reduce Deltas at k = "+k+"and Size:="+resultFilter.size()+"and newSize is:="+newSize);
					while(resultFilter.size()>newSize)
					resultFilter.remove(findMaxSolution (resultFilter));
					//System.out.println("After reduce Deltas at k = "+k+"and Size:="+resultFilter.size()+"and newSize is:="+newSize);
					for (int i = 0; i<previousFront.size();i++){
						for (int j = 0; j<previousFront.get(i).getNumberOfObjectives();j++) {
							temp.get(i).setObjective(j, Store.get(i)[j]);//;setObjectives(Store.get(m));
						}
					}
//					currentFront.clear();
//					currentFront.addAll(resultFilter);
//					return resultFilter;
				}
			}
			size = resultFilter.size();
		}
		System.out.println("After truncated at k = "+k+"and Size:="+resultFilter.size()+"and newSize is:="+newSize);
		return resultFilter;
	}
    protected int compare(S solution1, S solution2) {
		boolean dominate1 = false;
		boolean dominate2 = false;

		for (int i = 0; i < solution1.getNumberOfObjectives(); i++) {
			if (solution1.getObjective(i) < solution2.getObjective(i)) {
				dominate1 = true;

				if (dominate2) {
					return 0;
				}
			} else if (solution1.getObjective(i) > solution2.getObjective(i)) {
				dominate2 = true;

				if (dominate1) {
					return 0;
				}
			}
		}

		if (dominate1 == dominate2) {
			return 0;
		} else if (dominate1) {
			return -1;
		} else {
			return 1;
		}
	}
    @Override
    public List<S> getResult() {
        return getNonDominatedSolutions(getPopulation()) ;
    }

    protected Ranking<S> computeRanking(List<S> solutionList) {
        Ranking<S> ranking = new DominanceRanking<>() ;
        ranking.computeRanking(solutionList) ;

        return ranking ;
    }

    protected void addRankedSolutionsToPopulation(Ranking<S> ranking, int rank, List<S> population) {
        List<S> front ;

        front = ranking.getSubfront(rank);

        for (int i = 0 ; i < front.size(); i++) {
            population.add(front.get(i));
        }
    }

    protected List<S> getNonDominatedSolutions(List<S> solutionList) {
        return SolutionListUtils.getNondominatedSolutions(solutionList) ;
    }

    @Override public String getName() {
        return "NSGAV" ;
    }

    @Override public String getDescription() {
        return "Nondominated Sorting Genetic Algorithm version V" ;
    }

}
