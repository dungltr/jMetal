package org.uma.jmetal.algorithm.multiobjective.nsgav.utilsPopulation;
//import org.moeaframework.core.NondominatedPopulation;
//import org.moeaframework.core.Population;
//import org.moeaframework.core.Solution;
//import org.moeaframework.core.variable.EncodingUtils;
//import org.moeaframework.util.Vector;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.EnvironmentalSelection;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.ReferencePoint;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.SolutionUtils;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DominanceRanking;

import java.util.ArrayList;
import java.util.List;

public class utilsPopulation <S extends Solution<?>>{
    public void printPopulation(List<S> result) {
        System.out.println("Num of Solutions: "+ result.size());
        //double[][] matrixResult = new double [result.size()][result.get(0).getNumberOfObjectives()];
        // 2.2.4 Read solutions
        for (int m = 0; m < result.size(); m++) {
            Solution solution = result.get(m);
            double[] objectives = new double[solution.getNumberOfObjectives()];
            for(int i = 0; i< solution.getNumberOfObjectives(); i++)
                objectives[i] = -solution.getObjective(i);
            //Negate objectives to return them to their maximized form.
            //objectives = Vector.negate(objectives);//.negate(objectives);
            //2.2.5 Print results
            System.out.println("\n    Solution " + (m + 1) + ":");
            for (int i=0; i < objectives.length; i++)
                System.out.print("      Obj "+i+": " + -objectives[i]);
            //System.out.println("    Con 1: " + solution.getConstraint(0));

/*            for(int j=0;j<x.length;j++){
                System.out.print("      Var " + (j+1) + ":" + x[j]+"\n");
            }
*/
            //for (int j=0; j < objectives.length ;j++)
            //  matrixResult[m][j] = -objectives[j];
        }
    }
    public void printSolution(S solution) {
        System.out.println("Solutions: ");
        //       double[][] matrixResult = new double [result.size()][result.get(0).getNumberOfObjectives()];
        double[] objectives = new double[solution.getNumberOfObjectives()];
        for(int i = 0; i< solution.getNumberOfObjectives(); i++)
            objectives[i] = -solution.getObjective(i);

        //Negate objectives to return them to their maximized form.
        //objectives = Vector.negate(objectives);//.negate(objectives);
        //2.2.5 Print results
//            System.out.println("\n    Solution " + (m + 1) + ":");
        for (int i=0; i < objectives.length; i++)
            System.out.print("      Obj "+i+": " + -objectives[i]);
        //System.out.println("    Con 1: " + solution.getConstraint(0));

        //for(int j=0;j<x.length;j++){
        //    System.out.print("      Var " + (j+1) + ":" + x[j]+"\n");
        //}

//            for (int j=0; j < objectives.length ;j++)
//                matrixResult[m][j] = -objectives[j];
    }
    public double[] objectives (S solution){
        double[] objectives = new double[solution.getNumberOfObjectives()];
        for(int i = 0; i< solution.getNumberOfObjectives(); i++)
            objectives[i] = solution.getObjective(i);
        return objectives;
    }
    public static void printMatrix(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++)
            {
                System.out.printf("%.6f",a[i][j]);
                System.out.print(" ");
            }
            System.out.println();
        }
    }
    public static void printArray(double[] a) {
        int m = a.length;
        for (int i = 0; i < m; i++){
            System.out.printf("%.6f",a[i]);
            System.out.print(" ");
        }
    }

}

