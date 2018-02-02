package org.uma.jmetal.experiment;
import org.uma.jmetal.algorithm.multiobjective.nsgav.utilsnsgav.readWriteLatex;

import java.io.IOException;
import java.util.Scanner;

public class DungExperiments {
    public static int enterParameter(int variables){
        Scanner in = new Scanner(System.in);
        System.out.printf("Enter number:  ");
        try
        {
            variables = in.nextInt();
        }
        catch (java.util.InputMismatchException e)
        {
            System.out.println("Invalid Input, the default variables is:"+variables);
        }
        return variables;
    }
    public static void main(String[] args) throws IOException{
        int variables = 10;
        int max_objectives = 10;
	    int min_objectives = 3;
	    int gridPoint = 2;
	    String store = "store";
        //System.out.println("Enter variables (the default value is "+variables+")");
        //variables = enterParameter(variables);
        System.out.println("Enter min objectives (the default value is "+min_objectives+")");
        min_objectives = enterParameter(min_objectives);
        System.out.println("Enter max objectives (the default value is "+max_objectives+")");
        max_objectives = enterParameter(max_objectives);
        System.out.println("Enter gridPoint (the default value is "+gridPoint+")");
        gridPoint = enterParameter(gridPoint);
        System.out.println("Enter Store (the default value is "+store+")");
        store = readWriteLatex.enterStore(store);
        ConstraintProblemsStudy runExperiment = new ConstraintProblemsStudy();
        for (int i = min_objectives; i<=max_objectives;i++){
            runExperiment.ProblemsStudyRun(i+4,i, gridPoint,store);
        }
        readWriteLatex.generateLatex(store);
    }
}
