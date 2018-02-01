package org.uma.jmetal.experiment;
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
        //System.out.println("Enter variables (the default value is "+variables+")");
        //variables = enterParameter(variables);
        System.out.println("Enter min objectives (the default value is "+min_objectives+")");
        min_objectives = enterParameter(min_objectives);
        System.out.println("Enter max objectives (the default value is "+max_objectives+")");
        max_objectives = enterParameter(max_objectives);
        for (int i = min_objectives; i<=max_objectives;i++){
            ConstraintProblemsStudy.ProblemsStudyRun(i+4,i);
        }
    }
}
