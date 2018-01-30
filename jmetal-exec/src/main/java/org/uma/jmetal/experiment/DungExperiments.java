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
        int objectives = 10;
        System.out.println("Enter variables");
        variables = enterParameter(variables);
        System.out.println("Enter objectives");
        objectives = enterParameter(objectives);
        for (int i = 3; i<objectives;i++){
            ConstraintProblemsStudy.ProblemsStudyRun(variables,i);
        }
    }
}
