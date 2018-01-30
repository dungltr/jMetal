package org.uma.jmetal.experiment;

import java.io.IOException;

public class DungExperiments {
    public static void main (){
        int variables = 10;
        int objectives = 10;
        for (int i = 3; i<objectives;i++){
            try {
                ConstraintProblemsStudy constraintProblemsStudy = new ConstraintProblemsStudy();
                constraintProblemsStudy.set(variables,i);
                constraintProblemsStudy.ProblemsStudyRun();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }
}
