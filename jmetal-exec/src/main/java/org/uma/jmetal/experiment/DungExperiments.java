package org.uma.jmetal.experiment;

import java.io.IOException;
import org.uma.jmetal.experiment.ConstraintProblemsStudy;
public class DungExperiments {
    public static void main(String[] args) throws IOException{
    int variables = 10;
    int objectives = 10;
        for (int i = 3; i<objectives;i++){
            ConstraintProblemsStudy.ProblemsStudyRun(variables,i);
        }
    }
}
