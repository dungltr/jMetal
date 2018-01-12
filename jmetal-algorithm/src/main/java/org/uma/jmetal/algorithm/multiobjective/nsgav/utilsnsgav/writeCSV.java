package org.uma.jmetal.algorithm.multiobjective.nsgav.utilsnsgav;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public class writeCSV {
    private static final String COMMA_DELIMITER = ",";
    private static final String NEW_LINE_SEPARATOR = "\n";
    public void addArray2Csv(String filename, double[] tmp) throws IOException {
        Path filePath = Paths.get(filename);
        if (!Files.exists(filePath)) {
            Files.createFile(filePath);
        }
        String add = "";
        int i = 0;
        for (i = 0; tmp.length - 1 > i; i++)
            add = add + tmp[i] + COMMA_DELIMITER;
        if (tmp.length - 1 == i)
            add = add + tmp[i] + NEW_LINE_SEPARATOR;
        Files.write(filePath, add.getBytes(), StandardOpenOption.APPEND);
    }
}
