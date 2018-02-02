package org.uma.jmetal.algorithm.multiobjective.nsgav.utilsnsgav;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;

public class readWriteLatex {
    private static final int INDEPENDENT_RUNS = 50;
    public int PopulationSize = 500;
    public int MaxEvaluations = 10000;
    private int variable;
    private int objecitve;
    private int gridPoint;
    //String experimentName = "DTLZ"+variables+"V"+objecitves+"O"+this.gridPoint+"G"+this.PopulationSize+"P"+this.MaxEvaluations+"E"+INDEPENDENT_RUNS+"R";
    public static String referenceFronts = "referenceFronts";
    public static String homeFile = ReadFile.readhome("HOME_jMetal")+"/store2";//+"/"+experimentName;

    public static String[] listFolder(String path){
        File[] directories = new File(path).listFiles(File::isDirectory);
        String [] nameFolders = new String[directories.length];
        for (int i = 0; i< directories.length;i++){
            nameFolders[i] = directories[i].toString();
            nameFolders[i] = nameFolders[i].replace(path,"").replace("/","");

        }
        return nameFolders;
    }
    public static String[] listTexFile(String[] nameFolders){
        String[] nameFiles = listFolder(homeFile);
        for (int i = 0; i< nameFiles.length; i++){
            nameFiles[i] = homeFile + "/"+ nameFiles[i] + "/latex/" + nameFiles[i] +".tex";
            //System.out.println(nameFiles[i]);
        }
        return nameFiles;
    }
    public static String[] listComputeFile(String[] nameFolders){
        String[] nameFiles = listFolder(homeFile);
        for (int i = 0; i< nameFiles.length; i++){
            nameFiles[i] = homeFile + "/"+ nameFiles[i] + "/latex/" + "ComputeTime.tex";
            //System.out.println(nameFiles[i]);
        }
        return nameFiles;
    }
    /*
    public static List<String> generateTableHeader(String[] file) throws IOException {
        List<String> headers = new ArrayList<>();
        for (int i = 0; i < 8; i++){
            headers.add(Files.readAllLines(Paths.get(file[0])).get(i+11));
            //System.out.println(Files.readAllLines(Paths.get(file[0])).get(i+11));
        }
        return headers;
    }
    */
    public static List<String> generateTableHeader(String[] file, int index) throws IOException {
        List<String> headers = new ArrayList<>();
        for (int i = 0; i < 8; i++){
            headers.add(Files.readAllLines(Paths.get(file[0])).get(i+index));
            //System.out.println(Files.readAllLines(Paths.get(file[0])).get(i+11));
        }
        return headers;
    }
    public static List<String> generateTableBottom() throws IOException {
        List<String> bottom = new ArrayList<>();
        bottom.add("\\hline");
        bottom.add("\\end{tabular}\n\\end{scriptsize}\n\\end{table}\n");
        return bottom;
    }
    public static List<String> generateTableBody(String[] file, int line, int headerIndex) throws IOException {
        List<String> body = new ArrayList<>();
        for (int i = 0; i < file.length; i++){
            body.add(Files.readAllLines(Paths.get(file[i])).get(line+headerIndex+8));
            //System.out.println(Files.readAllLines(Paths.get(file[i])).get(headerIndex+8));
        }
        return body;
    }
    public static List<String> sortString(List<String> body){
        List<String> temp = new ArrayList<>();
        for (int i = 3; i < body.size() + 3; i++){
            for (int j = 0; j<body.size();j++){
                if(body.get(j).contains("& "+i+" &")){
                    temp.add(body.get(j));
                }
            }
            //System.out.println(temp.get(temp.size()-1));
        }
        return temp;
    }
    public static List<String> addListString(List<String> source, List<String> destination){
        List<String> temp = source;
        for (int j = 0; j< destination.size(); j++){
            temp.add(destination.get(j));
        }
        return temp;
    }
    public static List<String> fileTex(String[] listTexFile, int[] index, int headerIndex) throws IOException{
        List<String> headers = generateTableHeader(listTexFile,headerIndex);
        List<String> body = new ArrayList<>();
        for (int i = 0; i < index.length; i ++){
            List<String> temp = sortString(generateTableBody(listTexFile,index[i],headerIndex));
            for (int j = 0; j< temp.size(); j++){
                body.add(temp.get(j));
            }
        }
        List<String> bottom = generateTableBottom();
        headers = addListString(headers,body);
        headers = addListString(headers,bottom);
        printStringList(headers);
        return headers;
    }
    public static void printStringList(List<String> list){
        for (int i = 0; i< list.size(); i++) {
            System.out.println(list.get(i));
        }
    }
    public static void writeListStringToLatex(List<String> listString, String fileName){
        Path filePath = Paths.get(fileName);
            if (!Files.exists(filePath)) {
                try {
                    Files.createFile(filePath);
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
            String add = "";
            for(int i = 0; i < listString.size();i++){
                add = add + listString.get(i)+"\n";
            }
            try {
                Files.write(filePath, add.getBytes(), StandardOpenOption.APPEND);
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
    }
    public static void generateEPLatex(List<String> fileContain, String homeFile, int[] index){
        String fileName = ReadFile.readhome("HOME_jMetal") + "/EP.tex";
        writeListStringToLatex(fileContain, fileName);
    }
    public static void generateIGDLatex(List<String> fileContain, String homeFile, int[] index){
        String fileName = ReadFile.readhome("HOME_jMetal") + "/IGD.tex";
        writeListStringToLatex(fileContain, fileName);
    }
    public static void generateIGDPlusLatex(List<String> fileContain, String homeFile, int[] index){
        String fileName = ReadFile.readhome("HOME_jMetal") + "/IGDPlus.tex";
        writeListStringToLatex(fileContain, fileName);
    }
    public static void generateComputeLatex(List<String> fileContain, String homeFile, int[] index){
        String fileName = ReadFile.readhome("HOME_jMetal") + "/ComputeTime.tex";
        writeListStringToLatex(fileContain, fileName);
    }
    public static List<String> findSecond(List<String> fileContain){
        List<String> tempContain = fileContain;
        for(int i = 8; i < tempContain.size()-3; i++) {
            String temp = tempContain.get(i).substring(12).replace(" ", "").replace("\\\\", "");
            String[] numbers = temp.split("&");
            double min = Double.POSITIVE_INFINITY;
            int indexMin = 0;
            for (int j = 0; j < numbers.length; j++) {
                if (numbers[j].contains("\\cellcolor{gray95}")) {
                    numbers[j] = numbers[j].replace("\\cellcolor{gray95}", "");
                    min = Double.parseDouble(numbers[j]);
                    indexMin = j;
                }
                //System.out.println(numbers[j]);
            }
            double secondMin = Double.POSITIVE_INFINITY;
            int secondeIndex = 0;
            for (int j = 0; j < numbers.length; j++){
                if (j != indexMin) {
                    secondMin = Math.min(secondMin, Double.parseDouble(numbers[j]));
                    //System.out.println(secondMin);
                    if (secondMin == Double.parseDouble(numbers[j])) {
                        secondeIndex = j;
                    }
                }
            }
            //numbers[secondeIndex] = "\\cellcolor{gray25}" + numbers[secondeIndex];
            tempContain.set(i,tempContain.get(i).replace(numbers[secondeIndex],"\\cellcolor{gray25}" + numbers[secondeIndex]));
        }
        return tempContain;
    }
    public static void main(String[] args) throws IOException {
        String[] listTexFile = listTexFile(listFolder(homeFile));
        String[] listComputeFile = listComputeFile(listFolder(homeFile));
        int[] index = {0,2};
        int headerComputeIndex = 0;
        int headerEPIndex = 11;
        int headerIGDIndex = 43;
        int headerIGDPlusIndex = 75;
        List<String> fileContainEP = fileTex(listTexFile, index, headerEPIndex);
        List<String> fileContainIGD = fileTex(listTexFile, index, headerIGDIndex);
        List<String> fileContainIGDPlus = fileTex(listTexFile, index, headerIGDPlusIndex);
        List<String> fileContainCompute = fileTex(listComputeFile, index, headerComputeIndex);
        generateEPLatex(fileContainEP,homeFile,index);
        generateIGDLatex(fileContainIGD,homeFile,index);
        generateIGDPlusLatex(fileContainIGDPlus,homeFile,index);
        generateComputeLatex(findSecond(fileContainCompute),homeFile,index);
    }
    public static void generateLatex(String store) throws IOException {
        homeFile = ReadFile.readhome("HOME_jMetal")+"/"+store;
        String[] listTexFile = listTexFile(listFolder(homeFile));
        String[] listComputeFile = listComputeFile(listFolder(homeFile));
        int[] index = {0,2};
        int headerComputeIndex = 0;
        int headerEPIndex = 11;
        int headerIGDIndex = 43;
        int headerIGDPlusIndex = 75;
        List<String> fileContainEP = fileTex(listTexFile, index, headerEPIndex);
        List<String> fileContainIGD = fileTex(listTexFile, index, headerIGDIndex);
        List<String> fileContainIGDPlus = fileTex(listTexFile, index, headerIGDPlusIndex);
        List<String> fileContainCompute = fileTex(listComputeFile, index, headerComputeIndex);
        generateEPLatex(fileContainEP,homeFile,index);
        generateIGDLatex(fileContainIGD,homeFile,index);
        generateIGDPlusLatex(fileContainIGDPlus,homeFile,index);
        generateComputeLatex(findSecond(fileContainCompute),homeFile,index);
    }
}
