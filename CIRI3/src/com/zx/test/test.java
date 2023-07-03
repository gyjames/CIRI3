package com.zx.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class test {
public static void main(String[] args) throws IOException {
	
	/*int minMapqUni= 10, maxCircle = 200000, minCircle = 140, linear_range_size_min = 50000, strigency = 2 ,relExp = 0;
	boolean intronLable = false;
	boolean mLable = false;
	String mitochondrion = "chrM";
	SingleTest st = new SingleTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	String inputFile = "D:/Ciri/valid/aa.sam";
	String outputFile = "D:/Ciri/valid/CIRI3/aa_re.txt";
	//String inputFile = "D:\\zzh\\SRR3936695.allout.sam";
    //String outputFile = "D:\\zzh\\SRR3936695.allout_re.txt";
	String annotationFile = "D:\\BioInformation\\hg38\\gencode.v40.gtf";
	String faFile = "D:\\BioInformation\\hg38\\hg38.fa";
	//String UserGivecircRNA = "D:\\Ciri\\test\\tem.txt";
	String UserGivecircRNA = "";
	st.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	//MutTest mut = new MutTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//mut.CIRI3(inputFile, outputFile, annotationFile, faFile, 2, outputFileLog, UserGivecircRNA);*/
	BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File("D:/Ciri/valid/bb.txt"),true));
	BSJOut.write("aa");
	BSJOut.close();
}
}
