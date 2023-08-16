package com.zx.test;

import java.io.IOException;

import DE.BSJ_yes;

public class test {
public static void main(String[] args) throws IOException {
	
	int minMapqUni= 10, maxCircle = 200000, minCircle = 140, linear_range_size_min = 50000, strigency = 2 ,relExp = 0;
	boolean intronLable = false;
	boolean mLable = false;
	String mitochondrion = "chrM";	
	String inputFile = "D:\\Ciri\\test\\BAM\\2\\sample_infor.tsv";
	String outputFile = "D:\\Ciri\\test\\BAM\\2\\sample_infor_mut_re.txt";
	String annotationFile = "F";
	String faFile = "D:\\Ciri\\test\\chr1.fa";
	//String UserGivecircRNA = "D:\\Ciri\\test\\tem.txt";
	String UserGivecircRNA = "";
	//SingleTest st = new SingleTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//st.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	//MutTest mut = new MutTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//mut.CIRI3(inputFile, outputFile, annotationFile, faFile, 2,"");
	//FileTest ft = new FileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//ft.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	//MutFileTest mft = new MutFileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//mft.CIRI3(inputFile, outputFile, annotationFile, faFile, 2, UserGivecircRNA);
	//FileTsvTest ft = new FileTsvTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//ft.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	//MutTsvFileTest mft = new MutTsvFileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable);
	//mft.CIRI3(inputFile, outputFile, annotationFile, faFile, 2, UserGivecircRNA);
	
	String inforpath = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/infor.tsv";
	String genepath = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/Gene_Expression.txt";
	String output = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/result.txt";
	BSJ_yes by = new BSJ_yes();		
	by.DEBSJNoFile(inforpath, genepath, output);
	
}
}
