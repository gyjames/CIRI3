package com.zx.test;

import java.io.IOException;

public class test {
public static void main(String[] args) throws IOException {
	
	int minMapqUni= 10, maxCircle = 200000, minCircle = 140, linear_range_size_min = 50000, strigency = 2 ,relExp = 0;
	boolean intronLable = false;
	boolean mLable = false, spLable = false;
	String mitochondrion = "chrM";	
	String inputFile = "E:\\Ciri\\SATR\\STAR\\mut\\infor.tsv";
	//String inputFile = "E:/Ciri/SATR/STAR/Chimeric.out.junction,E:/Ciri/SATR/STAR/Aligned.out.sam,E:/Ciri/SATR/STAR/Unmapped.sam";
	String outputFile = "E:\\Ciri\\SATR\\STAR\\mut\\CIRI3_mut_user.txt";
	String annotationFile = "E:\\BioInformation\\hg38\\gencode_new.v40.gtf";
	String faFile = "E:\\BioInformation\\hg38\\hg38.fa";
	String UserGivecircRNA = "E:\\Ciri\\SATR\\STAR\\UserGivecircRNA.txt";
	//String UserGivecircRNA = "";
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
	
	/*String inforpath = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/infor.tsv";
	String genepath = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/Gene_Expression.txt";
	String output = "C:/Users/Administrator/OneDrive/桌面/data/DE/BSJ/With_RE/result.txt";
	BSJ_yes by = new BSJ_yes();		
	by.DEBSJNoFile(inforpath, genepath, output);*/
	//SingleSTARTest st = new SingleSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
	//st.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	//MutSTARTest mut = new MutSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
	//mut.CIRI3(inputFile, outputFile, annotationFile, faFile, 2,UserGivecircRNA);
	
	//FileSTARTest ft = new FileSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
	//ft.CIRI3(inputFile, outputFile, annotationFile, faFile, UserGivecircRNA);
	
	MutFileSTARTest mft = new MutFileSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
	mft.CIRI3(inputFile, outputFile, annotationFile, faFile, 2, UserGivecircRNA);
}
}
