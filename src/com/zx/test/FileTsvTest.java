package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.findcircrna.BamFindCircRNAScan1;
import com.zx.findcircrna.BamFindCircRNAScan2;
import com.zx.findcircrna.FindCircRNAScan1;
import com.zx.findcircrna.FindCircRNAScan2;
import com.zx.findcircrna.GetAnnotationInformation;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;
import com.zx.findcircrna.Summary;
import com.zx.hg38.Annotation;
import com.zx.hg38.AnnotationIntron;



public class FileTsvTest {
	private int minMapqUni,maxCircle,minCircle,linear_range_size_min,strigency,relExp;
	private boolean intronLable,mlable;
	private String mitochondrion;
	public FileTsvTest(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,int strigency,int relExp,String mitochondrion,boolean mlable) {
		super();
		this.minMapqUni = minMapqUni;
		this.maxCircle = maxCircle;
		this.minCircle = minCircle;
		this.linear_range_size_min = linear_range_size_min;
		this.intronLable = intronLable;
		this.strigency = strigency;
		this.relExp = relExp;
		this.mitochondrion = mitochondrion;
		this.mlable = mlable;

	}
    public boolean CIRI3(String samFileTsv,String outPutFile,String annotationFile,String faFile,String UserGivecircRNA) throws IOException {	    
    	String outPutBSJCountFile = outPutFile+".BSJ_Matrix";
        String outPutFSJCountFile = outPutFile+".FSJ_Matrix";
        String outPutcircRNAClass = outPutFile+".Score";
        String outputFileLog = outPutFile+".log";
    	long startTime = System.currentTimeMillis(); 
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");  
		BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
		System.out.println(df.format(System.currentTimeMillis())+" "+":CIRI3 start"); 
		fileLog.write(df.format(System.currentTimeMillis())+" "+":CIRI3 start"+"\n");
		HashMap<String, String> chrExonStartMap = new HashMap<String, String>(),chrExonEndMap = new HashMap<String, String>();
		HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
		HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
		HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>(),
				chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();;
		if (!annotationFile.equals("F")) {			
			GetAnnotationInformation GAI = new GetAnnotationInformation();	
			GAI.hand(annotationFile, intronLable);
			if (intronLable) {
				chrExonStartTranscriptMap = GAI.getChrExonStartTranscriptMap();
				chrExonEndTranscriptMap = GAI.getChrExonEndTranscriptMap();
				geneExonMap = GAI.getGeneExonMap();
				exonListMap = GAI.getExonListMap();
			}else {
				chrExonStartMap = GAI.getChrExonStartMap();
				chrExonEndMap = GAI.getChrExonEndMap();
				geneExonMap = GAI.getGeneExonMap();
				exonListMap = GAI.getExonListMap();
			}
			GAI=null;
			if(geneExonMap.size()==0) {
				System.out.println("please input formatted annotation file");
				fileLog.write(df.format(System.currentTimeMillis())+" "+":please input formatted annotation file"+"\n");
				return false;
			}
			System.out.println(df.format(System.currentTimeMillis())+" "+":Successfully imported comment files");	
			fileLog.write(df.format(System.currentTimeMillis())+" "+":Successfully imported comment files"+"\n");		
		}		
		//导入参考基因组文件	
		ReadFaFile RF = new ReadFaFile();
		RF.readFa(faFile);
		//获取最长染色体长度
		HashMap<String, Integer> chrLenMap = RF.getChrLenMap();	
		//获取Chr-Seq
		HashMap<String, String> chrTCGAMap = RF.getChrTCGAMap();
		RF = null;
		System.out.println(df.format(System.currentTimeMillis())+" "+":Successful import of reference genome files");  
		fileLog.write(df.format(System.currentTimeMillis())+" "+":Successful import of reference genome files"+"\n");
		
 		//第一遍扫描开始		
 		//记录哪些文件是Rnase
 		HashMap<String, ArrayList<String>> isRnaseFileMap = new HashMap<String, ArrayList<String>>();
 		//RNase：0：未进行RNase处理  1：进行RNase处理	2：是否进行RNase处理未知
 		isRnaseFileMap.put("0", new ArrayList<String>());
 		isRnaseFileMap.put("1", new ArrayList<String>());
 		isRnaseFileMap.put("2", new ArrayList<String>());
 		//存储文件名-circRNA信息
		//保存文件名list
 		String fileClass = "sam";
 		ArrayList<String> filePathListTem = new ArrayList<String>();
 		HashMap<String, Integer> fileSplitNumMap = new HashMap<String, Integer>();
 		BufferedReader samFileRead = new BufferedReader(new FileReader(new File(samFileTsv)));			
		String fileLine = samFileRead.readLine();
		while (fileLine != null) {
			if (fileLine.startsWith("#") || fileLine.equals("")) {
				fileLine = samFileRead.readLine();
				continue;
			}
			String[] lineArr = fileLine.split("\t");
			fileClass = lineArr[0].substring(lineArr[0].length()-3,lineArr[0].length());
			filePathListTem = isRnaseFileMap.get(lineArr[1]);	
			filePathListTem.add(lineArr[0]);
			isRnaseFileMap.put(lineArr[1], filePathListTem);
			fileSplitNumMap.put(lineArr[0], 1);
			fileLine = samFileRead.readLine();		
		}
		samFileRead.close();
		if (!fileClass.equals("sam") && !fileClass.equals("bam")) {
			System.out.println("Please enter the file that ends with sam or bam");
			return false;
		}
		int unRnaseFileNum = isRnaseFileMap.get("0").size();
	    int RnaseFileNum = isRnaseFileMap.get("1").size();   
	    ArrayList<String> filePathList = new ArrayList<String>();
	    ArrayList<String> fileNameList = new ArrayList<String>();
		for (int i = 0; i <= 2; i++) {
			filePathListTem = isRnaseFileMap.get(i+"");
			for (String samfile : filePathListTem) {
				filePathList.add(samfile);
				String[] arrTem = samfile.split("/");
				fileNameList.add(arrTem[arrTem.length-1]);
			}
			
		}
		filePathListTem = null;
				 		
		//第一遍扫描开始
		int seqLen = 0;
		long RNaseReadNum = 0;
		long unRNaseReadNum = 0;
		FindCircRNAScan1 scan1;
		if (fileClass.equals("sam")) {
			scan1 = new FindCircRNAScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
					chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);	
		}else {
			scan1 = new BamFindCircRNAScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
					chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);	
		}			
		for (int j = 0; j < filePathList.size(); j++) {
			String samFile = filePathList.get(j);
			System.out.println("Running:"+samFile);  
			scan1.findCircRNAScan1(samFile);
			//获取匹配的reads数目
			long matchNumTem = scan1.getReadNum();
			System.out.println(samFile+" Mapped_Reads"+" "+matchNumTem);  
			fileLog.write(samFile+" Mapped_Reads"+" "+matchNumTem+"\n");
			if (j < unRnaseFileNum) {
				unRNaseReadNum += matchNumTem;
			}else if (j < unRnaseFileNum + RnaseFileNum) {
				RNaseReadNum += matchNumTem;
			}
			scan1.setReadNum();
		}		
 		//获取最长read长度
		if (seqLen < scan1.getReadLen()) {
			seqLen = scan1.getReadLen();
		}
		scan1 = null;
 		System.out.println(df.format(System.currentTimeMillis())+" "+":First scan completed");  
 		fileLog.write(df.format(System.currentTimeMillis())+" "+":First scan completed"+"\n");
 		//导入文件
		HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSiteSet = new HashSet<String>();
		//存放每个sam文件第一遍扫描circRNA的ID Num
		HashMap<String, HashMap<String, String>> scan1IdTotalMap = new HashMap<String, HashMap<String, String>>();	
		HashMap<String, String> scan1IdMap;
		for (int i = 0; i < filePathList.size(); i++) {
			scan1IdMap = new HashMap<String, String>();
			String samFile = filePathList.get(i);
			BufferedReader BSJbr = new BufferedReader(new FileReader(new File(samFile+"BSJ1")));
			String line = BSJbr.readLine();
			while (line != null) {
				String[] BSJArr = line.split("\t",5);
				scan1IdMap.put(BSJArr[0], "");
				if(!chrCircSiteMap.containsKey(BSJArr[3])) {
					circSiteSet = new HashSet<String>();
					circSiteSet.add(BSJArr[4]);
					chrCircSiteMap.put(BSJArr[3], circSiteSet);
				}else {
					circSiteSet=chrCircSiteMap.get(BSJArr[3]);
					circSiteSet.add(BSJArr[4]);
					chrCircSiteMap.put(BSJArr[3], circSiteSet);
				}
				line = BSJbr.readLine();
			}
			BSJbr.close();
			scan1IdTotalMap.put(samFile, scan1IdMap);
		}
		
 		//第二遍扫描
 		//制作索引
		seqLen = seqLen-12;		
		HashMap<String, Integer> circFSJMap = new HashMap<String, Integer>();
		HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
		HashMap<Integer, ArrayList<SiteSort>> SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
		ArrayList<SiteSort> siteList1 = new ArrayList<SiteSort>();
		HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
		HashMap<Integer, ArrayList<SiteSort>> SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
		ArrayList<SiteSort> siteList2 = new ArrayList<SiteSort>();
		
		HashMap<String, byte[]> siteArrayMap1 = new HashMap<String, byte[]>();
		HashMap<String, byte[]> siteArrayMap2 = new HashMap<String, byte[]>();
		for (String chrKey : chrCircSiteMap.keySet()) {
			circSiteSet = chrCircSiteMap.get(chrKey);
			byte[] siteArray1 = new byte[(chrLenMap.get(chrKey)/seqLen)+1];
			byte[] siteArray2 = new byte[(chrLenMap.get(chrKey)/seqLen)+1];
			SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
			SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
			for (String siteInfor : circSiteSet) {
				String[] arr = siteInfor.split("\t");
				String temCircRNA = chrKey+"\t"+arr[0]+"\t"+arr[1];
				circFSJMap.put(temCircRNA, 0);
				int site1 = Integer.parseInt(arr[0])/seqLen;
				int site2 = Integer.parseInt(arr[1])/seqLen;
				siteArray1[site1] = 1;
				siteArray2[site2] = 1;
				////////////////////////////
				if (!SiteMap1.containsKey(site1)) {
					////////// map1
					siteList1 = new ArrayList<SiteSort>();
					siteList1.add(new SiteSort(Integer.parseInt(arr[0]),arr));
					SiteMap1.put(site1, siteList1);
				} else {
					////////// map1
					siteList1 = SiteMap1.get(site1);
					siteList1.add(new SiteSort(Integer.parseInt(arr[0]),arr));
					SiteMap1.put(site1, siteList1);
				}
				if (!SiteMap2.containsKey(site2)) {
					////////// map2
					siteList2 = new ArrayList<SiteSort>();
					siteList2.add(new SiteSort(Integer.parseInt(arr[1]),arr));
					SiteMap2.put(site2, siteList2);
				} else {
					////////// map2
					siteList2 = SiteMap2.get(site2);
					siteList2.add(new SiteSort(Integer.parseInt(arr[1]),arr));
					SiteMap2.put(site2, siteList2);
				}
			}
			chrSiteMap1.put(chrKey, SiteMap1);
			chrSiteMap2.put(chrKey, SiteMap2);
			siteArrayMap1.put(chrKey, siteArray1);
			siteArrayMap2.put(chrKey, siteArray2);
		}
		chrCircSiteMap = null;
		circSiteSet = null;			
		/////////////排序//////////////////////
		for (String chrKey : chrSiteMap1.keySet()) {
			SiteMap1 = chrSiteMap1.get(chrKey);
			for (Integer site : SiteMap1.keySet()) {
				siteList1 = SiteMap1.get(site);
				Collections.sort(siteList1);
				SiteMap1.put(site, siteList1);
			}
			chrSiteMap1.put(chrKey, SiteMap1);
		}
        for (String chrKey : chrSiteMap2.keySet()) {
        	SiteMap2 = chrSiteMap2.get(chrKey);
        	for (Integer site : SiteMap2.keySet()) {
            	siteList2 = SiteMap2.get(site);
    			Collections.sort(siteList2);
    			SiteMap2.put(site, siteList2);
    		} 
        	chrSiteMap2.put(chrKey, SiteMap2);
		}
        //整理矩阵输出
        //创建circRNA-文件矩阵
        int[][] BSJmatrix = new int[circFSJMap.keySet().size()][filePathList.size()];
		 //创建FSJ-文件矩阵
		int[][] FSJmatrix = new int[circFSJMap.keySet().size()][filePathList.size()];
 		//circRNA file_row
 		HashMap<String, Integer> circRowMap = new HashMap<String, Integer>();
 		int circNum = 0;
 		for (String circKey : circFSJMap.keySet()) {
 			circRowMap.put(circKey, circNum);
 			circNum++;
		}
		//第一遍比对结束，第二遍扫描开始
	    //存放第一遍扫描circRNA信息
 		FindCircRNAScan2 scan2;
 		if (fileClass.equals("sam")) {
 			scan2 = new FindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 					chrTCGAMap,seqLen,intronLable); 
		}else {
			scan2 = new BamFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
					chrTCGAMap,seqLen,intronLable); 
		} 		
 		for (int j = 0; j < filePathList.size(); j++) {
			String samFile = filePathList.get(j);
			System.out.println("Running:"+samFile);  
			scan2.findCircRNAScan2(samFile,scan1IdTotalMap.get(samFile));			
			//FSJ
			HashMap<String, Integer> circNewFSJMap = scan2.getCircFSJMap();
			for (String circKey : circNewFSJMap.keySet()) {				
				//合并信息
				int numAll = circFSJMap.get(circKey);
				int numTem = circNewFSJMap.get(circKey);
				circFSJMap.put(circKey, numAll+numTem);
				//输出FSJ矩阵
				FSJmatrix[circRowMap.get(circKey)][j] = numTem;
			}	
			//重置FSJ
			scan2.setFSJScan2List();
		}		
 		scan2 = null;
		chrSiteMap1 = null;
		SiteMap1 = null;
		siteList1 = null;
		chrSiteMap2 = null;
		SiteMap2 = null;
		siteList2 = null;
		siteArrayMap1 = null;
		siteArrayMap2 = null;
		scan1IdMap = null;
		scan1IdTotalMap = null;
        System.out.println(df.format(System.currentTimeMillis())+" "+":Second scan completed");  
        fileLog.write(df.format(System.currentTimeMillis())+" "+":Second scan completed"+"\n");
        //整理BSJMatrix
        HashMap<String, Integer> circMap = new HashMap<String, Integer>();
        for (int i = 0; i < filePathList.size(); i++) {
        	circMap.clear();
			String samFile = filePathList.get(i);
			BufferedReader BSJBr = new BufferedReader(new FileReader(new File(samFile+"BSJ1")));	
			String line = BSJBr.readLine();
			while (line != null ) {					
				String[] circLineArr = line.split("\t",7);				
				if (circLineArr[2].equals("1")) {
					String chrStartEnd = circLineArr[3] + "\t" + circLineArr[4] + "\t" + circLineArr[5];
					if (!circMap.containsKey(chrStartEnd)) {
						circMap.put(chrStartEnd, 1);
					} else {
						int temNum = circMap.get(chrStartEnd);
						circMap.put(chrStartEnd, temNum+1);
					}
				}
				line = BSJBr.readLine();					
			}
			BSJBr.close();
			for (String circKey : circMap.keySet()) {
				BSJmatrix[circRowMap.get(circKey)][i] = circMap.get(circKey);
			}
		}  
        
        //合并总结
  		Summary summary = new Summary(strigency,chrTCGAMap);
  		ArrayList<String> SummaryCircList = summary.summary(filePathList,fileSplitNumMap,circFSJMap,UserGivecircRNA);
  		HashMap<String, String> circTrueIdMap = summary.getCircMap();
  		chrTCGAMap = null;
  		circFSJMap = null;
  		summary = null;
  		fileSplitNumMap = null;
        //整理矩阵输出				
        //输出每个样本CircRNA对应BSJnum
        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
        BSJCount.write("circRNA_ID");
        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
        FSJCount.write("circRNA_ID");
        for (int i = 0; i < fileNameList.size(); i++) {
			BSJCount.write("\t"+fileNameList.get(i));
			FSJCount.write("\t"+fileNameList.get(i));
		}
        BSJCount.write("\n");
        FSJCount.write("\n");
		BufferedWriter circRNAClass = new BufferedWriter(new FileWriter(new File(outPutcircRNAClass)));
		circRNAClass.write("circRNA_ID	RnaseBSJ	UnRnaseBSJ	RnaseCPM	UnRnaseCPM	Score"+"\n");
		
		for (String circKey : circRowMap.keySet()) {
			String circId = circKey.replaceFirst("\t", ":").replace("\t", "|");
			if (circTrueIdMap.containsKey(circId)) {
				BSJCount.write(circId);
				FSJCount.write(circId);
				circRNAClass.write(circId);
				for (int j = 0; j < fileNameList.size(); j++) {
					BSJCount.write("\t"+BSJmatrix[circRowMap.get(circKey)][j]);
					FSJCount.write("\t"+FSJmatrix[circRowMap.get(circKey)][j]);
				}
				BSJCount.write("\n");
				FSJCount.write("\n");
				int rnaseBSJNum = 0;
				int unrnaseBSJNum = 0;				
				for (int j = 0; j < unRnaseFileNum; j++) {
					unrnaseBSJNum += BSJmatrix[circRowMap.get(circKey)][j];
				}
				for (int j = unRnaseFileNum; j < RnaseFileNum+unRnaseFileNum; j++) {
					rnaseBSJNum += BSJmatrix[circRowMap.get(circKey)][j];
				}
				if(unrnaseBSJNum == 0) {
					circRNAClass.write("\t"+rnaseBSJNum+"\t"+unrnaseBSJNum+"\t"+(rnaseBSJNum*2/((RNaseReadNum+0.0)/1000000))+"\t"+0+"\t"+"NA"+"\n");
				}else {
					double RnaseRatio = rnaseBSJNum*2/((RNaseReadNum+0.0)/1000000);
					double unRnaseRatio	= unrnaseBSJNum*2/((unRNaseReadNum+0.0)/1000000);
					double Ratio = RnaseRatio/unRnaseRatio;										
					circRNAClass.write("\t"+rnaseBSJNum+"\t"+unrnaseBSJNum+"\t"+String.format("%.2f", RnaseRatio )+"\t"+String.format("%.2f", unRnaseRatio )
					+"\t"+String.format("%.2f", Ratio )+"\n");						
				}			
			}
		}
		circRNAClass.close();
		BSJCount.close();	
		FSJCount.close();
		FSJmatrix = null;
		BSJmatrix = null;
		fileNameList = null;
		circTrueIdMap = null;
		circRowMap = null;
		//注释circRNA
		if (intronLable) {
			AnnotationIntron annotation = new AnnotationIntron();
			annotation.annotation(SummaryCircList, geneExonMap, chrExonStartTranscriptMap, chrExonEndTranscriptMap, exonListMap, outPutFile,relExp);
		}else {
			Annotation annotation = new Annotation();
			annotation.annotation(SummaryCircList, geneExonMap, chrExonStartMap, chrExonEndMap, exonListMap, outPutFile,relExp);
		}		
		System.out.println(df.format(System.currentTimeMillis())+" "+":Collation of circRNA completed"); 
		fileLog.write(df.format(System.currentTimeMillis())+" "+":Collation of circRNA completed"+"\n");
		//删除临时文件
 		for (int j = 0; j < filePathList.size(); j++) {
			String samFile = filePathList.get(j);
			new File(samFile+"BSJ1").delete();
 		}
		long endTime = System.currentTimeMillis(); 
		System.out.println("Program run time:" + (endTime - startTime) + "ms");
		fileLog.write("Program run time:" + (endTime - startTime) + "ms"+"\n");
		fileLog.close();
		return true;
	}
}
