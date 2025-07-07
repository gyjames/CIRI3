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

import com.zx.findcircrna.BamFindCircRNASTARScan1;
import com.zx.findcircrna.BamFindCircRNASTARScan2;
import com.zx.findcircrna.BamFindCircRNAScan2;
import com.zx.findcircrna.BamUserFindCircRNASTARScan2;
import com.zx.findcircrna.BamUserFindCircRNAScan2;
import com.zx.findcircrna.FindCircRNASTARScan1;
import com.zx.findcircrna.FindCircRNASTARScan2;
import com.zx.findcircrna.FindCircRNAScan2;
import com.zx.findcircrna.GetAnnotationInformation;
import com.zx.findcircrna.GetChimericOut;
import com.zx.findcircrna.GetUserCircRNA;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;
import com.zx.findcircrna.Summary;
import com.zx.findcircrna.UserFindCircRNASTARScan2;
import com.zx.findcircrna.UserFindCircRNAScan2;
import com.zx.hg38.Annotation;
import com.zx.hg38.AnnotationIntron;

public class FileSTARTest {
	private int minMapqUni,maxCircle,minCircle,linear_range_size_min,strigency,relExp;
	private boolean intronLable,mlable,spLable;
	private String mitochondrion;
	public FileSTARTest(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,int strigency,int relExp,String mitochondrion
			,boolean mlable,boolean spLable) {
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
		this.spLable = spLable;
	}

    public boolean CIRI3(String samFolder,String outPutFile,String annotationFile,String faFile,String UserGivecircRNA) throws IOException {
    	String outPutBSJCountFile = outPutFile+".BSJ_Matrix";
        String outPutFSJCountFile = outPutFile+".FSJ_Matrix";
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
 		//存储文件名-circRNA信息
 		//HashMap<String, ArrayList<String>> fileCircRNAMap = new HashMap<String, ArrayList<String>>();
		//保存文件名list
		HashMap<String, HashMap<String, String>> scan1IdTotalMap = new HashMap<String, HashMap<String, String>>();
		String fileClass = "sam";
 		ArrayList<String> filePathList = new ArrayList<String>();
 		ArrayList<String> fileNameList = new ArrayList<String>();
 		HashMap<String, Integer> fileSplitNumMap = new HashMap<String, Integer>();
 		BufferedReader samFileRead = new BufferedReader(new FileReader(new File(samFolder)));			
		String fileLine = samFileRead.readLine();
		while (fileLine != null) {
			if (fileLine.startsWith("#") || fileLine.equals("")) {
				fileLine = samFileRead.readLine();
				continue;
			}
			fileClass = fileLine.substring(fileLine.length()-3,fileLine.length());
			//String[] arrTem = fileLine.split(",");
			fileNameList.add(fileLine);			
			filePathList.add(fileLine);
			scan1IdTotalMap.put(fileLine, new HashMap<String, String>());
			fileSplitNumMap.put(fileLine, 1);
			fileLine = samFileRead.readLine();		
		}
		samFileRead.close();
		
		if (!fileClass.equals("sam") && !fileClass.equals("bam")) {
			System.out.println("Please enter the file that ends with sam or bam");
			return false;
		}
		//第一遍扫描开始
		int seqLen = 0;
		HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSiteSet = new HashSet<String>();
		HashMap<String, HashMap<String, Integer>> fileCircBSJMap = new HashMap<String, HashMap<String, Integer>>();
		HashMap<String, Integer> circBSJMap = new HashMap<String, Integer>();
		if (UserGivecircRNA.equals("")) {
			GetChimericOut getChiCirc = new GetChimericOut(maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
					chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable,spLable);
	 		for (int j = 0; j < filePathList.size(); j++) {
				String samFile = filePathList.get(j);
				System.out.println("Running:"+samFile);
				//输入文件
		        String[] samFileArr = samFile.split(",");
				String chimericPath = samFileArr[0];
				//String starSamPath = samFileArr[1];
				String unmappedSamPath = samFileArr[2];
				String outBSJPath = unmappedSamPath+"BSJ1";
				
				getChiCirc.getBSJ(chimericPath);
				HashMap<String, String> idCircMap = getChiCirc.getIdCircMap();
				BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(outBSJPath)));	
				for (String idKey : idCircMap.keySet()) {
					BSJOut.write(idCircMap.get(idKey)+"\n");
				}
				BSJOut.close();
				
				FindCircRNASTARScan1 scan1;
				if (unmappedSamPath.substring(unmappedSamPath.length()-3,unmappedSamPath.length()).equals("sam")) {
					scan1 = new FindCircRNASTARScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
							chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable,spLable);
					scan1.findCircRNAScan1(unmappedSamPath,idCircMap);
				}else if (unmappedSamPath.substring(unmappedSamPath.length()-3,unmappedSamPath.length()).equals("bam")){
					scan1 = new BamFindCircRNASTARScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
							chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable,spLable);
					scan1.findCircRNAScan1(unmappedSamPath,idCircMap);
				}else {
					System.out.println("Please enter the file that ends with sam or bam");
					return false;
				}
				//获取最长read长度
				if (seqLen < scan1.getReadLen()) {
					seqLen = scan1.getReadLen();
				}
				scan1 = null;
				idCircMap = null;
			}
	 		getChiCirc = null;
	 		
	 		System.out.println(df.format(System.currentTimeMillis())+" "+":First scan completed");  
	 		fileLog.write(df.format(System.currentTimeMillis())+" "+":First scan completed"+"\n");
	 		//导入文件		
			//存放每个sam文件第一遍扫描circRNA的ID Num	
			HashMap<String, String> scan1IdMap;
			for (int i = 0; i < filePathList.size(); i++) {
				scan1IdMap = new HashMap<String, String>();
				String samFile = filePathList.get(i);
				//输入文件
		        String[] samFileArr = samFile.split(",");
				String unmappedSamPath = samFileArr[2];
				BufferedReader BSJbr = new BufferedReader(new FileReader(new File(unmappedSamPath+"BSJ1")));
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
		}else {
			
			GetUserCircRNA guc = new GetUserCircRNA();
			chrCircSiteMap = guc.summaryUserCircRNA(UserGivecircRNA, chrTCGAMap);
			
			GetChimericOut getChiCirc = new GetChimericOut(maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
					chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable,spLable);
	 		for (int j = 0; j < filePathList.size(); j++) {
				String samFile = filePathList.get(j);
				//System.out.println("Running:"+samFile);
				//输入文件
		        String[] samFileArr = samFile.split(",");
				String chimericPath = samFileArr[0];
				//String starSamPath = samFileArr[1];
				String unmappedSamPath = samFileArr[2];
				String outBSJPath = unmappedSamPath+"BSJ1";
				
				getChiCirc.getBSJ(chimericPath);
				HashMap<String, String> idCircMap = getChiCirc.getIdCircMap();
				BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(outBSJPath)));	
				for (String idKey : idCircMap.keySet()) {
					BSJOut.write(idCircMap.get(idKey)+"\n");
				}
				BSJOut.close();
			}
	 		getChiCirc = null;
	 		//导入文件		
			//存放每个sam文件第一遍扫描circRNA的ID Num	
			HashMap<String, String> scan1IdMap;
			for (int i = 0; i < filePathList.size(); i++) {
				scan1IdMap = new HashMap<String, String>();
				circBSJMap = new HashMap<>();
				String samFile = filePathList.get(i);
				//输入文件
		        String[] samFileArr = samFile.split(",");
				String unmappedSamPath = samFileArr[2];
				BufferedReader BSJbr = new BufferedReader(new FileReader(new File(unmappedSamPath+"BSJ1")));
				String line = BSJbr.readLine();
				while (line != null) {
					String[] BSJArr = line.split("\t");
					scan1IdMap.put(BSJArr[0], "");
					String circKey = BSJArr[3]+"\t"+BSJArr[4]+"\t"+BSJArr[5];
					if(!circBSJMap.containsKey(circKey)) {
						circBSJMap.put(circKey, 1);
					}else {
						int temNum =  circBSJMap.get(circKey);
						circBSJMap.put(circKey, temNum+1);
					}
					line = BSJbr.readLine();
				}
				BSJbr.close();
				scan1IdTotalMap.put(samFile, scan1IdMap);
				fileCircBSJMap.put(samFile, circBSJMap);
			}
		}		
	    //第二遍扫描
	    //制作索引
		if (seqLen == 0) {
			seqLen = 500;
		}
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
		////////////
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
 		if (UserGivecircRNA.equals("")) {
 			for (int j = 0; j < filePathList.size(); j++) {
 				String samFile = filePathList.get(j);
 				System.out.println("Running:"+samFile);
 				
 				String[] samFileArr = samFile.split(",");
 				String starSamPath = samFileArr[1];
 				String unmappedSamPath = samFileArr[2];
 				String outBSJPath = unmappedSamPath+"BSJ1";
 				
 				FindCircRNAScan2 scan2 = null;
 		        FindCircRNASTARScan2 starScan2 = null;
 		        if (unmappedSamPath.substring(unmappedSamPath.length()-3,unmappedSamPath.length()).equals("sam")) {
 		        	scan2 = new FindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen,intronLable);  
 		        	scan2.setFSJScan2List();
 		    		scan2.findCircRNAScan2(unmappedSamPath,scan1IdTotalMap.get(samFile));
 				}else {
 					scan2 = new BamFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen,intronLable);  
 					scan2.setFSJScan2List();
 		    		scan2.findCircRNAScan2(unmappedSamPath,scan1IdTotalMap.get(samFile));
 				}	
 		        if (starSamPath.substring(starSamPath.length()-3,starSamPath.length()).equals("sam")) {
 		        	starScan2 = new FindCircRNASTARScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen,intronLable);  
 		        	starScan2.setFSJScan2List();
 		    		starScan2.findCircRNAScan2(starSamPath,scan1IdTotalMap.get(samFile),outBSJPath);
 				}else {
 					starScan2 = new BamFindCircRNASTARScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen,intronLable); 
 					starScan2.setFSJScan2List();
 					starScan2.findCircRNAScan2(starSamPath,scan1IdTotalMap.get(samFile),outBSJPath);
 				}
 		        //获取匹配的reads数目
				System.out.println(samFile+" Mapped_Reads"+" "+starScan2.getReadNum());  
				fileLog.write(samFile+" Mapped_Reads"+" "+starScan2.getReadNum()+"\n");
				starScan2.setReadNum();
				
 		        HashMap<String, Integer> circFSJNewTemMap = scan2.getCircFSJMap();
 		        HashMap<String, Integer> circFSJNewMap = starScan2.getCircFSJMap();
 		        for (String circKey : circFSJNewMap.keySet()) {
 					int num1 = circFSJNewMap.get(circKey);
 					int num2 = circFSJNewTemMap.get(circKey);
 					circFSJNewMap.put(circKey, num1+num2);
 				}
 				//FSJ
 				for (String circKey : circFSJNewMap.keySet()) {	
 					//System.out.println(circKey);
 					//合并信息
 					int numAll = circFSJMap.get(circKey);
 					int numTem = circFSJNewMap.get(circKey);
 					circFSJMap.put(circKey, numAll+numTem);
 					//输出FSJ矩阵
 					FSJmatrix[circRowMap.get(circKey)][j] = numTem;
 				}								
 				//重置FSJ
 				scan2 = null;
 				starScan2 = null;
 			}					
 			
 			chrSiteMap1 = null;
 			SiteMap1 = null;
 			siteList1 = null;
 			chrSiteMap2 = null;
 			SiteMap2 = null;
 			siteList2 = null;
 			siteArrayMap1 = null;
 			siteArrayMap2 = null;
 			scan1IdTotalMap = null;
 	        System.out.println(df.format(System.currentTimeMillis())+" "+":Second scan completed"); 
 	        fileLog.write(df.format(System.currentTimeMillis())+" "+":Second scan completed"+"\n");
 	        
 	        //整理BSJMatrix
 	        HashMap<String, Integer> circMap = new HashMap<String, Integer>();
 	        for (int i = 0; i < filePathList.size(); i++) {
 	        	circMap.clear();
 				String samFile = filePathList.get(i);
 				String[] samFileArr = samFile.split(",");
 				String unmappedSamPath = samFileArr[2];
 				BufferedReader BSJBr = new BufferedReader(new FileReader(new File(unmappedSamPath+"BSJ1")));	
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
 	        ArrayList<String> NewFilePathList = new ArrayList<String>();
 	       for (int i = 0; i < filePathList.size(); i++) {
				String samFile = filePathList.get(i);
				String[] samFileArr = samFile.split(",");
	 			String unmappedSamPath = samFileArr[2];
	 	        NewFilePathList.add(unmappedSamPath);
	 	        fileSplitNumMap.put(unmappedSamPath, 1);
			}  
 			Summary summary = new Summary(strigency,chrTCGAMap);
 			ArrayList<String> SummaryCircList = summary.summary(NewFilePathList,fileSplitNumMap,circFSJMap,UserGivecircRNA);
 			HashMap<String, String> circTrueIdMap = summary.getCircMap();
 			chrTCGAMap = null;
 			circFSJMap = null;
 			summary = null;
 			//输出每个样本CircRNA对应BSJnum
 	        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
 	        BSJCount.write("circRNA_ID"+"\t");
 	        for (String sampleName : fileNameList) {
 	        	BSJCount.write(sampleName+"\t");
 			}
 	        BSJCount.write("\n");
 	        //输出每个样本CircRNA对应FSJnum
 	        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
 	        FSJCount.write("circRNA_ID"+"\t");
 	        for (String sampleName : fileNameList) {
 	        	FSJCount.write(sampleName+"\t");
 			}
 	        FSJCount.write("\n");
 	        for (String circKey : circRowMap.keySet()) {
 				String circId = circKey.replaceFirst("\t", ":").replace("\t", "|");
 				if (circTrueIdMap.containsKey(circId)) {
 					BSJCount.write(circId);
 					FSJCount.write(circId);
 					for (int j = 0; j < fileNameList.size(); j++) {
 						BSJCount.write("\t"+BSJmatrix[circRowMap.get(circKey)][j]);
 						FSJCount.write("\t"+FSJmatrix[circRowMap.get(circKey)][j]);
 					}
 					BSJCount.write("\n");
 					FSJCount.write("\n");
 				}
 			}
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
		}else {
			UserFindCircRNAScan2 scan2 = null;
			UserFindCircRNASTARScan2 starScan2 = null;
 			for (int j = 0; j < filePathList.size(); j++) {
 				
 				String samFile = filePathList.get(j);
 				System.out.println("Running:"+samFile);
 				
 				String[] samFileArr = samFile.split(",");
 				String starSamPath = samFileArr[1];
 				String unmappedSamPath = samFileArr[2];

 		        if (unmappedSamPath.substring(unmappedSamPath.length()-3,unmappedSamPath.length()).equals("sam")) {
 		        	scan2 = new UserFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen);  
 		        	scan2.setFSJScan2List();
 		    		scan2.findCircRNAScan2(unmappedSamPath,scan1IdTotalMap.get(samFile));
 				}else {
 					scan2 = new BamUserFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen);  
 					scan2.setFSJScan2List();
 		    		scan2.findCircRNAScan2(unmappedSamPath,scan1IdTotalMap.get(samFile));
 				}	
 		        if (starSamPath.substring(starSamPath.length()-3,starSamPath.length()).equals("sam")) {
 		        	starScan2 = new UserFindCircRNASTARScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen);  
 		        	starScan2.setFSJScan2List();
 		    		starScan2.findCircRNAScan2(starSamPath,scan1IdTotalMap.get(samFile));
 				}else {
 					starScan2 = new BamUserFindCircRNASTARScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
 		    				chrTCGAMap,seqLen); 
 					starScan2.setFSJScan2List();
 					starScan2.findCircRNAScan2(starSamPath,scan1IdTotalMap.get(samFile));
 				}
 		        
 		        HashMap<String, Integer> circFSJNewTemMap = scan2.getCircFSJMap();
 	            HashMap<String, Integer> circFSJNewMap = starScan2.getCircFSJMap();
 	            
 	            HashMap<String, Integer> circBSJNewTemMap = scan2.getCircBSJMap();
 	            HashMap<String, Integer> circBSJNewMap = starScan2.getCircBSJMap();
 	            
 	            for (String circKey : circFSJNewMap.keySet()) {
 	    			int num1 = circFSJNewMap.get(circKey);
 	    			int num2 = circFSJNewTemMap.get(circKey);
 	    			circFSJNewMap.put(circKey, num1+num2);
 	    		}
 	            
 	            circBSJMap = fileCircBSJMap.get(samFile);
 	            for (String circKey : circBSJNewMap.keySet()) {
 	    			int num1 = circBSJNewMap.get(circKey);
 	    			int num2 = circBSJNewTemMap.get(circKey);
 	    			if(circBSJMap.containsKey(circKey)) {
 	    				circBSJNewMap.put(circKey, num1+num2+circBSJMap.get(circKey));
 	    			}else {
 	    				circBSJNewMap.put(circKey, num1+num2);
 	    			}
 	    			
 	    		}
 				//输出比对上的read
 				System.out.println(samFile+" Mapped_Reads"+" "+starScan2.getReadNum());  
				fileLog.write(samFile+" Mapped_Reads"+" "+starScan2.getReadNum()+"\n");
 				//FSJ
 				for (String circKey : circFSJNewMap.keySet()) {
 					FSJmatrix[circRowMap.get(circKey)][j] = circFSJNewMap.get(circKey);
 				}	
 				//BSJ
 				for (String circKey : circBSJNewMap.keySet()) {
 					BSJmatrix[circRowMap.get(circKey)][j] = circBSJNewMap.get(circKey);
 				}
 				//重置FSJ
 				scan2.setFSJScan2List();
 				scan2.setBSJScan2List();
 				starScan2.setFSJScan2List();
 				starScan2.setBSJScan2List();
 			}					
 			scan2 = null;
 			starScan2 = null;
 			chrSiteMap1 = null;
 			SiteMap1 = null;
 			siteList1 = null;
 			chrSiteMap2 = null;
 			SiteMap2 = null;
 			siteList2 = null;
 			siteArrayMap1 = null;
 			siteArrayMap2 = null;
 			scan1IdTotalMap = null;
 	        System.out.println(df.format(System.currentTimeMillis())+" "+":Second scan completed"); 
 	        fileLog.write(df.format(System.currentTimeMillis())+" "+":Second scan completed"+"\n");
 	        //输出每个样本CircRNA对应BSJnum
 	        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
 	        BSJCount.write("circRNA_ID"+"\t");
 	        for (String sampleName : fileNameList) {
 	        	BSJCount.write(sampleName+"\t");
 			}
 	        BSJCount.write("\n");
 	        //输出每个样本CircRNA对应FSJnum
 	        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
 	        FSJCount.write("circRNA_ID"+"\t");
 	        for (String sampleName : fileNameList) {
 	        	FSJCount.write(sampleName+"\t");
 			}
 	        FSJCount.write("\n");
 	        for (String circKey : circRowMap.keySet()) {
 				String circId = circKey.replaceFirst("\t", ":").replace("\t", "|");				
				BSJCount.write(circId);
				FSJCount.write(circId);
				for (int j = 0; j < fileNameList.size(); j++) {
					BSJCount.write("\t"+BSJmatrix[circRowMap.get(circKey)][j]);
					FSJCount.write("\t"+FSJmatrix[circRowMap.get(circKey)][j]);
				}
				BSJCount.write("\n");
				FSJCount.write("\n");			
 			}
 			BSJCount.close();	
 			FSJCount.close();
		} 
 		//删除临时文件
 		for (int j = 0; j < filePathList.size(); j++) {
			String samFile = filePathList.get(j);
			String[] samFileArr = samFile.split(",");
			String unmappedSamPath = samFileArr[2];
			new File(unmappedSamPath+"BSJ1").delete();
 		}
		long endTime = System.currentTimeMillis(); 
		System.out.println("Program run time:" + (endTime - startTime) + "ms");
		fileLog.write("Program run time:" + (endTime - startTime) + "ms"+"\n");
		fileLog.close();
		return true;
	}
}
