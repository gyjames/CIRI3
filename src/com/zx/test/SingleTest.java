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
import com.zx.findcircrna.GetUserCircRNA;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;
import com.zx.findcircrna.Summary;
import com.zx.findcircrna.UserFindCircRNAScan2;
import com.zx.hg19.Annotation;
import com.zx.hg19.AnnotationIntron;

public class SingleTest {
	
	private int minMapqUni,maxCircle,minCircle,linear_range_size_min,strigency,relExp;
	private boolean intronLable,mlable;
	private String mitochondrion;
	public SingleTest(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,int strigency,int relExp,String mitochondrion,boolean mlable) {
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
	public boolean CIRI3(String samFile,String outPutFile,String annotationFile,String faFile,String UserGivecircRNA) throws IOException {
		long startTime = System.currentTimeMillis(); 
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");  
		String outputFileLog = outPutFile+".log";
		String outPutBSJCountFile = outPutFile+".BSJ_Matrix";
        String outPutFSJCountFile = outPutFile+".FSJ_Matrix";
		BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
		System.out.println(df.format(System.currentTimeMillis())+" "+":CIRI3 start"); 
		fileLog.write(df.format(System.currentTimeMillis())+" "+":CIRI3 start"+"\n");
		HashMap<String, String> chrExonStartMap = new HashMap<String, String>(),chrExonEndMap = new HashMap<String, String>();
		HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
		HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
		HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>(),
				chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();
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
		int seqLen = 0;
		HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSiteSet = new HashSet<String>();
		HashMap<String, String> scan1IdMap = new HashMap<String, String>();
		if (UserGivecircRNA.equals("")) {
			//第一遍扫描开始
			FindCircRNAScan1 scan1;
			if (samFile.substring(samFile.length()-3,samFile.length()).equals("sam")) {
				scan1 = new FindCircRNAScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
						chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);
				scan1.findCircRNAScan1(samFile);
			}else if (samFile.substring(samFile.length()-3,samFile.length()).equals("bam")){
				scan1 = new BamFindCircRNAScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
						chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);
				scan1.findCircRNAScan1(samFile);
			}else {
				System.out.println("Please enter the file that ends with sam or bam");
				return false;
			}
			//获取最长read长度
			seqLen = scan1.getReadLen();
			//获取匹配的reads数目
			fileLog.write("Mapped_Reads"+" "+scan1.getReadNum()+"\n");
			scan1 = null;
			System.out.println(df.format(System.currentTimeMillis())+" "+":First scan completed");
			fileLog.write(df.format(System.currentTimeMillis())+" "+":First scan completed"+"\n");
			//导入文件			
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
			
		}else {
			GetUserCircRNA guc = new GetUserCircRNA();
			chrCircSiteMap = guc.summaryUserCircRNA(UserGivecircRNA, chrTCGAMap);
			seqLen = 500;
		}
		
		//制作索引，根据chr分组构建
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
		//     
		////////////
       //第一遍比对结束，第二遍扫描开始
        if (UserGivecircRNA.equals("")) {
        FindCircRNAScan2 scan2 = null;
        if (samFile.substring(samFile.length()-3,samFile.length()).equals("sam")) {
        	scan2 = new FindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
    				chrTCGAMap,seqLen,intronLable);  
    		scan2.findCircRNAScan2(samFile,scan1IdMap);
		}else {
			scan2 = new BamFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
    				chrTCGAMap,seqLen,intronLable);  
    		scan2.findCircRNAScan2(samFile,scan1IdMap);
		}	
		HashMap<String, Integer> circFSJNewMap = scan2.getCircFSJMap();
		scan2 = null;
		chrSiteMap1 = null;
		SiteMap1 = null;
		siteList1 = null;
		chrSiteMap2 = null;
		SiteMap2 = null;
		siteList2 = null;
		siteArrayMap1 = null;
		siteArrayMap2 = null;
		circFSJMap = null;
		scan1IdMap = null;
		System.out.println(df.format(System.currentTimeMillis())+" "+":Second scan completed");  
		fileLog.write(df.format(System.currentTimeMillis())+" "+":Second scan completed"+"\n");
		ArrayList<String> filePathList = new ArrayList<String>();
		filePathList.add(samFile);
		HashMap<String, Integer> fileSplitNumMap = new HashMap<String, Integer>();
		fileSplitNumMap.put(samFile, 1);
		Summary summary = new Summary(strigency,chrTCGAMap);
		ArrayList<String> SummaryCircList = summary.summary(filePathList,fileSplitNumMap,circFSJNewMap,UserGivecircRNA);
		chrTCGAMap = null;
		circFSJNewMap = null;
		summary = null;
		//new File(samFile+"BSJ1").delete();
	
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
            if (samFile.substring(samFile.length()-3,samFile.length()).equals("sam")) {
            	scan2 = new UserFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
        				chrTCGAMap,seqLen);  
        		scan2.findCircRNAScan2(samFile,scan1IdMap);
    		}else {
    			scan2 = new UserFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
        				chrTCGAMap,seqLen);  
        		scan2.findCircRNAScan2(samFile,scan1IdMap);
    		}	
    		HashMap<String, Integer> circFSJNewMap = scan2.getCircFSJMap();
    		HashMap<String, Integer> circBSJNewMap = scan2.getCircFSJMap();
    		//获取匹配的reads数目
			fileLog.write("Mapped_Reads"+" "+scan2.getReadNum()+"\n");
			scan2 = null;
			System.out.println(df.format(System.currentTimeMillis())+" "+":Scan completed");
			fileLog.write(df.format(System.currentTimeMillis())+" "+":Scan completed"+"\n");  		
    		BufferedWriter BSJBw = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
    		BufferedWriter FSJBw = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
    		BSJBw.write("circRNA"+"\t"+"BSJ"+"\n");
    		FSJBw.write("circRNA"+"\t"+"FSJ"+"\n");
    		for (String circKey : circBSJNewMap.keySet()) {
    			BSJBw.write(circKey+"\t"+circBSJNewMap.get(circKey)+"\n");
    			FSJBw.write(circKey+"\t"+circFSJNewMap.get(circKey)+"\n");
			}
    		BSJBw.close();
    		FSJBw.close();
		}
        long endTime = System.currentTimeMillis(); 
		System.out.println("Program run time:" + (endTime - startTime) + "ms");
		fileLog.write("Program run time:" + (endTime - startTime) + "ms"+"\n");
		fileLog.close();
		return true;		
}
}
