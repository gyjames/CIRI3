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
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.zx.findcircrna.BamToSam;
import com.zx.findcircrna.GetAnnotationInformation;
import com.zx.findcircrna.GetUserCircRNA;
import com.zx.findcircrna.MutFindCircRNAScan1;
import com.zx.findcircrna.MutFindCircRNAScan2;
import com.zx.findcircrna.MutUserFindCircRNAScan2;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;
import com.zx.findcircrna.Summary;
import com.zx.hg38.Annotation;
import com.zx.hg38.AnnotationIntron;

public class MutFileTest {
	private int minMapqUni,maxCircle,minCircle,linear_range_size_min,strigency,relExp,seqLen = 0,AllFileSplitNum;
	private long matchNum = 0;
	private boolean intronLable,mlable,isSam;
	private String mitochondrion,UserGivecircRNAG;
	public MutFileTest(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,int strigency,int relExp,String mitochondrion,boolean mlable) {
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
	public static String samFile;	
	public static ArrayList<String> filePathList = new ArrayList<String>();
	private static Lock lock = new ReentrantLock();
	private static HashMap<String, Integer> circFSJMap,circBSJMap,fileSplitNumMap;
	private static HashMap<String, String> chrTCGAMap;
	HashMap<String, String> chrExonStartMap = new HashMap<String, String>(),chrExonEndMap = new HashMap<String, String>();
	HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>(),
			chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();
	private static HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,chrSiteMap2;
	public static HashMap<String, byte[]> siteArrayMap1,siteArrayMap2;
	public boolean CIRI3(String samFolder,String outPutFile,String annotationFile,String faFile,int threads,String UserGivecircRNA) throws IOException {
		long startTime = System.currentTimeMillis(); 
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");  
		String outPutBSJCountFile = outPutFile+".BSJ_Matrix";
        String outPutFSJCountFile = outPutFile+".FSJ_Matrix";
        String outputFileLog = outPutFile+".log";
		BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
		System.out.println(df.format(System.currentTimeMillis())+" "+":CIRI3 start"); 
		fileLog.write(df.format(System.currentTimeMillis())+" "+":CIRI3 start"+"\n");
		UserGivecircRNAG = UserGivecircRNA;
		//线程池，设置线程数量
		ExecutorService poolExe = Executors.newFixedThreadPool(threads);				
		// 设置主程序和子程序锁
		final CyclicBarrier threadSub = new CyclicBarrier(threads+1);
		final CyclicBarrier threadMain = new CyclicBarrier(threads+1);
		AtomicInteger incr = new AtomicInteger(1);
		for (int i = 0; i < threads; i++) {
			Runnable runnable = new Runnable() {
				public void run() {
					try {
						threadMain.await();
						if (UserGivecircRNAG.equals("")) {
							//第一遍扫描									
					 		//存放每个sam文件第一遍扫描circRNA的ID Num						
							MutFindCircRNAScan1 scan1 = new MutFindCircRNAScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable,chrExonStartMap,
									chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);						
					 		for (int j = 0; j < filePathList.size()-1; j++) {
					 			while(true) {
									int threadNum= incr.getAndIncrement();
									if (threadNum > AllFileSplitNum) {
										break;
									}else {
										scan1.findCircRNAScan1(samFile,AllFileSplitNum,threadNum);																
									}
								}
					 			//获取最长read长度
								long matchNumTem = scan1.getReadNum();
								lock.lock();
								matchNum = matchNum + matchNumTem;	
								lock.unlock();
								scan1.setReadNum();
								threadSub.await();
								threadMain.await();					
							}	
					 		while(true) {
								int threadNum= incr.getAndIncrement();
								if (threadNum > AllFileSplitNum) {
									break;
								}else {
									scan1.findCircRNAScan1(samFile,AllFileSplitNum,threadNum);																
								}
							}
					 		//获取最长read长度
							int seqLenTem = scan1.getReadLen();
							long matchNumTem = scan1.getReadNum();
							lock.lock();
							if (seqLen < seqLenTem) {
								seqLen = seqLenTem;
							}
							matchNum = matchNum + matchNumTem;	
							lock.unlock();
							//清除第一遍的缓存	
							scan1 = null;
							threadSub.await();
							threadMain.await();						
							//第二遍扫描
							MutFindCircRNAScan2 scan2 = new MutFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
				    				chrTCGAMap,seqLen,intronLable);  
							HashMap<String, String> scan1IdMap = new HashMap<String, String>();
							for (int j = 0; j < filePathList.size(); j++) {
								while(true) {
									int threadNum = incr.getAndIncrement();
									if (threadNum > AllFileSplitNum) {
										break;
									}else {
										scan1IdMap.clear();
										BufferedReader BSJbr = new BufferedReader(new FileReader(new File(samFile+"BSJ"+threadNum)));
										String line = BSJbr.readLine();
										while (line != null) {
											String[] BSJArr = line.split("\t",2);
											scan1IdMap.put(BSJArr[0], "");
											line = BSJbr.readLine();
										}
										BSJbr.close();
										scan2.findCircRNAScan2(samFile,scan1IdMap,AllFileSplitNum,threadNum);																
									}								
								}
								HashMap<String, Integer> circFSJMapTem = scan2.getCircFSJMap();
								lock.lock();
								for (String circKey : circFSJMap.keySet()) {
									int num = circFSJMap.get(circKey);						
									int numNew = circFSJMapTem.get(circKey);
									circFSJMap.put(circKey, num+numNew);
								}	
								lock.unlock();
								scan2.setFSJScan2List();
								threadSub.await();
								threadMain.await();								
							}
							scan2 = null;
							scan1IdMap = null; 
						}else {
							//第二遍扫描
							MutUserFindCircRNAScan2 scan2 = new MutUserFindCircRNAScan2(minMapqUni,circFSJMap,linear_range_size_min,siteArrayMap1,siteArrayMap2,chrSiteMap1,chrSiteMap2,
				    				chrTCGAMap,seqLen);  
							HashMap<String, String> scan1IdMap = new HashMap<String, String>();
							for (int j = 0; j < filePathList.size(); j++) {
								while(true) {
									int threadNum = incr.getAndIncrement();
									if (threadNum > AllFileSplitNum) {
										break;
									}else {
										scan2.findCircRNAScan2(samFile,scan1IdMap,AllFileSplitNum,threadNum);																
									}								
								}
								HashMap<String, Integer> circFSJMapTem = scan2.getCircFSJMap();
								HashMap<String, Integer> circBSJMapTem = scan2.getCircBSJMap();
								lock.lock();
								//获取最长read长度
								long matchNumTem = scan2.getReadNum();		
								matchNum = matchNum + matchNumTem;		
								
								for (String circKey : circFSJMap.keySet()) {
									int num = circFSJMap.get(circKey);						
									int numNew = circFSJMapTem.get(circKey);
									circFSJMap.put(circKey, num+numNew);
									num = circBSJMap.get(circKey);						
									numNew = circBSJMapTem.get(circKey);
									circBSJMap.put(circKey, num+numNew);
								}	
								lock.unlock();
								scan2.setFSJScan2List();
								scan2.setBSJScan2List();
								threadSub.await();
								threadMain.await();								
							}
							scan2 = null;
							scan1IdMap = null; 
						}						
                      } catch (Exception e) {
						e.printStackTrace();
					}
				}
			};
			poolExe.execute(runnable); // 运动员开始任务
		}
		try {
			// 整理注释文件	
			HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
			HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
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
			chrTCGAMap = RF.getChrTCGAMap();
			RF = null;
			System.out.println(df.format(System.currentTimeMillis())+" "+":Successful import of reference genome files");
			fileLog.write(df.format(System.currentTimeMillis())+" "+":Successful import of reference genome files"+"\n");
						
			//保存文件名list
			filePathList = new ArrayList<String>();
	 		ArrayList<String> fileNameList = new ArrayList<String>();
	 		BufferedReader samFileRead = new BufferedReader(new FileReader(new File(samFolder)));			
			String fileLine = samFileRead.readLine();
			while (fileLine != null) {
				if (fileLine.startsWith("#") || fileLine.equals("")) {
					fileLine = samFileRead.readLine();
					continue;
				}
				String[] arrTem = fileLine.split("/");
				fileNameList.add(arrTem[arrTem.length-1]);			
				if (fileLine.substring(fileLine.length()-3,fileLine.length()).equals("sam")) {
					filePathList.add(fileLine);
				}else if (fileLine.substring(fileLine.length()-3,fileLine.length()).equals("bam")){
					isSam = true;
					samFile = fileLine.substring(0,fileLine.length()-3)+"sam";
					BamToSam bts = new BamToSam();
					bts.bamToBam(fileLine, samFile);
					filePathList.add(samFile);
				}else {
					System.out.println("Please enter the file that ends with sam or bam");
					return false;
				}
				fileLine = samFileRead.readLine();		
			}
			samFileRead.close();
			//存储文件的划分
			fileSplitNumMap = new HashMap<String, Integer>();
			HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
			HashSet<String> circSiteSet = new HashSet<String>();
			if (UserGivecircRNAG.equals("")) {			
				//存储文件名-circRNA信息
			    for (int i = 0; i < filePathList.size(); i++) {
			    	samFile = filePathList.get(i);
			    	matchNum = 0;
			    	//更新原子
					incr.set(1);
			    	//根据文件大小拆分文件
					File file = new File(samFile);
			    	long fileSizeGB = file.length()/1024/1024/1024;
			    	if (fileSizeGB > 200 && fileSizeGB > threads*10) {
			    		AllFileSplitNum = (int) fileSizeGB/10;
					}else {
						AllFileSplitNum = threads;
					}	
			    	fileSplitNumMap.put(samFile, AllFileSplitNum);
					System.out.println(df.format(System.currentTimeMillis())+" "+"Running:"+samFile);  
					fileLog.write(df.format(System.currentTimeMillis())+" "+"Running:"+samFile+"\n");
					threadSub.reset();//更新子线程
				    threadMain.await();//所有线程激活2
				    threadMain.reset();//更新主线程
					threadSub.await();//主线程关闭，子线程激活2  
					System.out.println(samFile+" Mapped_Reads"+" "+matchNum);  
					fileLog.write(samFile+" Mapped_Reads"+" "+matchNum+"\n");
			    }	
			    System.out.println(df.format(System.currentTimeMillis())+" "+":First scan completed");  
			    fileLog.write(df.format(System.currentTimeMillis())+" "+":First scan completed"+"\n");			    
			    //导入文件			
				//存放每个sam文件第一遍扫描circRNA的ID Num
				for (int i = 0; i < filePathList.size(); i++) {
					samFile = filePathList.get(i);
					for (int j = 1; j <= fileSplitNumMap.get(samFile); j++) {					
						BufferedReader BSJbr = new BufferedReader(new FileReader(new File(samFile+"BSJ"+j)));
						String line = BSJbr.readLine();
						while (line != null) {
							String[] BSJArr = line.split("\t",5);
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
					}				
				}			    
			}else {
				for (int i = 0; i < filePathList.size(); i++) {
			    	samFile = filePathList.get(i);
			    	//根据文件大小拆分文件
					File file = new File(samFile);
			    	long fileSizeGB = file.length()/1024/1024/1024;
			    	if (fileSizeGB > 200 && fileSizeGB > threads*10) {
			    		AllFileSplitNum = (int) fileSizeGB/10;
					}else {
						AllFileSplitNum = threads;
					}	
			    	fileSplitNumMap.put(samFile, AllFileSplitNum);
			    }	
				GetUserCircRNA guc = new GetUserCircRNA();
				chrCircSiteMap = guc.summaryUserCircRNA(UserGivecircRNA, chrTCGAMap);
				seqLen = 500;
			}
			//第二遍扫描
		    //制作索引
			seqLen = seqLen-12;		
			circFSJMap = new HashMap<String, Integer>();
			circBSJMap = new HashMap<String, Integer>();
			HashMap<String, Integer> circSumFSJMap = new HashMap<String, Integer>();
			chrSiteMap1 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
			HashMap<Integer, ArrayList<SiteSort>> SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
			ArrayList<SiteSort> siteList1 = new ArrayList<SiteSort>();
			chrSiteMap2 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
			HashMap<Integer, ArrayList<SiteSort>> SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
			ArrayList<SiteSort> siteList2 = new ArrayList<SiteSort>();
			
			siteArrayMap1 = new HashMap<String, byte[]>();
			siteArrayMap2 = new HashMap<String, byte[]>();
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
					circSumFSJMap.put(temCircRNA, 0);
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
			circBSJMap.putAll(circFSJMap);
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
	 		if (UserGivecircRNAG.equals("")) {
	 			//进行第二遍扫描   	  
		        for (int j = 0; j < filePathList.size(); j++) {
					 samFile = filePathList.get(j);
					 incr.set(1);
					 AllFileSplitNum = fileSplitNumMap.get(samFile);
	            	 System.out.println("Running:"+samFile);  
	            	 //初始化
	            	 for (String circId : circFSJMap.keySet()) {
	            		 circFSJMap.put(circId, 0);
					 }
	            	 threadSub.reset();//更新子线程
	 			     threadMain.await();//所有线程激活2
	 			     threadMain.reset();//更新主线程
	 				 threadSub.await();//主线程关闭，子线程激活2   
	 				//取出第二遍扫描后的circRNA信息，并与第一遍的合并
	 				for (String circKey : circFSJMap.keySet()) {
	 					int temNum = circFSJMap.get(circKey);					
	 					int sunNum = circSumFSJMap.get(circKey);
	 					circSumFSJMap.put(circKey, sunNum+temNum);
	 					FSJmatrix[circRowMap.get(circKey)][j] = temNum;
	 				}				
	 			}			
				System.out.println(df.format(System.currentTimeMillis())+" "+":Second scan completed");  
				fileLog.write(df.format(System.currentTimeMillis())+" "+":Second scan completed"+"\n");
				//清理内存
				chrSiteMap1 = null;
				SiteMap1 = null;
				siteList1 = null;
				chrSiteMap2 = null;
				SiteMap2 = null;
				siteList2 = null;
				siteArrayMap1 = null;
				siteArrayMap2 = null;
				circFSJMap = null;
				//整理BSJMatrix
		        HashMap<String, Integer> circMap = new HashMap<String, Integer>();
		        for (int i = 0; i < filePathList.size(); i++) {
		        	circMap.clear();
					String samFile = filePathList.get(i);
					 AllFileSplitNum = fileSplitNumMap.get(samFile);
					for (int j = 1; j <= AllFileSplitNum; j++) {
						BufferedReader BSJBr = new BufferedReader(new FileReader(new File(samFile+"BSJ"+j)));	
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
					}					
					for (String circKey : circMap.keySet()) {
						BSJmatrix[circRowMap.get(circKey)][i] = circMap.get(circKey);
					}
				}
		        //合并总结	
				Summary summary = new Summary(strigency,chrTCGAMap);
				ArrayList<String> SummaryCircList = summary.summary(filePathList,fileSplitNumMap,circSumFSJMap,UserGivecircRNA);
				HashMap<String, String> circTrueIdMap = summary.getCircMap();
				summary = null;
				circSumFSJMap = null;
				chrTCGAMap = null;
				//输出每个样本CircRNA对应BSJnum
		        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
		        BSJCount.write("circRNA_ID");
		        for (String sampleName : fileNameList) {
		        	BSJCount.write("\t"+sampleName);
				}
		        BSJCount.write("\n");
		      //输出每个样本CircRNA对应FSJnum
		        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
		        FSJCount.write("circRNA_ID");
		        for (String sampleName : fileNameList) {
		        	FSJCount.write("\t"+sampleName);
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
	 			//进行第二遍扫描   	  
		        for (int j = 0; j < filePathList.size(); j++) {
					 samFile = filePathList.get(j);
					 matchNum = 0;
					 incr.set(1);
					 AllFileSplitNum = fileSplitNumMap.get(samFile);
	            	 System.out.println("Running:"+samFile);  
	            	 //初始化
	            	 for (String circId : circFSJMap.keySet()) {
	            		 circFSJMap.put(circId, 0);
	            		 circBSJMap.put(circId, 0);
					 }       
	            	 threadSub.reset();//更新子线程
	 			     threadMain.await();//所有线程激活2
	 			     threadMain.reset();//更新主线程
	 				 threadSub.await();//主线程关闭，子线程激活2 
	 				System.out.println(samFile+" Mapped_Reads"+" "+matchNum);  
					fileLog.write(samFile+" Mapped_Reads"+" "+matchNum+"\n");
	 				//取出第二遍扫描后的circRNA信息，并与第一遍的合并
	 				for (String circKey : circFSJMap.keySet()) {
	 					int temNumFSJ = circFSJMap.get(circKey);	
	 					int temNumBSJ = circBSJMap.get(circKey);
	 					FSJmatrix[circRowMap.get(circKey)][j] = temNumFSJ;
	 					BSJmatrix[circRowMap.get(circKey)][j] = temNumBSJ;
	 				}
	 				
	 			}
		      //输出每个样本CircRNA对应BSJnum
		        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
		        BSJCount.write("circRNA_ID");
		        for (String sampleName : fileNameList) {
		        	BSJCount.write("\t"+sampleName);
				}
		        BSJCount.write("\n");
		      //输出每个样本CircRNA对应FSJnum
		        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
		        FSJCount.write("circRNA_ID");
		        for (String sampleName : fileNameList) {
		        	FSJCount.write("\t"+sampleName);
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
	       
			threadMain.await();//所有线程激活2
		} catch (Exception e) { 
			e.printStackTrace();
		}
		
		poolExe.shutdown();
		//删除临时文件
 		for (int j = 0; j < filePathList.size(); j++) {
			String samFile = filePathList.get(j);
			//删除sam文件
			if (isSam){
				new File(samFile).delete();			
			}
			AllFileSplitNum = fileSplitNumMap.get(samFile);
			for (int i = 1; i <= AllFileSplitNum; i++) {
				new File(samFile+"BSJ"+i).delete();
			}		
 		}
 		
		long endTime = System.currentTimeMillis(); 
		System.out.println("Program run time:" + (endTime - startTime) + "ms");
		fileLog.write("Program run time:" + (endTime - startTime) + "ms"+"\n");
		fileLog.close();
		return true;
	}
}
