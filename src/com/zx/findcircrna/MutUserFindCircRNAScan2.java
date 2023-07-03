package com.zx.findcircrna;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;

import com.zx.MultiThreaded.FileReader;

public class MutUserFindCircRNAScan2 extends UserFindCircRNAScan2{
	public MutUserFindCircRNAScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min,
			HashMap<String, byte[]> siteArrayMap1, HashMap<String, byte[]> siteArrayMap2,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2, HashMap<String, String> chrTCGAMap,
			int seqLen) throws IOException {
		super(minMapqUni, circFSJMap, linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1, chrSiteMap2, chrTCGAMap,
				seqLen);
		// TODO Auto-generated constructor stub
	}

	public void findCircRNAScan2(String samFile,HashMap<String, String> scan1IdMap,int threads, int threadNum) throws IOException {
		boolean matchLable = false;
		//存放第一遍扫描的circRNA id       
		FileInputStream fileIn = new FileInputStream(samFile);
		FileChannel fileChannel = fileIn.getChannel();
		long fileSize = fileChannel.size();
		// 判断从文件读取的起始位置
		long fileStart,fileEnd;
		if (threadNum == 1) {
			fileStart = 0;
			fileEnd = (long) (Math.floor((fileSize / threads) * (threadNum) / 1024.0) * 1024);		
		}else {
			fileStart = (long) (Math.floor((fileSize / threads) * (threadNum - 1) / 1024.0) * 1024)-1025;
			fileEnd = (long) (Math.floor((fileSize / threads) * (threadNum) / 1024.0) * 1024);		
		}
		//////////////////
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		//readsMap.put(0, new ArrayList<String[]>());
		//readsMap.put(1, new ArrayList<String[]>());
		String id = "",line = "";
		int oneRead = 0;
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
		//standMap.put(0, "");
		//standMap.put(1, "");
		//读取数据		
		FileReader fileReader = new FileReader(fileChannel, 1024, fileStart);
		line = fileReader.readline();
		// 判断来自第几个线程
		if (threadNum != 1) {			
			line = fileReader.readline();
			while (true) {
				if(fileChannel.position() <= fileStart + 1025) {
					String[] lineArr = line.split("\t",2);
					id = lineArr[0];
					line = fileReader.readline();
				}else {
					String[] lineArr = line.split("\t",2);
					if(!id.equals(lineArr[0])) {
						break;
					}
					line = fileReader.readline();
				}
				
			}
		}				
		if (threadNum == threads) {
			while (line != null) {
				String[] lineArr = line.split("\t",7);
				int readKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String[] serveInfor = {stand.stand5(lineArr[1]),lineArr[2],lineArr[3],lineArr[4],lineArr[5]};
				if (!id.equals(lineArr[0])) {
					//判断是否匹配上
					if(matchLable) {
						readNum++;
					}
					matchLable = false;
				    //判断是否含有BSJ
					
	                if(scan1IdMap.containsKey(id)) {
						
					}else {
						String circInfor = isBSJScan2.isCandidate(readsMap, standMap);
						if(circInfor != null) {
							String[] circLine = circInfor.split("\t");
							if(circLine[1].equals("1")) {
								String circKey = circLine[2]+"\t"+circLine[3]+"\t"+circLine[4];
								int bsjNum = circBSJMap.get(circKey);
								circBSJMap.put(circKey, bsjNum+1);
							}
						
						}
						
					}		
					//read
	                readsMap.clear();
					serveInforList = new ArrayList<String[]>();
					serveInforList.add(serveInfor);
					readsMap.put(readKey, serveInforList);
					//stand
					standMap.clear();
					String[] alignBackArr = lineArr[6].split("\t",5);
					standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
					id = lineArr[0];	
					oneRead = readKey;
					} else {
						if (readKey != oneRead) {
							String[] alignBackArr = lineArr[6].split("\t",5);
							standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
							serveInforList = new ArrayList<String[]>();
							serveInforList.add(serveInfor);
							readsMap.put(readKey, serveInforList);
							oneRead = readKey;						
						}else {	
							serveInforList = readsMap.get(readKey);
							serveInforList.add(serveInfor);
							readsMap.put(readKey, serveInforList);					
						}
					}

				line = fileReader.readline();
			}
		}else {
			while (fileChannel.position() <= fileEnd) {
				if (line.startsWith("@")) {
					line = fileReader.readline();
					continue;
				}
				String[] lineArr = line.split("\t",7);
				int readKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String[] serveInfor = {stand.stand5(lineArr[1]),lineArr[2],lineArr[3],lineArr[4],lineArr[5]};
				if (!id.equals(lineArr[0])) {
					//判断是否匹配上
					if(matchLable) {
						readNum++;
					}
					matchLable = false;
				    //判断是否含有BSJ
					
	                if(scan1IdMap.containsKey(id)) {
						
					}else {
						String circInfor = isBSJScan2.isCandidate(readsMap, standMap);
						if(circInfor != null) {
							String[] circLine = circInfor.split("\t");
							if(circLine[1].equals("1")) {
								String circKey = circLine[2]+"\t"+circLine[3]+"\t"+circLine[4];
								int bsjNum = circBSJMap.get(circKey);
								circBSJMap.put(circKey, bsjNum+1);
							}
						
						}
						
					}		
					//清空
	                readsMap.clear();
					serveInforList = new ArrayList<String[]>();
					serveInforList.add(serveInfor);
					readsMap.put(readKey, serveInforList);
					//stand
					standMap.clear();
					String[] alignBackArr = lineArr[6].split("\t",5);
					standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
					id = lineArr[0];	
					oneRead = readKey;
					} else {
						if (readKey != oneRead) {
							String[] alignBackArr = lineArr[6].split("\t",5);
							standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
							serveInforList = new ArrayList<String[]>();
							serveInforList.add(serveInfor);
							readsMap.put(readKey, serveInforList);
							oneRead = readKey;						
						}else {	
							serveInforList = readsMap.get(readKey);
							serveInforList.add(serveInfor);
							readsMap.put(readKey, serveInforList);					
						}
					}
				line = fileReader.readline();
			}
			//读取最后一个
			while (true) {			
				String[] lineArr = line.split("\t",7);
				int readKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String[] serveInfor = {stand.stand5(lineArr[1]),lineArr[2],lineArr[3],lineArr[4],lineArr[5]};
				if(id.equals(lineArr[0])) {
					if (readKey != oneRead) {
						String[] alignBackArr = lineArr[6].split("\t",5);
						standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
						serveInforList = new ArrayList<String[]>();
						serveInforList.add(serveInfor);
						readsMap.put(readKey, serveInforList);
						oneRead = readKey;						
					}else {	
						serveInforList = readsMap.get(readKey);
						serveInforList.add(serveInfor);
						readsMap.put(readKey, serveInforList);					
					}						
				}else {
					break;
				}
				line = fileReader.readline();
			}
		}
		fileReader.close();
		///判断是否匹配上
		if(matchLable) {
			readNum++;
		}
		matchLable = false;
	    //判断是否含有BSJ
		
        if(scan1IdMap.containsKey(id)) {
			
		}else {
			String circInfor = isBSJScan2.isCandidate(readsMap, standMap);
			if(circInfor != null) {
				String[] circLine = circInfor.split("\t");
				if(circLine[1].equals("1")) {
					String circKey = circLine[2]+"\t"+circLine[3]+"\t"+circLine[4];
					int bsjNum = circBSJMap.get(circKey);
					circBSJMap.put(circKey, bsjNum+1);
				}
			
			}
			
		}		
		
	}

}
