package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class UserFindCircRNASTARScan2 {
	IsBSJScan2Star isBSJScan2;
	IsStand stand = new IsStand();	
	long readNum = 0;
	HashMap<String, Integer> circBSJMap = new HashMap<String, Integer>();
	public UserFindCircRNASTARScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min, HashMap<String, byte[]> siteArrayMap1,
			HashMap<String, byte[]> siteArrayMap2, HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2, HashMap<String, String> chrTCGAMap, int seqLen) throws IOException {
		
		super();
		isBSJScan2 = new IsBSJScan2Star(minMapqUni,circFSJMap, linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1,chrSiteMap2,chrTCGAMap,seqLen);
		circBSJMap.putAll(circFSJMap);
	}	
    public HashMap<String, Integer> getCircFSJMap() throws IOException {
		return isBSJScan2.getCircFSJMap();
	}
    public void setFSJScan2List() {
    	isBSJScan2.setFSJScan2List();
	}
    public HashMap<String, Integer> getCircBSJMap(){
		return circBSJMap;
	}
    public void setBSJScan2List() {
		for (String circ : circBSJMap.keySet()) {
			circBSJMap.put(circ, 0);
		}
	}
    public long getReadNum(){
		return readNum;
	}
	public void setReadNum(){
		readNum = 0;
	}
	public void findCircRNAScan2(String input,HashMap<String, String> scan1IdMap) throws IOException {
		boolean matchLable = false;
		BufferedReader samReadScan2 = new BufferedReader(new FileReader(new File(input)));	
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		String id = "";
		int oneRead = 0;
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
		String line = samReadScan2.readLine();
		while (line != null) {
			if (line.startsWith("@")) {
				line = samReadScan2.readLine();
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
				//read
				readsMap.clear();
				serveInforList = new ArrayList<String[]>();
				serveInforList.add(serveInfor);
				readsMap.put(readKey, serveInforList);
				//stand
				standMap.clear();
				String[] alignBackArr = lineArr[6].split("\t",5);
				standMap.put(readKey, stand.stand5(lineArr[1])+alignBackArr[3].toUpperCase());
				//记录read的最大长度
				id = lineArr[0];	
				oneRead = readKey;
				}else {
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
			if(!lineArr[5].equals("*")) {
				matchLable = true;
			}
			line = samReadScan2.readLine();
		}
		samReadScan2.close();
		//判断是否匹配上
		if(matchLable) {
			readNum++;
		}
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
