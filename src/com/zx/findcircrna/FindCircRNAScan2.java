package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FindCircRNAScan2 {
	IsBSJScan2 isBSJScan2;
	IsStand stand = new IsStand();	
	public FindCircRNAScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min, HashMap<String, byte[]> siteArrayMap1,
			HashMap<String, byte[]> siteArrayMap2, HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2, HashMap<String, String> chrTCGAMap, int seqLen,boolean intronLable) throws IOException {
		
		super();
		if (intronLable) {
			isBSJScan2 = new IsBSJScan2Intron(minMapqUni,circFSJMap, linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1,chrSiteMap2,chrTCGAMap,seqLen);
		}else {
			isBSJScan2 = new IsBSJScan2(minMapqUni,circFSJMap, linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1,chrSiteMap2,chrTCGAMap,seqLen);
		}		
	}	
    public HashMap<String, Integer> getCircFSJMap() throws IOException {
		return isBSJScan2.getCircFSJMap();
	}
    public void setFSJScan2List() {
    	isBSJScan2.setFSJScan2List();
	}
	public void findCircRNAScan2(String input,HashMap<String, String> scan1IdMap) throws IOException {
		//存放第一遍扫描的circRNA id       
		BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(input+"BSJ1"),true));

		BufferedReader samReadScan2 = new BufferedReader(new FileReader(new File(input)));	
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		//readsMap.put(0, new ArrayList<String[]>());
		//readsMap.put(1, new ArrayList<String[]>());
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
			    //判断是否含有BSJ
				
                if(scan1IdMap.containsKey(id)) {
					
				}else {
					String circInfor = isBSJScan2.isCandidate(readsMap, standMap);
					if(circInfor != null) {
						BSJOut.write(id+"\t"+circInfor+"\n");
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

			line = samReadScan2.readLine();
		}
		samReadScan2.close();
		//判断是否含有BSJ
		if(scan1IdMap.containsKey(id)) {
			
		}else {
			String circInfor = isBSJScan2.isCandidate(readsMap, standMap);
			if(circInfor != null) {
				BSJOut.write(id+"\t"+circInfor+"\n");
			}
			
		}	
		BSJOut.close();
	}

}
