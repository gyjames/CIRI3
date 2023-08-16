package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FindCircRNAScan1 {
	IsBSJScan1 isBSJScan1;
	IsStand stand = new IsStand();
	int readLen = 0;
	long readNum = 0;
	public FindCircRNAScan1(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable
			,HashMap<String, String> chrExonStartMap,HashMap<String, String> chrExonEndMap,HashMap<String, String> chrTCGAMap,
			HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,HashMap<String, ArrayList<String>> chrExonEndTranscriptMap,String mitochondrion,
			boolean mlable) {
		super();
		isBSJScan1 = new IsBSJScan1(minMapqUni,maxCircle,minCircle,linear_range_size_min,intronLable
				,chrExonStartMap,chrExonEndMap,chrTCGAMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,mlable);
		
	}		
	public int getReadLen(){
		return readLen;
	}
	public long getReadNum(){
		return readNum;
	}
	public void setReadNum(){
		readNum = 0;
	}
	public void findCircRNAScan1(String samFile) throws IOException{	
		BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(samFile+"BSJ1")));		
		boolean matchLable = false;
		BufferedReader samReadScan1 = new BufferedReader(new FileReader(new File(samFile)));	
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		String id = "";
		int oneRead = 0,alignNum = 0,seqLen = 0;
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
		String line = samReadScan1.readLine();
		while (line != null) {
			if (line.startsWith("@")) {
				line = samReadScan1.readLine();
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
				if(alignNum > 2 || readKey == 1) {
					String circInfor = isBSJScan1.isBSJScan1(readsMap, standMap);
					if(circInfor != null) {
						BSJOut.write(id+"\t"+circInfor+"\n");
					}
							
				}
				alignNum = 0;
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
				seqLen = alignBackArr[3].length();
				if (seqLen > readLen) {
					readLen = seqLen;
				}
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

			if(!lineArr[5].equals("*")) {
				matchLable = true;
			}
			alignNum++;
			line = samReadScan1.readLine();
		}
		samReadScan1.close();
		//判断是否匹配上
		if(matchLable) {
			readNum++;
		}
		//判断是否含有BSJ		
		String circInfor = isBSJScan1.isBSJScan1(readsMap, standMap);
		if(circInfor != null) {
			BSJOut.write(id+"\t"+circInfor+"\n");
		}		
		BSJOut.close();
		
	}

	
}