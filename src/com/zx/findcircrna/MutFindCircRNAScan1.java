package com.zx.findcircrna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;

import com.zx.MultiThreaded.FileReader;

public class MutFindCircRNAScan1 extends FindCircRNAScan1{
	

	public MutFindCircRNAScan1(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
			boolean intronLable, HashMap<String, String> chrExonStartMap, HashMap<String, String> chrExonEndMap,
			HashMap<String, String> chrTCGAMap, HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,
			HashMap<String, ArrayList<String>> chrExonEndTranscriptMap, String mitochondrion, boolean mlable) {
		super(minMapqUni, maxCircle, minCircle, linear_range_size_min, intronLable, chrExonStartMap, chrExonEndMap, chrTCGAMap,
				chrExonStartTranscriptMap, chrExonEndTranscriptMap, mitochondrion, mlable);
		// TODO Auto-generated constructor stub
	}

	
	public void findCircRNAScan1(String samFile, int threads, int threadNum) throws IOException {		
		BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(samFile+"BSJ"+threadNum)));
		boolean matchLable = false;
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
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		//readsMap.put(0, new ArrayList<String[]>());
		//readsMap.put(1, new ArrayList<String[]>());
		String id = "",line = "";
		int oneRead = 0,alignNum = 0,seqLen = 0;
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
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
					if(alignNum > 2 || readKey == 1) {
						String circInfor = isBSJScan1.isBSJScan1(readsMap, standMap);
						if(circInfor != null) {
							BSJOut.write(id+"\t"+circInfor+"\n");
						}							
					}
					alignNum = 0;
					//清空
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
				if(!lineArr[5].equals("*")) {
					matchLable = true;
				}
				line = fileReader.readline();
			}
		}
		fileReader.close();
		//最后一个
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