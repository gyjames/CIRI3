package com.zx.findcircrna;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;

import com.zx.MultiThreaded.FileReader;

public class MutFileFindCircRNAScan2 extends FileFindCircRNAScan2{
	public MutFileFindCircRNAScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min,
			byte[] siteArray1, byte[] siteArray2, HashMap<Integer, ArrayList<SiteSort>> chrSiteMap1,
			HashMap<Integer, ArrayList<SiteSort>> chrSiteMap2, HashMap<String, String> chrTCGAMap, int seqLen,
			boolean intronLable) throws IOException {
		super(minMapqUni, circFSJMap, linear_range_size_min, siteArray1, siteArray2, chrSiteMap1, chrSiteMap2, chrTCGAMap,
				seqLen, intronLable);
		// TODO Auto-generated constructor stub
	}

	public void findCircRNAScan2(String input,HashMap<String, String> scan1IdMap,int segment, int segmentNum) throws IOException {
		readNum = 0;
		boolean lable = false;
		FileInputStream fileIn = new FileInputStream(input);
		FileChannel fileChannel = fileIn.getChannel();
		long fileSize = fileChannel.size();
		// 判断从文件读取的起始位置
		long start = (long) (Math.floor((fileSize / segment) * (segmentNum - 1) / 1024.0) * 1024);
		long end = (long) (Math.floor((fileSize / segment) * (segmentNum) / 1024.0) * 1024);
		//////////////////
		String id = "", oneRead = "";
		HashMap<Integer, ArrayList<String>> standMap = new HashMap<Integer, ArrayList<String>>();
		ArrayList<String> serveInforList = new ArrayList<String>();
		standMap.put(0, serveInforList);
		standMap.put(1, serveInforList);
		HashMap<Integer, String> readsMap = new HashMap<Integer, String>();
				
		// 判断来自第几个线程
		if (segmentNum == 1) {
			FileReader fileReader = new FileReader(fileChannel, 1024, start);
			String line = fileReader.readline();
			while (fileChannel.position() <= end) {
				if (line.startsWith("@")) {
					line = fileReader.readline();
					continue;
				}
				String lineArr[] = line.split("\t");
				int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3]
						+ "\t" + lineArr[4] + "\t" + lineArr[5];

				if (!id.equals(lineArr[0])) {
					//判断是否匹配上
					if(lable) {
						readNum++;
					}
					lable = false;
					if(scan1IdMap.containsKey(id)) {
						
					}else {
						isBSJScan2.isCandidate(standMap, readsMap);
					}
					standMap = new HashMap<Integer, ArrayList<String>>();
					serveInforList = new ArrayList<String>();
					serveInforList.add(serveInfor);
					standMap.put(standKey, serveInforList);
					readsMap.clear();
					readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
					id = lineArr[0];
					oneRead = standKey+"";
				} else {
					if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
						readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
						serveInforList = new ArrayList<String>();
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
						oneRead = standKey+"";

					} else {
						serveInforList = standMap.get(standKey);
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
					}
				}
				if(!lineArr[5].equals("*")) {
					lable = true;
				}
				line = fileReader.readline();
			}
			//读取最后一个
			while (true) {			
				String lineArr[] = line.split("\t");
				int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
				//id:0 stand:1 chr:2 pos:3 mapq:4 cigar:5 
				String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3] + 
						"\t" + lineArr[4] + "\t" + lineArr[5];
				if(id.equals(lineArr[0])) {
					if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
						readsMap.put(standKey, stand.stand5(lineArr[1])+lineArr[9].toUpperCase());
						serveInforList = new ArrayList<String>();
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
						oneRead = standKey+"";
						
					}else {					
							serveInforList = standMap.get(standKey);
							serveInforList.add(serveInfor);
							standMap.put(standKey, serveInforList);					
					}
					if(!lineArr[5].equals("*")) {
						lable = true;
					}
					line = fileReader.readline();	
				}else {
					break;
				}
				
			}
			//最后一个
			//判断是否匹配上
			if(lable) {
				readNum++;
			}
			if(scan1IdMap.containsKey(id)) {
				
			}else {
				isBSJScan2.isCandidate(standMap, readsMap);
			}	
			 fileReader.close();
		}else if (segmentNum == segment) {
			FileReader fileReader = new FileReader(fileChannel, 1024, start-1025);
			String line = fileReader.readline();
			line = fileReader.readline();
			while (true) {
				if(fileChannel.position() <= start) {
					String[] lineArr = line.split("\t");
					id = lineArr[0];
					line = fileReader.readline();
				}else {
					String[] lineArr = line.split("\t");
					if(!id.equals(lineArr[0])) {
						break;
					}
					line = fileReader.readline();
				}
				
			}
			while (line != null) {
				String[] lineArr = line.split("\t");
				int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3]
						+ "\t" + lineArr[4] + "\t" + lineArr[5];
				if (!id.equals(lineArr[0])) {
					//判断是否匹配上
					if(lable) {
						readNum++;
					}
					lable = false;
					if(scan1IdMap.containsKey(id)) {
						
					}else {
						isBSJScan2.isCandidate(standMap, readsMap);
					}
					standMap = new HashMap<Integer, ArrayList<String>>();
					serveInforList = new ArrayList<String>();
					serveInforList.add(serveInfor);
					standMap.put(standKey, serveInforList);
					readsMap.clear();
					readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
					id = lineArr[0];
					oneRead = standKey+"";
				} else {
					if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
						readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
						serveInforList = new ArrayList<String>();
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
						oneRead = standKey+"";

					} else {
						serveInforList = standMap.get(standKey);
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
					}
				}
				if(!lineArr[5].equals("*")) {
					lable = true;
				}
				line = fileReader.readline();

			}
			//判断是否匹配上
			if(lable) {
				readNum++;
			}
			if(scan1IdMap.containsKey(id)) {
				
			}else {
				isBSJScan2.isCandidate(standMap, readsMap);
			}
			fileReader.close();
		}else {
			//中间线程
			FileReader fileReader = new FileReader(fileChannel, 1024, start-1025);
			String line = fileReader.readline();
			line = fileReader.readline();
			while (true) {
				if(fileChannel.position() <= start) {
					String[] lineArr = line.split("\t");
					id = lineArr[0];
					line = fileReader.readline();
				}else {
					String[] lineArr = line.split("\t");
					if(!id.equals(lineArr[0])) {
						break;
					}
					line = fileReader.readline();
				}
				
			}
			while (fileChannel.position() <= end) {
				String[] lineArr = line.split("\t");
				int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
				String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3]
						+ "\t" + lineArr[4] + "\t" + lineArr[5];

				if (!id.equals(lineArr[0])) {
					//判断是否匹配上
					if(lable) {
						readNum++;
					}
					lable = false;
					if(scan1IdMap.containsKey(id)) {
						
					}else {
						isBSJScan2.isCandidate(standMap, readsMap);
					}
					standMap = new HashMap<Integer, ArrayList<String>>();
					serveInforList = new ArrayList<String>();
					serveInforList.add(serveInfor);
					standMap.put(standKey, serveInforList);
					readsMap.clear();
					readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
					id = lineArr[0];
					oneRead = standKey+"";
				} else {
					if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
						readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
						serveInforList = new ArrayList<String>();
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
						oneRead = standKey+"";

					} else {
						serveInforList = standMap.get(standKey);
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
					}
				}
				if(!lineArr[5].equals("*")) {
					lable = true;
				}
				line = fileReader.readline();

			}
			//读取最后一个
			while (true) {			
				String lineArr[] = line.split("\t");
				int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
				//id:0 stand:1 chr:2 pos:3 mapq:4 cigar:5 
				String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3] + 
						"\t" + lineArr[4] + "\t" + lineArr[5];
				if(id.equals(lineArr[0])) {
					if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
						readsMap.put(standKey, stand.stand5(lineArr[1])+lineArr[9].toUpperCase());
						serveInforList = new ArrayList<String>();
						serveInforList.add(serveInfor);
						standMap.put(standKey, serveInforList);
						oneRead = standKey+"";
						
					}else {					
							serveInforList = standMap.get(standKey);
							serveInforList.add(serveInfor);
							standMap.put(standKey, serveInforList);					
					}
					if(!lineArr[5].equals("*")) {
						lable = true;
					}
					line = fileReader.readline();	
				}else {
					break;
				}
				
			}
			//最后一个
			if(lable) {
				readNum++;
			}
			if(scan1IdMap.containsKey(id)) {
				
			}else {
				isBSJScan2.isCandidate(standMap, readsMap);
			}	
			 fileReader.close();
		}

		
	}

}
