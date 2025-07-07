package com.zx.findcircrna;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

public class BamUserFindCircRNASTARScan2 extends UserFindCircRNASTARScan2 {
	public BamUserFindCircRNASTARScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min,
			HashMap<String, byte[]> siteArrayMap1, HashMap<String, byte[]> siteArrayMap2,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2, HashMap<String, String> chrTCGAMap,
			int seqLen) throws IOException {
		super(minMapqUni, circFSJMap, linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1, chrSiteMap2, chrTCGAMap,
				seqLen);
		// TODO Auto-generated constructor stub
	}

	public void findCircRNAScan2(String input,HashMap<String, String> scan1IdMap) throws IOException {
		boolean matchLable = false;	
		//存放第一遍扫描的circRNA id       
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		String id = "";
		int oneRead = 0,readKey = 0,standKey = 0;
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
		// Open BAM file for reading
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(new File(input));
		CloseableIterator<SAMRecord> iterator = reader.iterator();		
		while (iterator.hasNext()) {
	    	SAMRecord record = iterator.next();
	    	String temid =record.getReadName();//simulate:1	        
	        String chr = record.getContig();//chr1
	        int start = record.getStart();//188535796
	        int MQ = record.getMappingQuality();//60
	        String cigar = record.getCigarString();//63H37M  
			if(record.getFirstOfPairFlag()) {
				readKey = 1;
	        }else {
	        	readKey = 0;
	        }	
			if(record.getReadNegativeStrandFlag()) {
				standKey = 1;
	        }else {
	        	standKey = 0;
	        }
			if(chr == null) {
				chr = "*";
			}
			String[] serveInfor = {standKey+"",chr,start+"",MQ+"",cigar};		
			if (!id.equals(temid)) {
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
				String seq = record.getReadString().toUpperCase();
				standMap.put(readKey, standKey+seq);
				//记录read的最大长度
				id = temid;	
				oneRead = readKey;
				}else {
					if (readKey != oneRead) {
						String seq = record.getReadString().toUpperCase();
						standMap.put(readKey, standKey+seq);
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
			if(!cigar.equals("*")) {
				matchLable = true;
			}
		}
		reader.close();
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
