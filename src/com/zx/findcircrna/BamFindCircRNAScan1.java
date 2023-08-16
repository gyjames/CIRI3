package com.zx.findcircrna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

public class BamFindCircRNAScan1  extends FindCircRNAScan1 {
	
	
	public BamFindCircRNAScan1(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
			boolean intronLable, HashMap<String, String> chrExonStartMap, HashMap<String, String> chrExonEndMap,
			HashMap<String, String> chrTCGAMap, HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,
			HashMap<String, ArrayList<String>> chrExonEndTranscriptMap, String mitochondrion, boolean mlable) {
		super(minMapqUni, maxCircle, minCircle, linear_range_size_min, intronLable, chrExonStartMap, chrExonEndMap, chrTCGAMap,
				chrExonStartTranscriptMap, chrExonEndTranscriptMap, mitochondrion, mlable);
		// TODO Auto-generated constructor stub
	}

	public void findCircRNAScan1(String samFile) throws IOException{		
		BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(samFile+"BSJ1")));		
		boolean matchLable = false;	
		HashMap<Integer, ArrayList<String[]>> readsMap = new HashMap<Integer, ArrayList<String[]>>();
		ArrayList<String[]> serveInforList;
		String id = "";
		int oneRead = 0,alignNum = 0,seqLen = 0,readKey = 0,standKey = 0;	
		HashMap<Integer, String> standMap = new HashMap<Integer, String>();
		// Open BAM file for reading
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(new File(samFile));		
		CloseableIterator<SAMRecord> iterator = reader.iterator();		
	    while (iterator.hasNext()) {
	    	SAMRecord record = iterator.next();
	    	String temid =record.getReadName();//simulate:1	        
	        String chr = record.getContig();//chr1
	        int start = record.getStart();//188535796
	        int MQ = record.getAlignmentStart();//60
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
	        String[] serveInfor = {standKey+"",chr,start+"",MQ+"",cigar};	
	        if (!id.equals(temid)) {
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
				String seq = record.getReadString().toUpperCase();
				standMap.put(readKey, standKey+seq);
				//记录read的最大长度
				seqLen = seq.length();
				if (seqLen > readLen) {
					readLen = seqLen;
				}
				id = temid;	
				oneRead = readKey;
				} else {
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
			alignNum++;
		}
	    reader.close();
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