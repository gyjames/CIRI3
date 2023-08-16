package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class ReadFaFile {
	
	HashMap<String, String> chrTCGAMap = new HashMap<String, String>();	
	HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
	ArrayList<String> hgList = new ArrayList<String>();	
	String chrName = "",chrTAGA;
	public void readFa(String faFile) throws IOException {
		BufferedReader hg = new BufferedReader(new FileReader(new File(faFile)));
		String HgLine = hg.readLine();
		while (HgLine != null) {
			if (HgLine.startsWith(">")) {		
			    chrTAGA = String.join("", hgList).toUpperCase();
			    chrLenMap.put(chrName, chrTAGA.length());
			    chrTCGAMap.put(chrName, chrTAGA);			    
			    chrName = HgLine.split("\t")[0].split(" ")[0].replace(">", "");
			    hgList.clear();
			    HgLine = hg.readLine();
				continue;
			}
			hgList.add(HgLine);
			HgLine = hg.readLine();
		}
		hg.close();
		chrTAGA = String.join("", hgList).toUpperCase();
		chrLenMap.put(chrName, chrTAGA.length());
		chrTCGAMap.put(chrName, chrTAGA);
	}
	public HashMap<String, String> getChrTCGAMap() {
		return chrTCGAMap;
	}

	public HashMap<String, Integer> getChrLenMap() {
		return chrLenMap;
	}
	
	
}
