package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class GetUserCircRNA {
	HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
    
	public HashMap<String, HashSet<String>> summaryUserCircRNA(String UserGivecircRNA, HashMap<String, String> chrTCGAMap) throws IOException {
		HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSiteSet = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(UserGivecircRNA)));
		String line = br.readLine();		
		while (line != null ) {
			String[] arr = line.split("\t");
			int start = Integer.valueOf(arr[1]);
			int end = Integer.valueOf(arr[2]);
			String TCGASeq = chrTCGAMap.get(arr[0]);
			String startSeq = TCGASeq.substring(start-3,start-1);
			String endSeq = TCGASeq.substring(end,end+2);
			String infor = arr[1]+"\t"+arr[2]+"\t"+"NA"+"\t"+startSeq+"\t"+endSeq;
			if(!chrCircSiteMap.containsKey(arr[0])) {
				circSiteSet = new HashSet<String>();
				circSiteSet.add(infor);
				chrCircSiteMap.put(arr[0], circSiteSet);
			}else {
				circSiteSet=chrCircSiteMap.get(arr[0]);
				circSiteSet.add(infor);
				chrCircSiteMap.put(arr[0], circSiteSet);
			}
			line = br.readLine();		
	  }
	br.close();
	return chrCircSiteMap;
	}
}
