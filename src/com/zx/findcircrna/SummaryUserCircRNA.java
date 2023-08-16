package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class SummaryUserCircRNA {
	HashSet<String> circSiteSet = new HashSet<String>();	
	ArrayList<String>  circScan1List = new ArrayList<String>();
	private String[][] bibasesMut = { { "AC", "AG", "GC", "AG", "AT", "AC", "AT", "AG", "AC"}, 
            { "CT", "GT", "CT", "GC", "GT", "AT", "GT", "AT", "GT"} };
    private String[][] bibases = { { "AC", "AG" }, { "CT", "GT" } };
    private String[] strandIndex = { "-", "+" };
    public HashSet<String> getCircSiteSet() {
		return circSiteSet;
	}	
    public ArrayList<String> getCircScan1List() {
		return circScan1List;
	}
	public void summaryUserCircRNA(String UserGivecircRNA, boolean intronLable,HashMap<String, String> chrTCGAMap,String mitochondrion,boolean mlable) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(UserGivecircRNA)));
		String line = br.readLine();
		if (intronLable) {
			while (line != null ) {
				String[] arr = line.split("\t");
				int start = Integer.valueOf(arr[1]);
				int end = Integer.valueOf(arr[2]);
				String TCGASeq = chrTCGAMap.get(arr[0]);
				String startSeq,endSeq;
				if(intronLable) {
					startSeq = TCGASeq.substring(start,start+2);
					endSeq = TCGASeq.substring(end-3,end-1);
				}else {
					startSeq = TCGASeq.substring(start-3,start-1);
					endSeq = TCGASeq.substring(end,end+2);
				}				
				String stand = "NA";
				if(mlable && arr[0].equals(mitochondrion)) {
					for (int i = 0; i <= 8; i++) {
						if (startSeq.equals(bibasesMut[0][i]) && endSeq.equals(bibasesMut[1][i]) ) {
							stand = strandIndex[i%2];
						}
					}
				}else {
					for (int i = 0; i <= 1; i++) {
						if (startSeq.equals(bibases[0][i]) && endSeq.equals(bibases[1][i]) ) {
							stand = strandIndex[i];
						}
					}
				}
				circSiteSet.add(arr[0]+"\t"+arr[1]+"\t"+arr[2]+"\t"+stand+"\t"+startSeq+"\t"+endSeq+"\t"+arr[3]);
				circScan1List.add(line+"\t"+stand);
				line = br.readLine();
			}
		}else {
			while (line != null ) {
				String[] arr = line.split("\t");
				int start = Integer.valueOf(arr[1]);
				int end = Integer.valueOf(arr[2]);
				String TCGASeq = chrTCGAMap.get(arr[0]);
				String startSeq = TCGASeq.substring(start-3,start-1);
				String endSeq = TCGASeq.substring(end,end+2);
				String stand = "NA";
				if(mlable && arr[0].equals(mitochondrion)) {
					for (int i = 0; i <= 8; i++) {
						if (startSeq.equals(bibasesMut[0][i]) && endSeq.equals(bibasesMut[1][i]) ) {
							stand = strandIndex[i%2];
						}
					}
				}else {
					for (int i = 0; i <= 1; i++) {
						if (startSeq.equals(bibases[0][i]) && endSeq.equals(bibases[1][i]) ) {
							stand = strandIndex[i];
						}
					}
				}
				circSiteSet.add(line+"\t"+stand+"\t"+startSeq+"\t"+endSeq);
				circScan1List.add(line+"\t"+stand);
				line = br.readLine();
			}
		}
		
	}
}
