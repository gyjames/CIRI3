package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GetChimericOut {
	MisdStar misd = new MisdStar(); 
	int maxCircle,minCircle;
	HashMap<String, String> chrTCGAMap;
	String mitochondrion;
	boolean mlable,spLable;
	private String[][] bibasesMut = { { "AC", "AG", "GC", "AG", "AT", "AC", "AT", "AG"}, 
            { "CT", "GT", "CT", "GC", "GT", "AT", "CT", "AT"} };
    private String[][] bibases = { { "AC", "AG" }, { "CT", "GT" } };
    private String[] strandIndex = { "-", "+" };
    HashMap<String, String> idCircMap = new HashMap<String, String>();
	public GetChimericOut(int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,
			HashMap<String, String> chrExonStartMap,HashMap<String, String> chrExonEndMap,HashMap<String, String> chrTCGAMap,
			HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,HashMap<String, ArrayList<String>> chrExonEndTranscriptMap,String mitochondrion,
			boolean mlable,boolean spLable) {
		super();
		this.maxCircle = maxCircle;
		this.minCircle = minCircle;
		this.chrTCGAMap = chrTCGAMap;	
		this.mitochondrion = mitochondrion;	
		this.mlable = mlable;	
		this.spLable = spLable;	
	}
	public HashMap<String, String> getIdCircMap() {
		return idCircMap;
	}
	public void getBSJ(String ChimericPath) throws IOException {
		Pattern pattern = Pattern.compile("^((?:\\d+[A-Za-z])+)(?=-?\\d+p)");
		idCircMap = new HashMap<String, String>();
		// TODO Auto-generated method stub
		BufferedReader br = new BufferedReader(new FileReader(new File(ChimericPath)));
		String line = br.readLine();
		while (line != null) {	
			String[] arr = line.split("\t");
			//染色体相同，比对方向相同
			if( (!arr[0].equals(mitochondrion) || mlable) && arr[0].equals(arr[3]) && arr[2].equals(arr[5]) && 
					(arr[2].equals("-") && Integer.valueOf(arr[1]) < Integer.valueOf(arr[4]) || arr[2].equals("+") && Integer.valueOf(arr[1]) > Integer.valueOf(arr[4]))) {
			
				int startSite1 = Integer.valueOf(arr[10]);
				int startSite2 = Integer.valueOf(arr[12]);
				String CIGARStr1 =  arr[11];
				String CIGARStr2 =  arr[13];
				String alignmentStr1 = "",alignmentStr2 = "",anoAlignmentStr = "";
				String alignmentSite1 = "",alignmentSite2 = "",anoAlignmentSite = "";
				String alignCIGAR1 = "", alignCIGAR2 = "";
				if (CIGARStr1.indexOf("p") >= 0 && CIGARStr2.indexOf("p") >= 0) {
					
					
				}else if (CIGARStr1.indexOf("p") >= 0) {
					alignCIGAR2 = CIGARStr2;
					ArrayList<String> CIGARiteList = misd.misd(CIGARStr1,startSite1);
					if (arr[2].equals("+")) {
						String[] pCIGATArr = CIGARStr1.split("p");
						alignCIGAR1 = pCIGATArr[1];
						anoAlignmentStr = CIGARiteList.get(0);
						anoAlignmentSite = CIGARiteList.get(1);
						alignmentStr1 = CIGARiteList.get(2);
						alignmentSite1 = CIGARiteList.get(3);
					}else {
						Matcher matcher = pattern.matcher(CIGARStr1);
						if (matcher.find()) {
							alignCIGAR1 = matcher.group(1);  // 输出：69S31M
					     }
						alignmentStr1 = CIGARiteList.get(0);
						alignmentSite1 = CIGARiteList.get(1);
						anoAlignmentStr = CIGARiteList.get(2);
						anoAlignmentSite = CIGARiteList.get(3);
					}
					CIGARiteList = misd.misd(CIGARStr2,startSite2);
					alignmentStr2 = CIGARiteList.get(0);
					alignmentSite2 = CIGARiteList.get(1);
				}else if (CIGARStr2.indexOf("p") >= 0) {
					alignCIGAR1 = CIGARStr1;
					ArrayList<String>  CIGARiteList = misd.misd(CIGARStr1,startSite1);
					alignmentStr2 = CIGARiteList.get(0);
					alignmentSite2 = CIGARiteList.get(1);
					CIGARiteList = misd.misd(CIGARStr2,startSite2);
					if (arr[2].equals("+")) {
						Matcher matcher = pattern.matcher(CIGARStr2);
						if (matcher.find()) {
							alignCIGAR2 = matcher.group(1);  // 输出：69S31M
					    }
						alignmentStr1 = CIGARiteList.get(0);
						alignmentSite1 = CIGARiteList.get(1);
						anoAlignmentStr = CIGARiteList.get(2);
						anoAlignmentSite = CIGARiteList.get(3);
					}else {
						String[] pCIGATArr = CIGARStr2.split("p");
						alignCIGAR2 = pCIGATArr[1];
						anoAlignmentStr = CIGARiteList.get(0);
						anoAlignmentSite = CIGARiteList.get(1);
						alignmentStr1 = CIGARiteList.get(2);
						alignmentSite1 = CIGARiteList.get(3);
					}
				}
				if (!alignmentStr1.equals("")) {
					String[] alignmentArr1 = alignmentStr1.split("\t");
					String[] alignmentArr2 = alignmentStr2.split("\t");
					String[] alignmentSiteArr1 = alignmentSite1.split("\t");
					String[] alignmentSiteArr2 = alignmentSite2.split("\t");
					//SM在前 MS在后  
					if (Integer.valueOf(alignmentArr1[0]) < Integer.valueOf(alignmentArr2[0])) {
						String[] alignmentArr =  alignmentArr1;
						alignmentArr1 = alignmentArr2;
						alignmentArr2 = alignmentArr;
						
						String[] alignmentSiteArr = alignmentSiteArr1;
						alignmentSiteArr1 = alignmentSiteArr2;
						alignmentSiteArr2 = alignmentSiteArr;
					}
					int circStart = Integer.valueOf(alignmentSiteArr1[0]);
					int circEnd = Integer.valueOf(alignmentSiteArr2[1]);
					int cirScale = circEnd - circStart + 1;
					if (cirScale > 0 && Math.abs(Integer.valueOf(alignmentArr1[0]) - Integer.valueOf(alignmentArr2[0]) - Integer.valueOf(alignmentArr2[1]))<=6 
							&& cirScale <= maxCircle && cirScale >= minCircle) {
						String chrTCGAStr = chrTCGAMap.get(arr[0]); 
						String signal1 = chrTCGAStr.substring(circStart-3,circStart-1);
						String signal2 = chrTCGAStr.substring(circEnd,circEnd+2);
						if (spLable) {
							for (int i = 0; i <= 7; i++) {
								if (bibasesMut[0][i].equals(signal1) && bibasesMut[1][i].equals(signal2)) {
									idCircMap.put(arr[9], arr[9]+"\t"+alignCIGAR1+";"+alignCIGAR2+"\t"+1+"\t"+arr[0]+"\t"+circStart+"\t"+circEnd+"\t"+strandIndex[i % 2]+"\t"+
											signal1+"\t"+signal2+"\t"+1);
									break;
								}
							}
						}else {
							for (int i = 0; i <= 1; i++) {
								if (bibases[0][i].equals(signal1) && bibases[1][i].equals(signal2)) {
									idCircMap.put(arr[9], arr[9]+"\t"+alignCIGAR1+";"+alignCIGAR2+"\t"+1+"\t"+arr[0]+"\t"+circStart+"\t"+circEnd+"\t"+strandIndex[i % 2]+"\t"+
											signal1+"\t"+signal2+"\t"+1);
									break;
								}
							}
						}
					}
				}
			}
			line = br.readLine();
		}
		br.close();
	}

}
