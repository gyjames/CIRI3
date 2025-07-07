package com.zx.findcircrna;

import java.util.ArrayList;

public class MisdStar {

	public ArrayList<String> misd(String CIGAR,int faStart) {
		//以faSite为准
        ArrayList<String> CIGARiteList = new ArrayList<String>(); 
		String[] classNumStr = CIGAR.split("M|S|I|D|H|p|N");
		String classCIGAR = CIGAR.replaceAll("[^(a-zA-Z)]", "");
		int faSite = faStart-1;
		int readM = 0;
		int S1 = 0,S2 = 0;
		for (int i = 0; i <= classCIGAR.length()-1; i++) {
	    	int gap = Integer.valueOf(classNumStr[i]);
	    	if (classCIGAR.charAt(i) == 'M') {
	    		faSite = faSite + gap;
	    		readM = readM + gap;
			}else if (classCIGAR.charAt(i) == 'S') {
				if (i==0) {
					S1 = gap;
				}else {
					S2 = gap;
				}
			}else if (classCIGAR.charAt(i) == 'D') {
				faSite = faSite + gap;
				readM = readM + gap;
			}else if (classCIGAR.charAt(i) == 'N') {
				faSite = faSite + gap;
			}else if (classCIGAR.charAt(i) == 'I') {
				
			}else if (classCIGAR.charAt(i) == 'p') {
				CIGARiteList.add(S1+"\t"+readM+"\t"+S2);
				CIGARiteList.add(faStart+"\t"+faSite);
				//更新
				faSite = faSite + gap; 
				faStart = faSite+1;
				S1 = 0;
				readM = 0;	
			}
		}
		CIGARiteList.add(S1+"\t"+readM+"\t"+S2);
		CIGARiteList.add(faStart+"\t"+faSite);
		return CIGARiteList;
	}
	
	public ArrayList<String> misdScan2(String CIGAR,int faStart) {
		//以faSite为准
        ArrayList<String> CIGARiteList = new ArrayList<String>(); 
        ArrayList<String> FaiteList = new ArrayList<String>(); 
		String[] classNumStr = CIGAR.split("M|S|I|D|H|p|N");
		String classCIGAR = CIGAR.replaceAll("[^(a-zA-Z)]", "");
		int faSite = faStart-1;
		int readM = 0;
		int S1 = 0,S2 = 0;
		for (int i = 0; i <= classCIGAR.length()-1; i++) {
	    	int gap = Integer.valueOf(classNumStr[i]);
	    	if (classCIGAR.charAt(i) == 'M') {
	    		faSite = faSite + gap;
	    		readM = readM + gap;
			}else if (classCIGAR.charAt(i) == 'S') {
				if (i==0) {
					S1 = gap;
				}else {
					S2 = gap;
				}
			}else if (classCIGAR.charAt(i) == 'D') {
				faSite = faSite + gap;
				readM = readM + gap;
			}else if (classCIGAR.charAt(i) == 'N') {
				FaiteList.add(faStart+"-"+faSite);
				//更新
				faSite = faSite + gap; 
				faStart = faSite+1;
			}
		}
		FaiteList.add(faStart+"-"+faSite);
		CIGARiteList.add(S1+"\t"+readM+"\t"+S2);
		CIGARiteList.add(String.join("\t", FaiteList));
		return CIGARiteList;
	}
	 /*public static void main(String[] args) {
		 ArrayList<String> CIGARiteList = misd("37S63M-63p98M",148038789);
		 for (String aa : CIGARiteList) {
			System.out.println(aa);
		}
	}*/
}
