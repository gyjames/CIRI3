package com.zx.findcircrna;

public class Misd {
	
	
	
	
public int[] misd(String oldCIGAR,int seqLength) {
	int[] CIGARite = {0,0,0,0};
	String CIGAR = oldCIGAR.replaceAll("H", "S");
	String[] classNumStr = CIGAR.split("M|S|I|D|H");
	String classCIGAR = CIGAR.replaceAll("[^(a-zA-Z)]", "");
	if (classNumStr.length==1) {
		if (classCIGAR.equalsIgnoreCase("M")) {
			CIGARite[3]=  seqLength;//[0, 0, 0, seqLength];
		}else {
			CIGARite[3]=  -1;//[0, 0, 0, -1];
		}
		
	}else {		
		if (classNumStr.length==2) {
			if (classCIGAR.equalsIgnoreCase("MS")) {
				CIGARite[0]=  1;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  Integer.valueOf(classNumStr[0])-1;
				CIGARite[3]=  Integer.valueOf(classNumStr[0]);//[1, $counts[0], $counts[0]-1, $counts[0]];
			}else if (classCIGAR.equalsIgnoreCase("SM")) {
				CIGARite[0]=  -1;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[3]=  Integer.valueOf(classNumStr[1]);//[-1, $counts[0], 0, $counts[1]];
			}else {
				CIGARite[3]=  -2;//[0, undef, undef, -2];
			}
		}else if (classNumStr.length==3) {
			if (classCIGAR.substring(0,1).equalsIgnoreCase("S") && classCIGAR.substring(2,3).equalsIgnoreCase("S")) {
				CIGARite[0]=  10;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  Integer.valueOf(classNumStr[2]);
				CIGARite[3]=  Integer.valueOf(classNumStr[1]);//[10, $counts[0], $counts[2], $counts[1]];
			} else if (classCIGAR.substring(0,1).equalsIgnoreCase("M") && classCIGAR.substring(1,2).equalsIgnoreCase("D") && classCIGAR.substring(2,3).equalsIgnoreCase("M")) {
				CIGARite[3]=  seqLength+Integer.valueOf(classNumStr[1]);//[0, 0, 0, $read_length+$counts[1]];
			} else if (classCIGAR.substring(0,1).equalsIgnoreCase("M") && classCIGAR.substring(1,2).equalsIgnoreCase("I") && classCIGAR.substring(2,3).equalsIgnoreCase("M")) {
				CIGARite[3]=  seqLength-Integer.valueOf(classNumStr[1]);//[0, 0, 0, $read_length-$counts[1]];
			} else {
				CIGARite[3]=  -2;//[0, undef, undef, -2];
			}
		}else if (classCIGAR.substring(0,1).equalsIgnoreCase("M") && classCIGAR.substring(classCIGAR.length()-1).equalsIgnoreCase("S")) {
			int M_sum = 0,D_sum = 0;
			for (int i = 0; i <= classNumStr.length-1; i++) {
				if (classCIGAR.substring(i,i+1).equalsIgnoreCase("M")) {
					M_sum += Integer.valueOf(classNumStr[i]);
				}else if (classCIGAR.substring(i,i+1).equalsIgnoreCase("D")) {
					D_sum += Integer.valueOf(classNumStr[i]);
				}
			}
			if (D_sum==0) {
				CIGARite[0]=  1;
				CIGARite[1]=  seqLength-Integer.valueOf(classNumStr[classCIGAR.length()-1]);
				CIGARite[2]=  M_sum -1;
				CIGARite[3]=  M_sum;//[1, $read_length-$counts[-1], $M_sum-1, $M_sum];
			} else {
				CIGARite[0]=  1;
				CIGARite[1]=  seqLength-Integer.valueOf(classNumStr[classCIGAR.length()-1]);
				CIGARite[2]=  M_sum +D_sum-1;
				CIGARite[3]=  M_sum+D_sum;//[1, $read_length-$counts[-1], $M_sum+$D_sum-1, $M_sum+$D_sum];
			}
			
		}else if (classCIGAR.substring(0,1).equalsIgnoreCase("S") && classCIGAR.substring(classCIGAR.length()-1).equalsIgnoreCase("M")) {
			int M_sum = 0,D_sum = 0;
			for (int i = 1; i < classNumStr.length; i++) {
				if (classCIGAR.substring(i,i+1).equalsIgnoreCase("M")) {
					M_sum += Integer.valueOf(classNumStr[i]);
				}else if (classCIGAR.substring(i,i+1).equalsIgnoreCase("D")) {
					D_sum += Integer.valueOf(classNumStr[i]);
				}
			}
			if (D_sum==0) {
				CIGARite[0]=  -1;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  0;
				CIGARite[3]=  M_sum;//[-1, $counts[0], 0, $M_sum];
			} else {
				CIGARite[0]=  -1;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  0;
				CIGARite[3]=  M_sum+D_sum;//[-1, $counts[0], 0, $M_sum+$D_sum];
			}
		}else if (classCIGAR.substring(0,1).equalsIgnoreCase("M") && classCIGAR.substring(classCIGAR.length()-1).equalsIgnoreCase("M")) {
			int M_sum = 0,D_sum = 0;
			for (int i = 0; i < classNumStr.length; i++) {
				if (classCIGAR.substring(i,i+1).equalsIgnoreCase("M")) {
					M_sum += Integer.valueOf(classNumStr[i]);
				}else if (classCIGAR.substring(i,i+1).equalsIgnoreCase("D")) {
					D_sum += Integer.valueOf(classNumStr[i]);
				}
			}
			if (D_sum==0) {
				CIGARite[3]=  M_sum;//[0, 0, 0, $M_sum];
			} else {
				CIGARite[3]=  M_sum+D_sum;//[0, 0, 0, $M_sum+$D_sum];
			}
		} else if ( classCIGAR.substring(0,1).equalsIgnoreCase("S") && classCIGAR.substring(classCIGAR.length()-1).equalsIgnoreCase("S")) {
			int M_sum = 0,D_sum = 0;
			for (int i = 1; i < classNumStr.length-1; i++) {
				if (classCIGAR.substring(i,i+1).equalsIgnoreCase("M")) {
					M_sum += Integer.valueOf(classNumStr[i]);
				}else if (classCIGAR.substring(i,i+1).equalsIgnoreCase("D")) {
					D_sum += Integer.valueOf(classNumStr[i]);
				}
			}
			if (D_sum==0) {
				CIGARite[0]=  10;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  Integer.valueOf(classNumStr[classCIGAR.length()-1]);
				CIGARite[3]=  M_sum;//[10, $counts[0], $counts[-1], $M_sum];
			} else {
				CIGARite[0]=  10;
				CIGARite[1]=  Integer.valueOf(classNumStr[0]);
				CIGARite[2]=  Integer.valueOf(classNumStr[classCIGAR.length()-1]);
				CIGARite[3]=  M_sum+D_sum;//[10, $counts[0], $counts[-1], $M_sum+$D_sum];
			}
		} else {
			CIGARite[3]=  -2;//[0, undef, undef, -2];
		} 
	}
	return CIGARite;
}
}
