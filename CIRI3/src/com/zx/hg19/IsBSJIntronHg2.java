package com.zx.hg19;

public class IsBSJIntronHg2 extends IsBSJHg2{
	int initial_size2 = 5;
	public IsBSJIntronHg2(int linear_range_size_min,int minMapqUni) {
		super(linear_range_size_min, minMapqUni);
	}
	public String isBSJHg2Intron(String circLineArr[], String chrTAGA) {
	    //Variables are assigned according to inputted variables
	    // Existence of str3 can be determined by number of inputted variables	
	    String str,pem_null_range_seq = "",initial_seq,circ_range_seq;
		int startSite = Integer.valueOf(circLineArr[3]);
		int endSite = Integer.valueOf(circLineArr[4]);
		int quant = Integer.valueOf(circLineArr[12]);
		circ_range_seq = chrTAGA.substring(startSite-1,endSite);			
		int len_str = circLineArr[5].length();
		int circRangeLen = circ_range_seq.length();
		str = circLineArr[5];
		//Alignment style: xS/HyM or xS/HyMzS/H
		if (circLineArr[2].equals("sm")) {			
			if(str.length()<5) {
				if (!circ_range_seq.substring(circRangeLen-len_str,circRangeLen).equals(str)) {
					return "0";
				}else {
					return "2";	
				}
			}else {
				int Initial = 1;
				initial_seq = str.substring(len_str-initial_size2,len_str);
				if (circ_range_seq.substring(circRangeLen-initial_size2,circRangeLen).equals(initial_seq)) {
					Initial = 0;
				}else {
					aligner.setSeq(circ_range_seq.substring(circRangeLen-initial_size2,circRangeLen),initial_seq);
					String[] alignment = aligner.getAlignment();
					if (aligner.getAlignmentScore()>=3 && alignment[0].length()>=4) {
						
					}else {
						return "0";
					}
				}
				boolean lable = true;
				String linear_range;
				if (endSite - startSite +1 >= linear_range_size_min) {
					if (2*startSite >= endSite +2) {	
						linear_range = chrTAGA.substring(2*startSite-endSite -2, startSite-1);
					}else {
						linear_range = chrTAGA.substring(0, startSite-1);
					}
				}else {
                    if(startSite -1 >= linear_range_size_min) {
                    	linear_range = chrTAGA.substring(startSite-linear_range_size_min-1, startSite-1);
                    }else {
                    	linear_range = chrTAGA.substring(0, startSite-1);
					}										
				}
				aligner.setSeq(circ_range_seq.substring(circRangeLen-len_str,circRangeLen),str);
				String[] alignment = aligner.getAlignment();
				if (aligner.getAlignmentScore()>=len_str-2*(Initial+len_str/10) && alignment[1].length() >= len_str -2) {
					if (circLineArr[0].equals("1")) {
						pem_null_range_seq = linear_range;
					}
					lable = false;
				}	
				if(lable) {												
					//Comparison of detected seeds in the two regions
					//Seed is searched iteratively in the descending order of length
					int tag = IIC1_1.isInCircRNA1_1(len_str, str, circ_range_seq, linear_range);	
					if (tag == 1) {
						if (quant >= minMapqUni) {
							if (circLineArr[0].equals("1")) {
								pem_null_range_seq = linear_range;
							}
						}else {
							return "2";
						}												
					}else {
						return "0";
					}								
				}					
			}							
		}else {		
			if(str.length()<5) {
				if (!circ_range_seq.substring(0,circRangeLen).equals(str)) {
					return "0";
				}else {
					return "2";	
				}
			}else {
				int Initial = 1;
				initial_seq = str.substring(0,initial_size2);
				if (circ_range_seq.substring(0,initial_size2).equals(initial_seq)) {
					Initial = 0;
				}else {
					aligner.setSeq(circ_range_seq.substring(0,initial_size2),initial_seq);
					String[] alignment = aligner.getAlignment();
					if (aligner.getAlignmentScore()>=3 && alignment[0].length()>=4) {
						
					}else {
						return "0";
					}
				}
				boolean lable = true;
				String linear_range;
				if (endSite - startSite +1 >= linear_range_size_min) {	
					if(2*endSite-startSite+1 > chrTAGA.length()) {
                    	linear_range = chrTAGA.substring(endSite,chrTAGA.length());
                    }else {
                    	linear_range = chrTAGA.substring(endSite,2*endSite-startSite+1);
					}						
				}else {
                    if(endSite+linear_range_size_min>chrTAGA.length()) {
                    	linear_range = chrTAGA.substring(endSite,chrTAGA.length());
                    }else {
                    	linear_range = chrTAGA.substring(endSite,endSite+linear_range_size_min);
					}									
				}
				aligner.setSeq(circ_range_seq.substring(0,len_str),str);
				String[] alignment = aligner.getAlignment();
				if (aligner.getAlignmentScore() >= str.length()-2*(Initial+len_str/10)  && alignment[1].length() >= len_str -2) {
					if (circLineArr[0].equals("1")) {
						pem_null_range_seq = linear_range;
					}
					lable = false;
				}
				if(lable) {												
					//Comparison of detected seeds in the two regions
					//Seed is searched iteratively in the descending order of length
					int tag = IIC1_2.isInCircRNA1_2(len_str, str, circ_range_seq, linear_range);		
					if (tag == 1) {
						if (quant >= minMapqUni) {
							if (circLineArr[0].equals("0")) {
								pem_null_range_seq = linear_range;
							}
						}else {
							return "2";
						}												
					}else {
						return "0";
					}								
				}
			}																										
		}
		if (!circLineArr[7].equals("*")) {
			if (IIC2.isInCircRNA2(circLineArr[7], circ_range_seq) == 0) {
				return "0";
			}
		}			
		if (circLineArr[6].length()>5) {
			return IIC3.isInCircRNA3(circLineArr[6],circLineArr[8],circ_range_seq,pem_null_range_seq)+"2";		
		}
		return "12";				
	}
}
