package com.zx.hg38;

import java.util.ArrayList;

import com.zx.findcircrna.DistanceLoci;

import smith.SmithWaterman;

public class IsBSJStarHg2 {
	int linear_range_size_min,minMapqUni,initial_size1 = 5 + 2;
	SmithWaterman aligner = new SmithWaterman( 1, -1, -1);
	DistanceLoci distanceLoci = new DistanceLoci();
	IsInCircRNA1_1 IIC1_1 = new IsInCircRNA1_1();
	IsInCircRNA1_2 IIC1_2 = new IsInCircRNA1_2();
	IsInCircRNA2 IIC2 = new IsInCircRNA2();
	IsInCircRNA3 IIC3 = new IsInCircRNA3();
	ArrayList<Integer> locusListTem = new  ArrayList<Integer> ();
	ArrayList<Integer> locus2ListTem = new  ArrayList<Integer> ();
	public IsBSJStarHg2(int linear_range_size_min,int minMapqUni) {
		super();
		this.linear_range_size_min = linear_range_size_min;
		this.minMapqUni = minMapqUni;
	}
	//存放结果	
	public String isBSJHg2(String circLineArr[], String chrTAGA) {
		//Variables are assigned according to inputted variables
		// Existence of str3 can be determined by number of inputted variables					
		    String judgeTag = "3";
		    //2:CIRI2 3:CIRI3 
		    String str,pem_null_range_seq = "",circ_range_seq,initial_seq;
			int startSite = Integer.valueOf(circLineArr[3]);
			int endSite = Integer.valueOf(circLineArr[4]);
			int quant = Integer.valueOf(circLineArr[12]);
			if (startSite-3 < 0 && endSite+2 >chrTAGA.length()) {
				circ_range_seq = chrTAGA.substring(0,chrTAGA.length());	
			}else if (startSite-3<0) {
				circ_range_seq = chrTAGA.substring(0,endSite+2);
			}else if (endSite+2 >chrTAGA.length()) {
				circ_range_seq = chrTAGA.substring(startSite-3,chrTAGA.length());
			}else {
				circ_range_seq = chrTAGA.substring(startSite-3,endSite+2);
			}
			int circRangeLen = circ_range_seq.length();
			int len_str = circLineArr[5].length()+2;
			//Alignment style: xS/HyM or xS/HyMzS/H
			if (circLineArr[2].equals("sm")) {				
				str = circLineArr[5]+circLineArr[11];
				if(len_str<7) {
					if (!circ_range_seq.substring(circRangeLen-len_str,circRangeLen).equals(str)) {
						return "0";
					}else {
						return "2";	
					}
				}else {
					if (circ_range_seq.substring(circRangeLen-4,circRangeLen-2).equals(str.substring(len_str-4,len_str-2)) || 
							circ_range_seq.substring(circRangeLen-7,circRangeLen-4).equals(str.substring(len_str-7,len_str-4))	) {						
						boolean lable = true;
						String linear_range;
						if (endSite - startSite +5 >= linear_range_size_min) {
							if (2*startSite >= endSite +6) {	
								linear_range = chrTAGA.substring(2*startSite-endSite -6, startSite-1);
							}else {
								linear_range = chrTAGA.substring(0, startSite-1);
							}
						}else {
	                        if(startSite - 1 >= linear_range_size_min) {
	                        	linear_range = chrTAGA.substring(startSite-linear_range_size_min-1, startSite-1);
	                        }else {
	                        	linear_range = chrTAGA.substring(0, startSite-1);
							}										
						}						
						//先用Smithwaterman	1~10错一个 11~20错两个	
						if(circRangeLen < len_str) {
							return "2";
						}
						aligner.setSeq(circ_range_seq.substring(circRangeLen-len_str,circRangeLen),str);
						aligner.getAlignment();
						//if (aligner.getAlignmentScore() >= len_str -2-(len_str-3)/10*2) {
						if (aligner.getAlignmentScore() >= len_str-(len_str-2)/10*2) {
							if (quant >= minMapqUni) {
								initial_seq = str.substring(len_str-initial_size1,len_str);
								if (circ_range_seq.substring(circRangeLen-initial_size1,circRangeLen).equals(initial_seq)) {								
									int tag = IIC1_1.isInCircRNA1_1(len_str, str, circ_range_seq, linear_range);	
									if (tag == 1) {
										judgeTag = "2";
									}
								}
							}						
							if (circLineArr[0].equals("1")) {
								pem_null_range_seq = linear_range;
							}
							lable = false;																			
						}
						if(lable) {
							initial_seq = str.substring(len_str-initial_size1,len_str);
							if (circ_range_seq.substring(circRangeLen-initial_size1,circRangeLen).equals(initial_seq)) {								
								//Comparison of detected seeds in the two regions
								//Seed is searched iteratively in the descending order of length
								int tag = IIC1_1.isInCircRNA1_1(len_str, str, circ_range_seq, linear_range);	
								if (tag == 1) {
									if (quant >= minMapqUni) {
										judgeTag = "2";
										if (circLineArr[0].equals("1")) {
											pem_null_range_seq = linear_range;
										}
									}else {
										return "2";
									}												
								}else {
									return "0";
								}
							}else {
								return "0";
							}
						}
					}else {
						return "0";
					}
					
					
				}
				
			}else {
				str = circLineArr[10]+circLineArr[5];
				if(len_str<7) {
					if (!circ_range_seq.substring(0,len_str).equals(str)) {
						return "0";
					}else {
						return "2";	
					}
				}else {
					if (circ_range_seq.substring(2,4).equals(str.substring(2,4)) || 
							circ_range_seq.substring(4,7).equals(str.substring(4,7))) {
						boolean lable = true;
						String linear_range;
						if (endSite - startSite +5 >= linear_range_size_min) {	
							if(2*endSite-startSite+5>chrTAGA.length()) {
	                        	linear_range = chrTAGA.substring(endSite,chrTAGA.length());
	                        }else {
	                        	linear_range = chrTAGA.substring(endSite,2*endSite-startSite+5);
							}
							
						}else {
	                        if(endSite+linear_range_size_min>chrTAGA.length()) {
	                        	linear_range = chrTAGA.substring(endSite,chrTAGA.length());
	                        }else {
	                        	linear_range = chrTAGA.substring(endSite,endSite+linear_range_size_min);
							}									
						}
						//先用Smithwaterman	1~10错一个 11~20错两个		
						if(circ_range_seq.length() < len_str) {
							return "2";
						}
						aligner.setSeq(circ_range_seq.substring(0,len_str),str);
						aligner.getAlignment();
						//if (aligner.getAlignmentScore() >= len_str -2-(len_str-3)/10*2) {
						if (aligner.getAlignmentScore() >= len_str  - (len_str-2)/10*2) {
							initial_seq = str.substring(0,initial_size1);
							if (circ_range_seq.substring(0,initial_size1).equals(initial_seq)) {								
								 int tag = IIC1_2.isInCircRNA1_2(len_str, str, circ_range_seq, linear_range);	
						    	 if (tag == 1) {
						    		 if (quant >= minMapqUni) {
						    			 judgeTag = "2";
										}	 										
									}									
							}							
							if (circLineArr[0].equals("0")) {
								pem_null_range_seq = linear_range;
							}
								lable = false;													
						}
						if(lable) {
							initial_seq = str.substring(0,initial_size1);
							if (circ_range_seq.substring(0,initial_size1).equals(initial_seq)) {								
								 int tag = IIC1_2.isInCircRNA1_2(len_str, str, circ_range_seq, linear_range);	
						    	 if (tag == 1) {
						    		 if (quant >= minMapqUni) {
						    			 judgeTag = "2";
						    			 if (circLineArr[0].equals("0")) {
												pem_null_range_seq = linear_range;
											}
										}else {
											return "2";
										}	 										
									}else {
										return "0";
									}									
							}else {
								return "0";
							}	
						}
					}else {
						return "0";
					}
					
											
				}
					
			}
			if (!circLineArr[7].equals("*")) {
				if (IIC2.isInCircRNA2(circLineArr[7], circ_range_seq) == 0) {
					return "0";
				}
			}
			//只是多了这个地方
			if (circLineArr[8].equals("2")) {
				return 1+judgeTag;		
			}else if (circLineArr[6].length()>5) {
				return IIC3.isInCircRNA3(circLineArr[6],circLineArr[8],circ_range_seq,pem_null_range_seq)+judgeTag;					
			}
			return 1+judgeTag;		
	}
}
