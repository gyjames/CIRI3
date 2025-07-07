package com.zx.findcircrna;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import com.zx.hg38.IsBSJHg1;
import com.zx.hg38.IsBSJIntronHg1;

public class IsBSJScan1 {
	int minMapqUni,maxCircle,minCircle;
	Misd misd = new Misd();
	CompRev compRev = new CompRev();
	IsBSJHg1 isBSJHg1;
	HashMap<String, String> chrTCGAMap;
	String mitochondrion;
	boolean mlable;
	public IsBSJScan1(int minMapqUni, int maxCircle, int minCircle,int linear_range_size_min,boolean intronLable,
			HashMap<String, String> chrExonStartMap,HashMap<String, String> chrExonEndMap,HashMap<String, String> chrTCGAMap,
			HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,HashMap<String, ArrayList<String>> chrExonEndTranscriptMap,String mitochondrion,
			boolean mlable,boolean spLable) {
		super();
		this.minMapqUni = minMapqUni;
		this.maxCircle = maxCircle;
		this.minCircle = minCircle;
		this.chrTCGAMap = chrTCGAMap;	
		this.mitochondrion = mitochondrion;	
		this.mlable = mlable;	
		if (intronLable) {
			isBSJHg1 = new IsBSJIntronHg1(linear_range_size_min, chrExonStartMap, chrExonEndMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,spLable);
		}else {
			isBSJHg1 = new IsBSJHg1(linear_range_size_min, chrExonStartMap, chrExonEndMap,chrExonStartTranscriptMap,chrExonEndTranscriptMap,mitochondrion,spLable);
		}		
	}
	
	public String isBSJScan1(HashMap<Integer, ArrayList<String[]>> readsMap,HashMap<Integer, String> standMap) throws IOException {
		int isParied = readsMap.keySet().size()-1;		
		for (int n = 0; n <= isParied; n++) {
			int seqLen = standMap.get(n).length() - 1;
			for (int i = 0; i < readsMap.get(n).size()-1; i++) {
				for (int j = i+1; j < readsMap.get(n).size(); j++) {
					String[] aligment1 = readsMap.get(n).get(i);
					String[] aligment2 = readsMap.get(n).get(j);				
					if((!aligment1[1].equals(mitochondrion) || mlable) && aligment1[1].equals(aligment2[1]) && aligment1[0].equals(aligment2[0]) && 
							!aligment1[4].equals("*") && !aligment2[4].equals("*") && 
							(Integer.valueOf(aligment1[3]) >= minMapqUni ||  Integer.valueOf(aligment2[3]) >= minMapqUni)) {
						
						aligment1[4] = aligment1[4].replaceAll("H", "S");
						aligment2[4] = aligment2[4].replaceAll("H", "S");
						int[] CIGARite1 = misd.misd(aligment1[4], seqLen);
						int[] CIGARite2 = misd.misd(aligment2[4], seqLen);
						//对CIGAR排序,升序
						if (CIGARite1[0] > CIGARite2[0]) {
							String[] aligmentTem = aligment1;
							aligment1 = aligment2;
							aligment2 = aligmentTem;
							int[] CIGARiteTem = CIGARite1;
							CIGARite1 = CIGARite2;
							CIGARite2 = CIGARiteTem;
						}
						int alingStart1 = Integer.valueOf(aligment1[2]);
						int alingStart2 = Integer.valueOf(aligment2[2]);
						int alingMQ1 = Integer.valueOf(aligment1[3]);
						int alingMQ2 = Integer.valueOf(aligment2[3]);
						String str3 = "*",str2,str1,str4;
						//Two-segmment
						//CIGAR values reflecting potential PCC signal in the form of upstream xS/HyM and downstream xMyS/H, where x and y represent the number of mapping (M), soft clipping (S) or hard clipping (H) bases.
						if (CIGARite1[0]*CIGARite2[0] == -1) {
							int cirScale = CIGARite1[0]*alingStart1+CIGARite1[2]+ CIGARite2[0]*alingStart2+CIGARite2[2];
							if (cirScale > 0 && Math.abs(CIGARite1[1]-CIGARite2[1])<=6 && cirScale <= maxCircle && cirScale >= minCircle) {
								//adjust the putative boundaris of the potential BSJ junction
								int end_adjustment1 = (int)((CIGARite1[1]*CIGARite1[0]+CIGARite2[1]*CIGARite2[0])/2);
								int end_adjustment2 = 	CIGARite1[1]*CIGARite1[0]+CIGARite2[1]*CIGARite2[0]-end_adjustment1;
								if (Math.abs(end_adjustment1)>4) {
									continue;
								}
								//extract the corresponding sequence of each segment in the read
								if (aligment1[0].equalsIgnoreCase(standMap.get(n).substring(0,1)) ) {
									str2 = standMap.get(n).substring(1).substring(0,CIGARite1[1]+end_adjustment1);
									str1 = standMap.get(n).substring(1).substring(CIGARite1[1]+end_adjustment1);
								}else {
									str2 = compRev.compRev(standMap.get(n).substring(1)).substring(0,CIGARite1[1]+end_adjustment1);
									str1 = compRev.compRev(standMap.get(n).substring(1)).substring(CIGARite1[1]+end_adjustment1);
								}
								
								int str4_ok = 0;
                                if (isParied == 1) {
                                	if (!aligment1[0].equals(standMap.get(1-n).substring(0,1)) ) {
    									str4 = standMap.get(1-n).substring(1) ;
    								}else {
    									str4 = compRev.compRev(standMap.get(1-n).substring(1)) ;
    								}																	
    								for (int k = 0; k < readsMap.get(1-n).size(); k++) {
    									String[] aligmentAno = readsMap.get(1-n).get(k);
    									int[] CIGARiteAno = misd.misd(aligmentAno[4], seqLen);
    									if (aligmentAno[1].equals(aligment1[1]) && Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
    										if(!aligmentAno[0].equals(aligment1[0]) && Integer.valueOf(aligmentAno[2])>=alingStart1+end_adjustment1-6 &&
    												Integer.valueOf(aligmentAno[2])+CIGARiteAno[3]<=alingStart2+CIGARite2[3]-end_adjustment2+6) {
    											str4_ok = 1;
    											break;
    										}else {
    											str4_ok = -1;
    											break;
    										}
    									}
    									
    								}
								}else {
									str4_ok = 1;
									str4 = "";
								}								
								int aligmentQ1 = 0,aligmentQ2 =0,sumQ = 0;
								if (alingMQ1>=minMapqUni && alingMQ2>=minMapqUni) {
									aligmentQ1 = 1;
									aligmentQ2 = 1;	
									sumQ = 1;
								}else if (alingMQ1>=minMapqUni) {
									aligmentQ1 = 1;
								}else if (alingMQ2>=minMapqUni) {
									aligmentQ2 = 1;
								}
								String[] circCandidate = {aligment1[0],aligment1[1],str1,str2,str3,str4,aligmentQ1+"",aligmentQ2+"",str4_ok+"",
										(alingStart1+end_adjustment1)+"",(alingStart2+CIGARite2[3]-1-end_adjustment2)+"",end_adjustment1+"",end_adjustment2+""};										
								String chrTAGA = chrTCGAMap.get(aligment1[1]);
								String circInfor = isBSJHg1.isBSJHg1(circCandidate, chrTAGA,sumQ);
							
								if(circInfor != null) {
									return aligment1[4]+";"+aligment2[4]+"\t"+circInfor;
								}
							}						
							//Multiple-segmment
							//CIGAR values reflecting potential PCC signal in the form of xS/HyMzS/H and corresponding (x + y)S/HzM and/or xM(y + z)S/H, where x, y and z represent the number of mapping (M), soft clipping (S) or hard clipping (H) bases.											
						}else if (Math.abs(CIGARite1[0]*CIGARite2[0]) == 10) {
							//for potential PCC signal in the form of xS/HyMzS/H and corresponding (x + y)S/HzM
							if (CIGARite1[0] == -1) {
								int cirScale = alingStart2+CIGARite2[3]-1-alingStart1;
								if (cirScale > 0 && Math.abs(seqLen-CIGARite2[2]-CIGARite1[1])<=6 && cirScale <= maxCircle && cirScale >= minCircle) {
									//adjust the putative boundaris of the potential BSJ junction									
									int end_adjustment1 = (int)((CIGARite2[1]+CIGARite2[3]-CIGARite1[1])/2);
									int end_adjustment2 = 	CIGARite2[1]+CIGARite2[3]-CIGARite1[1]-end_adjustment1;	
									if (Math.abs(end_adjustment1)>4) {
										continue;
									}
									//extract the corresponding sequence of each segment in the read
									if (aligment1[0].equals(standMap.get(n).substring(0,1)) ) {
										str1 = standMap.get(n).substring(1).substring(CIGARite1[1]+end_adjustment1);
										str2 = standMap.get(n).substring(1).substring(CIGARite2[1],CIGARite1[1]+end_adjustment1);
										str3 = standMap.get(n).substring(1).substring(0,CIGARite2[1]);
									}else {
										str1 = compRev.compRev(standMap.get(n).substring(1)).substring(CIGARite1[1]+end_adjustment1);
										str2 = compRev.compRev(standMap.get(n).substring(1)).substring(CIGARite2[1],CIGARite1[1]+end_adjustment1);
										str3 = compRev.compRev(standMap.get(n).substring(1)).substring(0,CIGARite2[1]);												
									}
									
									int str4_ok = 0;	
									if (isParied == 1) {
										if (!aligment1[0].equals(standMap.get(1-n).substring(0,1)) ) {
											str4 = standMap.get(1-n).substring(1) ;
										}else {
											str4 = compRev.compRev(standMap.get(1-n).substring(1)) ;
										}
																		
										for (int k = 0; k < readsMap.get(1-n).size(); k++) {
											String[] aligmentAno = readsMap.get(1-n).get(k);
											int[] CIGARiteAno = misd.misd(aligmentAno[4],seqLen);
											if (aligmentAno[1].equals(aligment1[1]) && Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												if(!aligmentAno[0].equals(aligment1[0]) && Integer.valueOf(aligmentAno[2])>=alingStart1+end_adjustment1-6 &&
														Integer.valueOf(aligmentAno[2])+CIGARiteAno[3]<=alingStart2+CIGARite2[3]-end_adjustment2+6) {
													str4_ok = 1;
													break;
												}else {
													str4_ok = -1;
													break;
												}
											}
											
										}
									}else {
										str4_ok = 1;
										str4 = "";
									}									
									//pass the info of the read to gtag_pem_repeat for further determination
									int aligmentQ1 = 0,aligmentQ2 =0, sumQ = 0;
									if (alingMQ1>=minMapqUni && alingMQ2>=minMapqUni) {
										aligmentQ1 = 1;
										aligmentQ2 = 1;
										sumQ = 1;
									}else if (alingMQ1>=minMapqUni) {
										aligmentQ1 = 1;
									}else if (alingMQ2>=minMapqUni) {	
										aligmentQ2 = 1;										
									}
									String[] circCandidate = {aligment1[0],aligment1[1],str1,str2,str3,str4,aligmentQ1+"",aligmentQ2+"",str4_ok+"",
											(alingStart1+end_adjustment1)+"",(alingStart2+CIGARite2[3]-1-end_adjustment2)+"",end_adjustment1+"",end_adjustment2+""};	
									String chrTAGA = chrTCGAMap.get(aligment1[1]);
									String circInfor = isBSJHg1.isBSJHg1(circCandidate, chrTAGA,sumQ);
									if(circInfor != null) {
										return aligment1[4]+";"+aligment2[4]+"\t"+circInfor;
									}
								}
								//id:0 stand:1 chr:2 pos:3 mapq:4 cigar:5 
								//for potential PCC signal in the form of xS/HyMzS/H and corresponding xM(y + z)S/H
							}else {
								int cirScale = alingStart1+CIGARite1[3]-1-alingStart2;
								if (cirScale > 0 && Math.abs(CIGARite1[1]-CIGARite2[1])<=6 && cirScale <= maxCircle && cirScale >= minCircle ) {
									//adjust the putative boundaris of the potential BSJ junction
									int end_adjustment1 = (int)((CIGARite1[1]-CIGARite2[1])/2);
									int end_adjustment2 = 	CIGARite1[1]-CIGARite2[1]-end_adjustment1;	
									if (Math.abs(end_adjustment1)>4) {
										continue;
									}
									//extract the corresponding sequence of each segment in the read
									if (aligment1[0].equals(standMap.get(n).substring(0,1)) ) {
										str2 = standMap.get(n).substring(1).substring(0,CIGARite1[1]-end_adjustment2);
										str1 = standMap.get(n).substring(1).substring(CIGARite1[1]-end_adjustment2,seqLen-CIGARite2[2]);
										str3 = standMap.get(n).substring(1).substring(seqLen-CIGARite2[2]);
									}else {
										str2 = compRev.compRev(standMap.get(n).substring(1)).substring(0,CIGARite1[1]-end_adjustment2);
										str1 = compRev.compRev(standMap.get(n).substring(1)).substring(CIGARite1[1]-end_adjustment2,seqLen-CIGARite2[2]);
										str3 = compRev.compRev(standMap.get(n).substring(1)).substring(seqLen-CIGARite2[2]);																						
									}
									int str4_ok = 0;	
									if (isParied == 1) {
										if (!aligment1[0].equals(standMap.get(1-n).substring(0,1)) ) {
											str4 = standMap.get(1-n).substring(1) ;
										}else {
											str4 = compRev.compRev(standMap.get(1-n).substring(1)) ;
										}
										for (int k = 0; k < readsMap.get(1-n).size(); k++) {
											String[] aligmentAno = readsMap.get(1-n).get(k);
											int[] CIGARiteAno = misd.misd(aligmentAno[4], seqLen);
											if (aligmentAno[1].equals(aligment1[1]) && Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												if(!aligmentAno[0].equals(aligment1[0]) && Integer.valueOf(aligmentAno[2])>=alingStart2+end_adjustment1-6 &&
														Integer.valueOf(aligmentAno[2])+CIGARiteAno[3]<=alingStart1+CIGARite1[3]-end_adjustment2+6) {
													str4_ok = 1;
													break;
												}else {
													str4_ok = -1;
													break;
												}
											}
											
										}
									}else {
										str4_ok = 1;
										str4 = "";
									}								
									//id:0 stand:1 chr:2 pos:3 mapq:4 cigar:5 
									//pass the info of the read to gtag_pem_repeat for further determination
									int aligmentQ1 = 0,aligmentQ2 =0, sumQ = 0;
									if (alingMQ1>=minMapqUni && alingMQ2>=minMapqUni) {
										aligmentQ1 = 1;
										aligmentQ2 = 1;
										sumQ = 1;
									}else if (alingMQ2>=minMapqUni) {
										aligmentQ1 = 1;
									}else if (alingMQ1>=minMapqUni) {
										aligmentQ2 = 1;
									}
									String[] circCandidate = {aligment1[0],aligment1[1],str1,str2,str3,str4,aligmentQ1+"",aligmentQ2+"",str4_ok+"",
											(alingStart2+end_adjustment1)+"",(alingStart1+CIGARite1[3]-1-end_adjustment2)+"",end_adjustment1+"",end_adjustment2+""};	
									String chrTAGA = chrTCGAMap.get(aligment1[1]);
									String circInfor = isBSJHg1.isBSJHg1(circCandidate, chrTAGA,sumQ);
									if(circInfor != null) {
										return aligment1[4]+";"+aligment2[4]+"\t"+circInfor;
									}
								}
							}								
						}
					}
		          }
					
					}
				}
		return null;
	}
}
