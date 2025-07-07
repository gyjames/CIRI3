package com.zx.findcircrna;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.hg38.IsBSJStarHg2;

public class IsBSJScan2Star {
	MisdStar misd = new MisdStar(); 
	CompRev compRev = new CompRev();
	int minMapqUni, seqLen;
	IsBSJStarHg2 isBSJHg2;
	GetFSJClass GFC;
	HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,chrSiteMap2;
	HashMap<String, String> chrTCGAMap;
	HashMap<String, Integer> circFSJNewMap = new HashMap<String, Integer>();
	HashSet<String> temFSJId = new HashSet<String>();
	HashMap<String, byte[]> siteArrayMap1,siteArrayMap2;
	public IsBSJScan2Star(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min, HashMap<String, byte[]> siteArrayMap1,
			HashMap<String, byte[]> siteArrayMap2, HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,
			HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2, HashMap<String, String> chrTCGAMap, int seqLen) throws IOException {
		super();
		this.minMapqUni = minMapqUni;
		this.siteArrayMap1 = siteArrayMap1;
		this.siteArrayMap2 = siteArrayMap2;
		this.chrSiteMap1 = chrSiteMap1;
		this.chrSiteMap2 = chrSiteMap2;
		this.chrTCGAMap = chrTCGAMap;
		this.seqLen = seqLen;
		isBSJHg2 = new IsBSJStarHg2(linear_range_size_min,minMapqUni);
		GFC = new GetFSJClass(seqLen);
		circFSJNewMap.putAll(circFSJMap);
	}
	public HashMap<String, Integer> getCircFSJMap() throws IOException {
		return circFSJNewMap;		
	}
	public void setFSJScan2List() {
		for (String circ : circFSJNewMap.keySet()) {
			circFSJNewMap.put(circ, 0);
		}
	}
	public String isCandidate(HashMap<Integer, ArrayList<String[]>> readsMap, HashMap<Integer, String> standMap) throws IOException {
		// 初始化记录FSJ的Set
		int isParied = readsMap.keySet().size()-1;		
		temFSJId.clear();
		for (int n = 0; n <= isParied; n++) {
			if (standMap.get(n) == null) {
				return null;
			}
			int temSeqLen = standMap.get(n).length() - 1;
			for (int i = 0; i < readsMap.get(n).size(); i++) {
				String[] aligment= readsMap.get(n).get(i);
				if (aligment[4].equals("*") || !siteArrayMap1.containsKey(aligment[1])) {
					continue;
				}
				HashMap<Integer, ArrayList<SiteSort>> siteMap1 = chrSiteMap1.get(aligment[1]);
				HashMap<Integer, ArrayList<SiteSort>> siteMap2 = chrSiteMap2.get(aligment[1]);				
				byte[] siteArray1 = siteArrayMap1.get(aligment[1]);
				byte[] siteArray2 = siteArrayMap2.get(aligment[1]);
				
				int startSite = Integer.valueOf(aligment[2]);
				ArrayList<String> CIGARiteList = misd.misdScan2(aligment[4],startSite);
				String[] alignmentArr = CIGARiteList.get(0).split("\t");
				String[] alignmentSiteArr = CIGARiteList.get(1).split("\t");
				int S1 = Integer.valueOf(alignmentArr[0]);
				int S2 = Integer.valueOf(alignmentArr[2]);
				int M = Integer.valueOf(alignmentArr[1]);
				ArrayList<int[]> faSiteList = new ArrayList<int[]>();
				for (String siteStr : alignmentSiteArr) {
					String[] siteStrArr = siteStr.split("-");
					faSiteList.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
				}
				
				int faStart = faSiteList.get(0)[0];
				int faEnd = faSiteList.get(faSiteList.size()-1)[1];
				
				if (S1 == 0 && S2 == 0) {
					for (int j = 0; j < faSiteList.size(); j++) {
						int[] siteArr = faSiteList.get(j);
						int startTem = siteArr[0];
						int endTem = siteArr[1];
						if (j == 0) {
							startTem = startTem + 6;
						}
						if (j == faSiteList.size() - 1) {
							endTem = endTem - 6;
						}
						//取出这个范围内的circRNA
						temFSJId.addAll(GFC.getFSJ(startTem, endTem,aligment[1],0,siteArray1,siteArray2,siteMap1,siteMap2));
					}
				} else {	
					//SM
					if (S1 != 0 && S2 == 0) {
						int newNumSite = faStart;// read在参考基因组匹配的位置
						// adjust the putative boundaris according to the candidate BSJ junction(s)
						int num1 = (newNumSite - 6) / seqLen;
						int num2 = (newNumSite + 6) / seqLen;
						if (num1 > 0 && siteArray1[num1] == 1) {
							ArrayList<SiteSort> siteList1 = siteMap1.get(num1);
							for (int k = siteList1.size() - 1; k >= 0; k--) {
								if (siteList1.get(k).getSite() - newNumSite >= (-6)) {
									if (siteList1.get(k).getSite() - newNumSite <= 6) {
										String[] temChrSite = siteList1.get(k).getLength();
										int bias = siteList1.get(k).getSite() - newNumSite;
										
										// extract the corresponding sequence of key segment in the read
										String str, str2;
										if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
											if (S1 + bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(0, S1 + bias);
											}
										} else {
											if (S1 + bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1)).substring(0,S1 + bias);
											}
										}
											// extract the corresponding sequence of its paired read
										int str2_ok = 0;
										if (isParied == 1) {
											if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
												str2 = standMap.get(1 - n).substring(1);
											} else {
												str2 = compRev.compRev(standMap.get(1 - n).substring(1));
											}												
											for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
												String[] aligmentAno = readsMap.get(1 - n).get(m);	
												if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
													
													int startSiteAno = Integer.valueOf(aligmentAno[2]);
													ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
													String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
													String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
													int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
													int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
													//int MAno = Integer.valueOf(alignmentAnoArr[1]);
													ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
													for (String siteStr : alignmentSiteAnoArr) {
														String[] siteStrArr = siteStr.split("-");
														faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
													}
													int faStartAno = faSiteListAno.get(0)[0];
													int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
													
													if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
															&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
														if(S1Ano == 0 && S2Ano == 0) {
															str2_ok = 2;
															
														}else {
															str2_ok = 1;
														}
														break;
													}else {
														str2_ok = -1;
														break;
													}
												}

											}
										}else {
											str2_ok = 1;
											str2 = "";
										}
										//如果str2_ok = 2，则不用看paried read
										String[] circCandidate = {aligment[0],aligment[1],"sm",temChrSite[0],temChrSite[1],str,str2,"*",
												str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};											
										String chrTAGA = chrTCGAMap.get(aligment[1]);										
										String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
										if (tag.equals("0")) {
											temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
										} else if (!tag.equals("2")) {
											String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
													+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);											
											return circInfor;											
										}									
									}

								} else {
									break;
								}
							}
						}
						if (num2 != num1 && siteArray1[num2] == 1) {
							ArrayList<SiteSort> siteList1 = siteMap1.get(num2);
							for (int k = 0; k < siteList1.size(); k++) {
								if (siteList1.get(k).getSite() - newNumSite <= 6) {
									String[] temChrSite = siteList1.get(k).getLength();								
									int bias = siteList1.get(k).getSite() - newNumSite;
									// extract the corresponding sequence of key segment in the read
									String str, str2;
									if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
										if (S1 + bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(0, S1 + bias);
										}
									} else {
										if (S1 + bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													S1 + bias);
										}
									}
									// extract the corresponding sequence of its paired read
									int str2_ok = 0;
									if (isParied == 1) {
										if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
											str2 = standMap.get(1 - n).substring(1);
										} else {
											str2 = compRev.compRev(standMap.get(1 - n).substring(1));
										}												
										for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
											String[] aligmentAno = readsMap.get(1 - n).get(m);	
											if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												
												int startSiteAno = Integer.valueOf(aligmentAno[2]);
												ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
												String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
												String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
												int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
												int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
												//int MAno = Integer.valueOf(alignmentAnoArr[1]);
												ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
												for (String siteStr : alignmentSiteAnoArr) {
													String[] siteStrArr = siteStr.split("-");
													faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
												}
												int faStartAno = faSiteListAno.get(0)[0];
												int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
												
												if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
														&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
													if(S1Ano == 0 && S2Ano == 0) {
														str2_ok = 2;
														
													}else {
														str2_ok = 1;
													}
													break;
												}else {
													str2_ok = -1;
													break;
												}
											}

										}
									}else {
										str2_ok = 1;
										str2 = "";
									}
									String[] circCandidate = {aligment[0],aligment[1],"sm",temChrSite[0],temChrSite[1],str,str2,"*",
											str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
									String chrTAGA = chrTCGAMap.get(aligment[1]);
									String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
									if (tag.equals("0")) {
										temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
									} else if (!tag.equals("2")) {
										String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
												+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);											
										return circInfor;										
									}

								} else {
									break;
								}
							}
						}

						// MS
					} else if (S1 == 0 && S2 != 0) {
						int newNumSite = faEnd;
						// adjust the putative boundaris according to the candidate BSJ junction(s)
						int num1 = (newNumSite - 6) / seqLen;
						int num2 = (newNumSite + 6) / seqLen;
						if (siteArray2[num1] == 1) {
							ArrayList<SiteSort> siteList2 = siteMap2.get(num1);
							for (int k = siteList2.size() - 1; k >= 0; k--) {
								if (siteList2.get(k).getSite() - newNumSite >= (-6)) {
									if (siteList2.get(k).getSite() - newNumSite <= 6) {
										String[] temChrSite = siteList2.get(k).getLength();
										int bias = siteList2.get(k).getSite() - newNumSite;
										// extract the corresponding sequence of key segment in the read
										String str, str2;
										if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
											if (M + bias >= temSeqLen) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(M + bias);
											}

										} else {
											if (M + bias >= temSeqLen) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1))
														.substring(M + bias);
											}

										}							
											// extract the corresponding sequence of its paired read
										int str2_ok = 0;
										if (isParied == 1) {
											if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
												str2 = standMap.get(1 - n).substring(1);
											} else {
												str2 = compRev.compRev(standMap.get(1 - n).substring(1));
											}
											for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
												String[] aligmentAno = readsMap.get(1 - n).get(m);
												if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
													
													int startSiteAno = Integer.valueOf(aligmentAno[2]);
													ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
													String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
													String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
													int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
													int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
													//int MAno = Integer.valueOf(alignmentAnoArr[1]);
													ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
													for (String siteStr : alignmentSiteAnoArr) {
														String[] siteStrArr = siteStr.split("-");
														faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
													}
													int faStartAno = faSiteListAno.get(0)[0];
													int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
													
													if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
															&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
														if(S1Ano == 0 && S2Ano == 0) {
															str2_ok = 2;
															
														}else {
															str2_ok = 1;
														}
														break;
													}else {
														str2_ok = -1;
														break;
													}
												}

											}
										}else {
											str2_ok = 1;
											str2 = "";
										}
										String[] circCandidate = {aligment[0],aligment[1],"ms",temChrSite[0],temChrSite[1],str,str2,"*",
												str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
											
										String chrTAGA = chrTCGAMap.get(aligment[1]);
										String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
										if (tag.equals("0")) {
											temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
										} else if (!tag.equals("2")) {
											String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
													+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);												
											return circInfor;											
										}

									}

								} else {
									break;
								}
							}
						}
						if (num2 != num1 && num2 < siteArray2.length && siteArray2[num2] == 1) {
							ArrayList<SiteSort> siteList2 = siteMap2.get(num2);
							for (int k = 0; k < siteList2.size(); k++) {
								if (siteList2.get(k).getSite() - newNumSite <= 6) {
									String[] temChrSite = siteList2.get(k).getLength();
									int bias = siteList2.get(k).getSite() - newNumSite;
									// extract the corresponding sequence of key segment in the read
									String str, str2;
									if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
										if (M + bias >= temSeqLen) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(M + bias);
										}

									} else {
										if (M + bias >= temSeqLen) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1))
													.substring(M + bias);
										}

									}
									
										// extract the corresponding sequence of its paired read
									int str2_ok = 0;
									if (isParied == 1) {
										if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
											str2 = standMap.get(1 - n).substring(1);
										} else {
											str2 = compRev.compRev(standMap.get(1 - n).substring(1));
										}
										for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
											String[] aligmentAno = readsMap.get(1 - n).get(m);	
											if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												
												int startSiteAno = Integer.valueOf(aligmentAno[2]);
												ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
												String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
												String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
												int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
												int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
												//int MAno = Integer.valueOf(alignmentAnoArr[1]);
												ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
												for (String siteStr : alignmentSiteAnoArr) {
													String[] siteStrArr = siteStr.split("-");
													faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
												}
												int faStartAno = faSiteListAno.get(0)[0];
												int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
												
												if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
														&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
													if(S1Ano == 0 && S2Ano == 0) {
														str2_ok = 2;
														
													}else {
														str2_ok = 1;
													}
													break;
												}else {
													str2_ok = -1;
													break;
												}
											}

										}
									}else {
										str2_ok = 1;
										str2 = "";
									}											
									String[] circCandidate = {aligment[0],aligment[1],"ms",temChrSite[0],temChrSite[1],str,str2,"*",
											str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
									String chrTAGA = chrTCGAMap.get(aligment[1]);
									String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
									if (tag.equals("0")) {
										temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
									} else if (!tag.equals("2")) {
										String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
												+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);											
										return circInfor;										
									}
									

								} else {
									break;
								}
							}
						}

					} else if (S1 != 0 && S2 != 0) {
						String str, str2, str3;
						int newNumSite = faStart;
						int num1 = (newNumSite - 6) / seqLen;
						int num2 = (newNumSite + 6) / seqLen;
						if (num1 > 0 && siteArray1[num1] == 1) {
							ArrayList<SiteSort> siteList1 = siteMap1.get(num1);
							for (int k = siteList1.size() - 1; k >= 0; k--) {
								// adjust the putative boundaris according to the candidate BSJ junction(s)
								if (siteList1.get(k).getSite() - newNumSite >= (-6)) {
									if (siteList1.get(k).getSite() - newNumSite <= 6) {
										String[] temChrSite = siteList1.get(k).getLength();
										int bias = siteList1.get(k).getSite() - newNumSite;
										// extract the corresponding sequence of key segment in the read
										if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
											if (S1 + bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(0, S1 + bias);
											}

											str3 = standMap.get(n).substring(1)
													.substring(temSeqLen - S2);
										} else {
											if (S1 + bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
														S1 + bias);
											}
                                           
											str3 = compRev.compRev(standMap.get(n).substring(1))
													.substring(temSeqLen - S2);
										}
									
											// extract the corresponding sequence of its paired read
										int str2_ok = 0;
										if (isParied == 1) {
											if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
												str2 = standMap.get(1 - n).substring(1);
											} else {
												str2 = compRev.compRev(standMap.get(1 - n).substring(1));
											}
											for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
												String[] aligmentAno = readsMap.get(1 - n).get(m);	
												if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
													
													int startSiteAno = Integer.valueOf(aligmentAno[2]);
													ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
													String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
													String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
													int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
													int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
													//int MAno = Integer.valueOf(alignmentAnoArr[1]);
													ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
													for (String siteStr : alignmentSiteAnoArr) {
														String[] siteStrArr = siteStr.split("-");
														faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
													}
													int faStartAno = faSiteListAno.get(0)[0];
													int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
													
													if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
															&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
														if(S1Ano == 0 && S2Ano == 0) {
															str2_ok = 2;
															
														}else {
															str2_ok = 1;
														}
														break;
													}else {
														str2_ok = -1;
														break;
													}
												}

											}
										}else {
											str2_ok = 1;
											str2 = "";
										}	
										String[] circCandidate = {aligment[0],aligment[1],"sm",temChrSite[0],temChrSite[1],str,str2,str3,
												str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
										String chrTAGA = chrTCGAMap.get(aligment[1]);
										String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
										if (tag.equals("0")) {
											temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
										} else if (!tag.equals("2")) {
											String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
													+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);												
											return circInfor;											
										}

											
										
									}
								}
							}
						}
						if (num2 != num1 && siteArray1[num2] == 1) {
							ArrayList<SiteSort> siteList1 = siteMap1.get(num2);
							for (int k = 0; k < siteList1.size(); k++) {

								if (siteList1.get(k).getSite() - newNumSite <= 6) {
									String[] temChrSite = siteList1.get(k).getLength();
									
									int bias = siteList1.get(k).getSite() - newNumSite;
									// extract the corresponding sequence of key segment in the read
									if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
										if (S1 + bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(0, S1 + bias);
										}

										str3 = standMap.get(n).substring(1).substring(temSeqLen - S2);
									} else {
										if (S1 + bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													S1+ bias);
										}

										str3 = compRev.compRev(standMap.get(n).substring(1))
												.substring(temSeqLen - S2);
									}
									
										// extract the corresponding sequence of its paired read
									int str2_ok = 0;
									if (isParied == 1) {
										if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
											str2 = standMap.get(1 - n).substring(1);
										} else {
											str2 = compRev.compRev(standMap.get(1 - n).substring(1));
										}
										for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
											String[] aligmentAno = readsMap.get(1 - n).get(m);	
											if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												
												int startSiteAno = Integer.valueOf(aligmentAno[2]);
												ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
												String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
												String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
												int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
												int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
												//int MAno = Integer.valueOf(alignmentAnoArr[1]);
												ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
												for (String siteStr : alignmentSiteAnoArr) {
													String[] siteStrArr = siteStr.split("-");
													faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
												}
												int faStartAno = faSiteListAno.get(0)[0];
												int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
												
												if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
														&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
													if(S1Ano == 0 && S2Ano == 0) {
														str2_ok = 2;
														
													}else {
														str2_ok = 1;
													}
													break;
												}else {
													str2_ok = -1;
													break;
												}
											}

										}
									}else {
										str2_ok = 1;
										str2 = "";
									}											
									String[] circCandidate = {aligment[0],aligment[1],"sm",temChrSite[0],temChrSite[1],str,str2,str3,
											str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
									String chrTAGA = chrTCGAMap.get(aligment[1]);
									String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
									if (tag.equals("0")) {
										temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
									} else if (!tag.equals("2")) {
										String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
												+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);	
										
										return circInfor;
										
									}


								} else {
									break;
								}
							}
						}

						newNumSite = faEnd;
						// adjust the putative boundaris according to the candidate BSJ junction(s)
						num1 = (newNumSite - 6) / seqLen;
						num2 = (newNumSite + 6) / seqLen;
						if (siteArray2[num1] == 1) {
							ArrayList<SiteSort> siteList2 = siteMap2.get(num1);
							for (int k = siteList2.size() - 1; k >= 0; k--) {
								if (siteList2.get(k).getSite() - newNumSite >= (-6)) {
									if (siteList2.get(k).getSite() - newNumSite <= 6) {
										String[] temChrSite = siteList2.get(k).getLength();
										int bias = siteList2.get(k).getSite() - newNumSite;
										// extract the corresponding sequence of key segment in the read
										// extract the corresponding sequence of key segment in the read
										if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
											if (S2 - bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1)
														.substring(temSeqLen - S2 + bias);
											}
											str3 = standMap.get(n).substring(1).substring(0, S1);
										} else {
											if (S2 - bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1))
														.substring(temSeqLen -S2 + bias);
											}

											str3 = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													S1);
										}
										
											// extract the corresponding sequence of its paired read
										int str2_ok = 0;
										if (isParied == 1) {
											if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
												str2 = standMap.get(1 - n).substring(1);
											} else {
												str2 = compRev.compRev(standMap.get(1 - n).substring(1));
											}
											for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
												String[] aligmentAno = readsMap.get(1 - n).get(m);	
												if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
													
													int startSiteAno = Integer.valueOf(aligmentAno[2]);
													ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
													String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
													String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
													int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
													int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
													//int MAno = Integer.valueOf(alignmentAnoArr[1]);
													ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
													for (String siteStr : alignmentSiteAnoArr) {
														String[] siteStrArr = siteStr.split("-");
														faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
													}
													int faStartAno = faSiteListAno.get(0)[0];
													int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
													
													if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
															&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
														if(S1Ano == 0 && S2Ano == 0) {
															str2_ok = 2;
															
														}else {
															str2_ok = 1;
														}
														break;
													}else {
														str2_ok = -1;
														break;
													}
												}

											}
										}else {
											str2_ok = 1;
											str2 = "";
										}
										String[] circCandidate = {aligment[0],aligment[1],"ms",temChrSite[0],temChrSite[1],str,str2,str3,
												str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
										String chrTAGA = chrTCGAMap.get(aligment[1]);
										String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
										if (tag.equals("0")) {
											temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
										} else if (!tag.equals("2")) {
											String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
													+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);	
											
											return circInfor;
											
										}											
										
									}

								} else {
									break;
								}
							}
						}
						if (num2 != num1 && num2 < siteArray2.length && siteArray2[num2] == 1) {
							ArrayList<SiteSort> siteList2 = siteMap2.get(num2);
							for (int k = 0; k < siteList2.size(); k++) {

								if (siteList2.get(k).getSite() - newNumSite <= 6) {
									String[] temChrSite = siteList2.get(k).getLength();
									
									int bias = siteList2.get(k).getSite() - newNumSite;
									// extract the corresponding sequence of key segment in the read
									// extract the corresponding sequence of key segment in the read
									if (aligment[0].equals(standMap.get(n).substring(0, 1))) {
										if (S2 - bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1)
													.substring(temSeqLen - S2 + bias);
										}

										str3 = standMap.get(n).substring(1).substring(0,S1);
									} else {
										if (S2 - bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1))
													.substring(temSeqLen - S2 + bias);
										}
										str3 = compRev.compRev(standMap.get(n).substring(1)).substring(0,S1);
									}
									
										// extract the corresponding sequence of its paired read
									int str2_ok = 0;
									if (isParied == 1) {
										if (!standMap.get(1 - n).substring(0, 1).equals(aligment[0])) {
											str2 = standMap.get(1 - n).substring(1);
										} else {
											str2 = compRev.compRev(standMap.get(1 - n).substring(1));
										}
										for (int m = 0; m < readsMap.get(1 - n).size(); m++) {
											String[] aligmentAno = readsMap.get(1 - n).get(m);	
											if (aligmentAno[1].equals(aligment[1])&& Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												
												int startSiteAno = Integer.valueOf(aligmentAno[2]);
												ArrayList<String> CIGARiteListAno = misd.misdScan2(aligmentAno[4],startSiteAno);
												String[] alignmentAnoArr = CIGARiteListAno.get(0).split("\t");
												String[] alignmentSiteAnoArr = CIGARiteListAno.get(1).split("\t");
												int S1Ano = Integer.valueOf(alignmentAnoArr[0]);
												int S2Ano = Integer.valueOf(alignmentAnoArr[2]);
												//int MAno = Integer.valueOf(alignmentAnoArr[1]);
												ArrayList<int[]> faSiteListAno = new ArrayList<int[]>();
												for (String siteStr : alignmentSiteAnoArr) {
													String[] siteStrArr = siteStr.split("-");
													faSiteListAno.add(new int[]{Integer.valueOf(siteStrArr[0]),Integer.valueOf(siteStrArr[1])});
												}
												int faStartAno = faSiteListAno.get(0)[0];
												int faEndAno = faSiteListAno.get(faSiteListAno.size()-1)[1];
												
												if (!aligmentAno[0].equals(aligment[0])&& faStartAno >= Integer.valueOf(temChrSite[0]) - 6
														&& faEndAno <= Integer.valueOf(temChrSite[1]) + 6) {
													if(S1Ano == 0 && S2Ano == 0) {
														str2_ok = 2;
														
													}else {
														str2_ok = 1;
													}
													break;
												}else {
													str2_ok = -1;
													break;
												}
											}

										}
									}else {
										str2_ok = 1;
										str2 = "";
									}
									String[] circCandidate = {aligment[0],aligment[1],"ms",temChrSite[0],temChrSite[1],str,str2,str3,
											str2_ok+"",temChrSite[2],temChrSite[3],temChrSite[4],aligment[3]};
										
									String chrTAGA = chrTCGAMap.get(aligment[1]);
									String tag = isBSJHg2.isBSJHg2(circCandidate, chrTAGA);
									if (tag.equals("0")) {
										temFSJId.add(aligment[1] + "\t" + temChrSite[0] + "\t"+ temChrSite[1]);
									} else if (!tag.equals("2")) {
										String circInfor = aligment[4]+ "\t"+ tag.substring(0,tag.length()-1) + "\t" + aligment[1] + "\t" + temChrSite[0]+"\t"+temChrSite[1] + "\t" + temChrSite[2]
												+ "\t" + temChrSite[3]+ "\t" + temChrSite[4]+ "\t" +tag.substring(tag.length()-1);	
										
										return circInfor;
										
									}
								} else {
									break;
								}
							}
						}
					} //

					for (int j = 0; j < faSiteList.size(); j++) {
						int[] siteArr = faSiteList.get(j);
						int startTem = siteArr[0];
						int endTem = siteArr[1];
						if (j == 0) {
							startTem = startTem + 6;
						}
						if (j == faSiteList.size() - 1) {
							endTem = endTem - 6;
						}
						//取出这个范围内的circRNA
						temFSJId.addAll(GFC.getFSJ(startTem, endTem,aligment[1],0,siteArray1,siteArray2,siteMap1,siteMap2));
					}
				}
			}
		}
		
		for (String circKey : temFSJId) {			
			int temNum = circFSJNewMap.get(circKey);
			circFSJNewMap.put(circKey, temNum + 1);				
		}
		return null;	
	}
}
