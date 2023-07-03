package com.zx.findcircrna;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.hg19.IsBSJHg2;

public class IsBSJScan2 {
	Misd misd = new Misd();
	CompRev compRev = new CompRev();
	int minMapqUni, seqLen;
	IsBSJHg2 isBSJHg2;
	GetFSJClass GFC;
	HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1,chrSiteMap2;
	HashMap<String, String> chrTCGAMap;
	HashMap<String, Integer> circFSJNewMap = new HashMap<String, Integer>();
	HashSet<String> temFSJId = new HashSet<String>();
	HashMap<String, byte[]> siteArrayMap1,siteArrayMap2;
	public IsBSJScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min, HashMap<String, byte[]> siteArrayMap1,
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
		isBSJHg2 = new IsBSJHg2(linear_range_size_min,minMapqUni);
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
				if (aligment[4].equals(temSeqLen + "M")) {
					int newNumSite = Integer.valueOf(aligment[2]);// read在参考基因组匹配的位置
					int startTem = newNumSite + 6;
					int endTem = newNumSite + temSeqLen - 7;
					//取出这个范围内的circRNA
					temFSJId.addAll(GFC.getFSJ(startTem, endTem,aligment[1],0,siteArray1,siteArray2,siteMap1,siteMap2));
				} else {					
					aligment[4] = aligment[4].replaceAll("H", "S");
					int[] CIGARite = misd.misd(aligment[4], temSeqLen);
					if (CIGARite[0] == -1) {
						int newNumSite = Integer.valueOf(aligment[2]);// read在参考基因组匹配的位置
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
											if (CIGARite[1] + bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(0, CIGARite[1] + bias);
											}
										} else {
											if (CIGARite[1] + bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
														CIGARite[1] + bias);
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
													int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
													if (!aligmentAno[0].equals(aligment[0])&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
															&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
														str2_ok = 1;
														break;
													} else {
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
										if (CIGARite[1] + bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(0, CIGARite[1] + bias);
										}
									} else {
										if (CIGARite[1] + bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													CIGARite[1] + bias);
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
												int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
												if (!aligmentAno[0].equals(aligment[0])&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
														&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3] - 1 <= Integer.valueOf(temChrSite[1]) + 6) {
													str2_ok = 1;
													break;
												} else {
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

						// Alignment style: xMyS/H
					} else if (CIGARite[0] == 1) {
						int newNumSite = Integer.valueOf(aligment[2]) + CIGARite[3] - 1;
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
											if (CIGARite[1] + bias >= temSeqLen) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(CIGARite[1] + bias);
											}

										} else {
											if (CIGARite[1] + bias >= temSeqLen) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1))
														.substring(CIGARite[1] + bias);
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
													int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
													if (!aligmentAno[0].equals(aligment[0])&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
															&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
														str2_ok = 1;
														break;
													} else {
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
										if (CIGARite[1] + bias >= temSeqLen) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(CIGARite[1] + bias);
										}

									} else {
										if (CIGARite[1] + bias >= temSeqLen) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1))
													.substring(CIGARite[1] + bias);
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
											if (aligmentAno[1].equals(aligment[1]) && Integer.valueOf(aligmentAno[3]) >= minMapqUni) {
												int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
												if (!aligmentAno[0].equals(aligment[0])&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
																&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
													str2_ok = 1;
													break;
												} else {
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

					} else if (CIGARite[0] == 10) {
						String str, str2, str3;
						int newNumSite = Integer.valueOf(aligment[2]);
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
											if (CIGARite[1] + bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1).substring(0, CIGARite[1] + bias);
											}

											str3 = standMap.get(n).substring(1)
													.substring(temSeqLen - CIGARite[2]);
										} else {
											if (CIGARite[1] + bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
														CIGARite[1] + bias);
											}
                                           
											str3 = compRev.compRev(standMap.get(n).substring(1))
													.substring(temSeqLen - CIGARite[2]);
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
													int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
													if (!aligmentAno[0].equals(aligment[0])
															&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
															&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
														str2_ok = 1;
														break;
													} else {
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
										if (CIGARite[1] + bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1).substring(0, CIGARite[1] + bias);
										}

										str3 = standMap.get(n).substring(1).substring(temSeqLen - CIGARite[2]);
									} else {
										if (CIGARite[1] + bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													CIGARite[1] + bias);
										}

										str3 = compRev.compRev(standMap.get(n).substring(1))
												.substring(temSeqLen - CIGARite[2]);
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
												int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
												if (!aligmentAno[0].equals(aligment[0])
														&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
														&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
													str2_ok = 1;
													break;
												} else {
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

						newNumSite = Integer.valueOf(aligment[2]) + CIGARite[3] - 1;
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
											if (CIGARite[2] - bias <= 0) {
												continue;
											} else {
												str = standMap.get(n).substring(1)
														.substring(temSeqLen - CIGARite[2] + bias);
											}
											str3 = standMap.get(n).substring(1).substring(0, CIGARite[1]);
										} else {
											if (CIGARite[2] - bias <= 0) {
												continue;
											} else {
												str = compRev.compRev(standMap.get(n).substring(1))
														.substring(temSeqLen - CIGARite[2] + bias);
											}

											str3 = compRev.compRev(standMap.get(n).substring(1)).substring(0,
													CIGARite[1]);
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
													int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
													if (!aligmentAno[0].equals(aligment[0])
															&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
															&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
														str2_ok = 1;
														break;
													} else {
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
										if (CIGARite[2] - bias <= 0) {
											continue;
										} else {
											str = standMap.get(n).substring(1)
													.substring(temSeqLen - CIGARite[2] + bias);
										}

										str3 = standMap.get(n).substring(1).substring(0, CIGARite[1]);
									} else {
										if (CIGARite[2] - bias <= 0) {
											continue;
										} else {
											str = compRev.compRev(standMap.get(n).substring(1))
													.substring(temSeqLen - CIGARite[2] + bias);
										}
										str3 = compRev.compRev(standMap.get(n).substring(1)).substring(0,CIGARite[1]);
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
												int[] CIGARiteAno = misd.misd(aligmentAno[4],temSeqLen);
												if (!aligmentAno[0].equals(aligment[0])
														&& Integer.valueOf(aligmentAno[2]) >= Integer.valueOf(temChrSite[0]) - 6
														&& Integer.valueOf(aligmentAno[2]) + CIGARiteAno[3]- 1 <= Integer.valueOf(temChrSite[1]) + 6) {
													str2_ok = 1;
													break;
												} else {
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

					int newNumSite = Integer.valueOf(aligment[2]);
					int startTem = newNumSite + 6;
					int endTem = newNumSite + CIGARite[3] - 7;
					temFSJId.addAll(GFC.getFSJ(startTem, endTem,aligment[1],CIGARite[0],siteArray1,siteArray2,siteMap1,siteMap2));
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
