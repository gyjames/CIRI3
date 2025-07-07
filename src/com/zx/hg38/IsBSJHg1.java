package com.zx.hg38;

import java.util.ArrayList;
import java.util.HashMap;

public class IsBSJHg1 {
	int linear_range_size_min,initial_size1 = 5 + 2;
	String mitochondrion;
	boolean spLable;
	IsInCircRNA1_2 IIC1_2 = new IsInCircRNA1_2();
	IsInCircRNA2 IIC2 = new IsInCircRNA2();
	IsInCircRNA3 IIC3 = new IsInCircRNA3();
	IndexCompare indexCompare = new IndexCompare();
	HashMap<String, String> chrExonStartMap,chrExonEndMap;
	String[] shiftArr = null;
	String[] str_new = new String[2];
	public IsBSJHg1(int linear_range_size_min,HashMap<String, String> chrExonStartMap,HashMap<String, String> chrExonEndMap,
			HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,HashMap<String, ArrayList<String>> chrExonEndTranscriptMap,
			String mitochondrion,boolean spLable) {
		super();
		this.linear_range_size_min = linear_range_size_min;
		this.chrExonStartMap = chrExonStartMap;
		this.chrExonEndMap = chrExonEndMap;
		this.mitochondrion = mitochondrion;
		this.spLable = spLable;
	}
	public String isBSJHg1(String circLineArr[], String chrTAGA,int sumQ){
			String end_string1, end_string2, str_adjustment, initial_seq1, initial_seq2,
					circ_range_seq = null, linear_range = null;
			int tmp_site1, tmp_site2, adjt_bp, diff_adjt, site1_new = 0, site2_new = 0;				
				//////////////////////////
				int site1 = Integer.valueOf(circLineArr[9]);
				int site2 = Integer.valueOf(circLineArr[10]);
				int end_adjt1 = Integer.valueOf(circLineArr[11]);
				int end_adjt2 = Integer.valueOf(circLineArr[12]);
				int total_adjustment = end_adjt1 + end_adjt2;
				int chrTAGALen = chrTAGA.length();
				// Extract end strings from reference genome according to putative BSJ as well
				// as adjustment
				if (end_adjt2 >= 0) {
					tmp_site1 = site1 - end_adjt1 - 1;
					tmp_site2 = site2 - end_adjt1 - 1;
					adjt_bp = 2 + end_adjt1 + end_adjt2;
					if (site1 - end_adjt1 - 4 >= 0) {
						end_string1 = chrTAGA.substring(site1 - end_adjt1 - 4, end_adjt2 + site1);
					} else {
						StringBuffer remChart = new StringBuffer();
						for (int n = site1 - end_adjt1 - 4; n < 0; n++) {
							remChart.append("N");
						}
						end_string1 = remChart.toString() + chrTAGA.substring(0, end_adjt2 + site1);
					}
					if (3 + end_adjt2 + site2 > chrTAGALen) {
						end_string2 = chrTAGA.substring(site2 - end_adjt1 - 1, chrTAGALen);
					} else {
						end_string2 = chrTAGA.substring(site2 - end_adjt1 - 1, 3 + end_adjt2 + site2);
					}

				} else {
					tmp_site1 = site1 + end_adjt2 - 1;
					tmp_site2 = site2 + end_adjt2 - 1;
					adjt_bp = 2 - end_adjt1 - end_adjt2;
					if (site1 + end_adjt2 - 4 >= 0) {
						end_string1 = chrTAGA.substring(site1 + end_adjt2 - 4, site1 - end_adjt1);
					} else {
						StringBuffer remChart = new StringBuffer();
						for (int n = site1 + end_adjt2 - 4; n < 0; n++) {
							remChart.append("N");
						}
						end_string1 = remChart.toString() + chrTAGA.substring(0, site1 - end_adjt1);
					}
					if (3 + site2 - end_adjt1 > chrTAGALen) {
						end_string2 = chrTAGA.substring(site2 + end_adjt2 - 1, chrTAGALen);
					} else {
						end_string2 = chrTAGA.substring(site2 + end_adjt2 - 1, 3 + site2 - end_adjt1);
					}

				}
				// Pass end strings to index_compare to find possible splicing signals and the
				// corresponding strand	
				HashMap<Integer, String> indexStrandMap;
				if (circLineArr[1].equals(mitochondrion) || spLable) {
					indexStrandMap = indexCompare.indexCompareChrM(end_string1, end_string2);
				}else {
					indexStrandMap = indexCompare.indexCompare(end_string1, end_string2);
				}				
				for (int i = 0; i <= adjt_bp; i++) {
					String startKey = circLineArr[1] +"\t"+ (tmp_site1 + i);
					String endKey = circLineArr[1] +"\t"+ (tmp_site2 + i);
					if (!indexStrandMap.containsKey(i) && chrExonStartMap.containsKey(startKey)
							&& chrExonEndMap.containsKey(endKey)) {
						String geneStand = chrExonStartMap.get(startKey);
						if (geneStand.equals(chrExonEndMap.get(endKey))) {
							String[] geneStandArr = geneStand.split("\t");
							indexStrandMap.put(i,
									i+ "\t"+geneStandArr[1]+ "\t"+ chrTAGA.substring(tmp_site1 + i - 3, tmp_site1 + i - 1) + "\t"
											+ chrTAGA.substring(tmp_site2 + i, tmp_site2 + i + 2));
						}
					}
				}
				if (indexStrandMap.keySet().size() > 0) {
					// 1:strand 2:chr 3:str1 4:str2 5:str3 6:str4 7:0 8:1 9:str4_ok 10:site1
					// 11:site2 12:end_adjt1 13:end_adjt2 45H31M 46M30S
					int junc_ok = 0;
					for (Integer shift : indexStrandMap.keySet()) {
						shiftArr = indexStrandMap.get(shift).split("\t");
						if (end_adjt2 >= 0) {
							diff_adjt = Integer.valueOf(shiftArr[0]) - 1 - end_adjt1;
						} else {
							diff_adjt = Integer.valueOf(shiftArr[0]) - 1 + total_adjustment - end_adjt1;
						}
						site1_new = site1 + diff_adjt;
						site2_new = site2 + diff_adjt;
						if (diff_adjt >= 0) {
							str_adjustment = circLineArr[2].substring(0, diff_adjt);
							str_new[1] = circLineArr[3] + str_adjustment;
							str_new[0] = circLineArr[2].substring(diff_adjt);
						} else if (diff_adjt < 0) {
							str_adjustment = circLineArr[3].substring(circLineArr[3].length() + diff_adjt);
							str_new[0] = str_adjustment + circLineArr[2];
							str_new[1] = circLineArr[3].substring(0, circLineArr[3].length() + diff_adjt);
						}
						str_new[0] = shiftArr[2] + str_new[0];
						str_new[1] = str_new[1] + shiftArr[3];
						initial_seq1 = str_new[0].substring(0, initial_size1);
						initial_seq2 = str_new[1].substring(str_new[1].length() - initial_size1,
								str_new[1].length());

						if (site1_new - 3 < 0 && site2_new + 2 > chrTAGALen) {
							circ_range_seq = chrTAGA.substring(0, chrTAGALen);
						} else if (site1_new - 3 < 0) {
							circ_range_seq = chrTAGA.substring(0, site2_new + 2);
						} else if (site2_new + 2 > chrTAGALen) {
							circ_range_seq = chrTAGA.substring(site1_new - 3, chrTAGALen);
						} else {
							circ_range_seq = chrTAGA.substring(site1_new - 3, site2_new + 2);
						}
						if (circ_range_seq.substring(0, initial_seq1.length()).equals(initial_seq1)
								&& circ_range_seq.substring(circ_range_seq.length() - initial_seq2.length(),
										circ_range_seq.length()).equals(initial_seq2)) {
			
							junc_ok = 1;
							for (int i = 0; i <= 1; i++) {
								if (Integer.valueOf(circLineArr[6 + i]) != 1) {
									int len_str = str_new[i].length();
									if (i == 1 && site2_new - site1_new + 5 >= linear_range_size_min) {
										if (2 * site1_new >= site2_new + 6) {
											linear_range = chrTAGA.substring(2 * site1_new - site2_new - 6,
													site1_new - 1);
										} else {
											linear_range = chrTAGA.substring(0, site1_new - 1);
										}
									} else if (i == 1) {
										if (site1_new >= linear_range_size_min + 1) {
											linear_range = chrTAGA.substring(
													site1_new - linear_range_size_min - 1, site1_new - 1);
										} else {
											linear_range = chrTAGA.substring(0, site1_new - 1);
										}
									} else if (i == 0 && site2_new - site1_new + 5 >= linear_range_size_min) {
										if (2 * site2_new - site1_new + 5 > chrTAGALen) {
											linear_range = chrTAGA.substring(site2_new, chrTAGALen);
										} else {
											linear_range = chrTAGA.substring(site2_new,
													2 * site2_new - site1_new + 5);
										}
									} else {
										if (site2_new + linear_range_size_min > chrTAGALen) {
											linear_range = chrTAGA.substring(site2_new, chrTAGALen);
										} else {
											linear_range = chrTAGA.substring(site2_new,
													site2_new + linear_range_size_min);
										}

									}
									// Comparison of detected seeds in the two regions
									// Seed is searched iteratively in the descending order of length
									circLineArr[6 + i] = IIC1_2.isInCircRNA1_2(len_str, str_new[i], circ_range_seq, linear_range)+"";
									
								}
							}
						}				
					if (circLineArr[6].equals("1") && circLineArr[7].equals("1") && junc_ok == 1) {						
						if (!circLineArr[4].equals("*")) {
							if (IIC2.isInCircRNA2(circLineArr[4], circ_range_seq) == 0) {
								return null;
							}							
						}
						// Paired end mapping signal is also detected by multiple seed matching
						int tag =1;
						if (circLineArr[5].length() > 5) {
							String pem_null_range_seq = "";
							if (circLineArr[0].equals("1")
									&& site1_new - site1_new + 5 >= linear_range_size_min) {
								if (2 * site1_new >= site2_new + 6) {
									pem_null_range_seq = chrTAGA.substring(2 * site1_new - site2_new - 6,
											site1_new - 1);
								} else {
									pem_null_range_seq = chrTAGA.substring(0, site1_new - 1);
								}
							} else if (circLineArr[0].equals("1")) {
								if (site1_new >= linear_range_size_min + 1) {
									pem_null_range_seq = chrTAGA
											.substring(site1_new - linear_range_size_min - 1, site1_new - 1);
								} else {
									pem_null_range_seq = chrTAGA.substring(0, site1_new - 1);
								}
							} else if (circLineArr[0].equals("0")
									&& site2_new - site1_new + 5 > linear_range_size_min) {
								if (2 * site2_new - site1_new + 5 > chrTAGALen) {
									pem_null_range_seq = chrTAGA.substring(site2_new, chrTAGALen);
								} else {
									pem_null_range_seq = chrTAGA.substring(site2_new,
											2 * site2_new - site1_new + 5);
								}

							} else {
								if (site2_new + linear_range_size_min > chrTAGALen) {
									pem_null_range_seq = chrTAGA.substring(site2_new, chrTAGALen);
								} else {
									pem_null_range_seq = chrTAGA.substring(site2_new,
											site2_new + linear_range_size_min);
								}

							}
							tag = IIC3.isInCircRNA3(circLineArr[5],circLineArr[8],circ_range_seq,pem_null_range_seq);
											
						} 
						String circInfor =tag + "\t" + circLineArr[1]+"\t"+site1_new + "\t" + site2_new + "\t" + shiftArr[1]+"\t"+shiftArr[2]+ "\t" + shiftArr[3] 
								+"\t" +sumQ;
						return circInfor;
					    }
				    }
				}
				return null;			
	}

}
