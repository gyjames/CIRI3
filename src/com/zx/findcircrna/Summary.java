package com.zx.findcircrna;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import smith.SmithWaterman;

public class Summary {
	int strigency;
	SmithWaterman aligner = new SmithWaterman( 1, -1, -3);
	HashMap<String, String> circIdMap = new HashMap<String, String>();
	HashMap<String, String> chrTCGAMap;
	public Summary(int strigency,HashMap<String, String> chrTCGAMap) {
		super();
		this.strigency = strigency;
		this.chrTCGAMap = chrTCGAMap;
	}	
	public HashMap<String, String> getCircMap() {
		return circIdMap;
	}
	public ArrayList<String> summary(ArrayList<String> filePathList,HashMap<String, Integer> fileSplitNumMap,HashMap<String, Integer> circFSJNewMap,String UserGivecircRNA)
			throws IOException {
		ArrayList<String> SummaryCircList = new ArrayList<String>();
		HashMap<String, HashSet<String>> circMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSet = new HashSet<String>();	
		HashMap<String, HashSet<String>> circStartMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> circEndMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circStartSet = new HashSet<String>();
		HashSet<String> circEndSet = new HashSet<String>();
		//统计第一遍扫描的结果		
		for (String samFile : filePathList) {
			int threads = fileSplitNumMap.get(samFile);
			for (int i = 1; i <= threads; i++) {
				BufferedReader BSJBr = new BufferedReader(new FileReader(new File(samFile+"BSJ"+i)));	
				String line = BSJBr.readLine();
				while (line != null ) {					
					String[] circLineArr = line.split("\t",8);
					String chrStartEnd = circLineArr[3] + "\t" + circLineArr[4] + "\t" + circLineArr[5];
					if (!circMap.containsKey(chrStartEnd)) {
						circSet = new HashSet<String>();
						circSet.add(line);
						circMap.put(chrStartEnd, circSet);
					} else {
						circSet = circMap.get(chrStartEnd);
						circSet.add(line);
						circMap.put(chrStartEnd, circSet);
					}
					//strat
					String chrStart = circLineArr[3] + "\t" + circLineArr[4];
					if (!circStartMap.containsKey(chrStart)) {
						circStartSet = new HashSet<String>();
						circStartSet.add(chrStartEnd);
						circStartMap.put(chrStart, circStartSet);
					} else {
						circStartSet = circStartMap.get(chrStart);
						circStartSet.add(chrStartEnd);
						circStartMap.put(chrStart, circStartSet);
					}
					line = BSJBr.readLine();					
				}
				BSJBr.close();
				new File(samFile+"BSJ"+i).delete();	
			}
		}		
		String[] circRNAArr = null;	
		String chr = "";
		for (String chrStart : circStartMap.keySet()) {
			chr = chrStart.split("\t")[0];
			circStartSet = circStartMap.get(chrStart);	
			if (circStartSet.size()>1) {
				int[] score = new int[circStartSet.size()];
				int[] endSite = new int[circStartSet.size()];
				int[] CIGARS = new int[circStartSet.size()];
				String[] circNameArr = new String[circStartSet.size()];
				int index = 0;
				for (String chrStartEnd : circStartSet) {
					circSet = circMap.get(chrStartEnd);
					circNameArr[index] = chrStartEnd;
					int len = 0;
					for (String circRNA : circSet) {
						circRNAArr = circRNA.split("\t");
						if(circRNAArr[1].indexOf(";")>=0) {
							endSite[index] = Integer.valueOf(circRNAArr[5]);
							if (circRNAArr[2].equals("1") && circRNAArr[9].equals("1")) {
								score[index]++;
						    }
							String[] CIGARArr = circRNAArr[1].split(";");
							String temCIGAR = CIGARArr[1].replaceAll("[^(a-zA-Z)]", "");
							String[] CIGARNumArr = CIGARArr[1].split("S|M|D|I");
							int matchNum = 0;
							for (int k = temCIGAR.indexOf("M"); k <= temCIGAR.lastIndexOf("M"); k++) {
								if(temCIGAR.substring(k,k+1).equals("M") || temCIGAR.substring(k,k+1).equals("D")) {
									matchNum += Integer.valueOf(CIGARNumArr[k]);
								}										
							}
							if (matchNum > len) {
								len = matchNum;
							}
						}				
					}
					CIGARS[index] = len;
					index++;
				}
				ArrayList<Integer> upScoreList = new ArrayList<Integer>();
				ArrayList<Integer> zeroScoreList = new ArrayList<Integer>();
				for (int i = 0; i < score.length; i++) {
					if (score[i]>0) {
						upScoreList.add(i);
					}else{
						zeroScoreList.add(i);
					}
				}
				ArrayList<Integer> diffList = new ArrayList<>();
				for (Integer indexUp : upScoreList) {
					diffList = new ArrayList<>();
					for (Integer indexZero : zeroScoreList) {
						String chrTCGA = chrTCGAMap.get(chr);
						String seq1 = chrTCGA.substring(endSite[indexUp]-CIGARS[indexZero],endSite[indexUp]);
						String seq2 = chrTCGA.substring(endSite[indexZero]-CIGARS[indexZero],endSite[indexZero]);
						aligner.setSeq(seq1,seq2);
						String[] alignment = aligner.getAlignment();
						if (aligner.getAlignmentScore() >= seq1.length() - 2 -(seq1.length()-1)/10*2  && alignment[1].length() > seq1.length() - 2) {
							
							HashSet<String> circSet1 = circMap.get(circNameArr[indexUp]);	
							HashSet<String> circSet2 = circMap.get(circNameArr[indexZero]);
							circSet1.addAll(circSet2);
							circMap.put(circNameArr[indexUp], circSet1);
							circMap.remove(circNameArr[indexZero]);
							diffList.add(indexZero);
						}
					}
					zeroScoreList.removeAll(diffList);
				}
			}
		}
		
		for (String chrStartEnd : circMap.keySet()) {
			String[] chrStartEndArr = chrStartEnd.split("\t");
			//end
			String chrEnd = chrStartEndArr[0] + "\t" + chrStartEndArr[2];
			if (!circEndMap.containsKey(chrEnd)) {
				circEndSet = new HashSet<String>();
				circEndSet.add(chrStartEnd);
				circEndMap.put(chrEnd, circEndSet);
			} else {
				circEndSet = circEndMap.get(chrEnd);
				circEndSet.add(chrStartEnd);
				circEndMap.put(chrEnd, circEndSet);
			}
		}
		
		for (String chrEnd : circEndMap.keySet()) {
			chr = chrEnd.split("\t")[0];
			circEndSet = circEndMap.get(chrEnd);
			if (circEndSet.size()>1) {
				int[] score = new int[circEndSet.size()];
				int[] startSite = new int[circEndSet.size()];
				int[] CIGARS = new int[circEndSet.size()];
				String[] circNameArr = new String[circEndSet.size()];
				int index = 0;
				for (String chrStartEnd : circEndSet) {
					circSet = circMap.get(chrStartEnd);
					circNameArr[index] = chrStartEnd;
					int len = 0;
					for (String circRNA : circSet) {
						circRNAArr = circRNA.split("\t");
						if(circRNAArr[1].indexOf(";")>=0) {
							startSite[index] = Integer.valueOf(circRNAArr[4]);
							if (circRNAArr[2].equals("1") && circRNAArr[9].equals("1")) {
								score[index]++;
						    }
							String[] CIGARArr = circRNAArr[1].split(";");
							String temCIGAR = CIGARArr[0].replaceAll("[^(a-zA-Z)]", "");
							String[] CIGARNumArr = CIGARArr[0].split("S|M|D|I");
							int matchNum = 0;
							for (int k = temCIGAR.indexOf("M"); k <= temCIGAR.lastIndexOf("M"); k++) {
								if(temCIGAR.substring(k,k+1).equals("M") || temCIGAR.substring(k,k+1).equals("D")) {
									matchNum += Integer.valueOf(CIGARNumArr[k]);
								}										
							}
							if (matchNum > len) {
								len = matchNum;
							}
						}							
					}
					CIGARS[index] = len;
					index++;
					
				}
				ArrayList<Integer> upScoreList = new ArrayList<Integer>();
				ArrayList<Integer> zeroScoreList = new ArrayList<Integer>();
				for (int i = 0; i < score.length; i++) {
					if (score[i]>0) {
						upScoreList.add(i);
					}else{
						zeroScoreList.add(i);
					}
				}
				ArrayList<Integer> diffList = new ArrayList<>();
				for (Integer indexUp : upScoreList) {
					diffList = new ArrayList<>();
					for (Integer indexZero : zeroScoreList) {
						String chrTCGA = chrTCGAMap.get(chr);
						String seq1 = chrTCGA.substring(startSite[indexUp]-CIGARS[indexZero],startSite[indexUp]);
						String seq2 = chrTCGA.substring(startSite[indexZero]-CIGARS[indexZero],startSite[indexZero]);
						aligner.setSeq(seq1,seq2);
						String[] alignment = aligner.getAlignment();
						if (aligner.getAlignmentScore() >= seq1.length() - 2 -(seq1.length()-1)/10*2 && alignment[1].length() > seq1.length() - 2) {
							HashSet<String> circSet1 = circMap.get(circNameArr[indexUp]);	
							HashSet<String> circSet2 = circMap.get(circNameArr[indexZero]);
							circSet1.addAll(circSet2);
							circMap.put(circNameArr[indexUp], circSet1);
							circMap.remove(circNameArr[indexZero]);
							diffList.add(indexZero);
						}
					}
					zeroScoreList.removeAll(diffList);
				}
			}
		}
		
		HashSet<String> circIdSet3 = new HashSet<String>();
     	HashSet<String> CIGARSet3 = new HashSet<String>();
		HashSet<String> falseCIGARSet3 = new HashSet<String>();
		HashSet<String> circIdSet = new HashSet<String>();
     	HashSet<String> CIGARSet = new HashSet<String>();
		HashSet<String> falseCIGARSet = new HashSet<String>();
		
		
		for (String chrStartEnd : circMap.keySet()) {
			String[] chrStartEndArr = chrStartEnd.split("\t");
			circSet = circMap.get(chrStartEnd);
			int TPReads = 0, nonReads = 0, FPReads = 0,TPReads3 = 0, nonReads3 = 0, FPReads3 = 0,tag = 0;
			Integer[] CIGARCount3 = { 0, 0, 0 };
			CIGARSet.clear();
			circIdSet.clear();
			falseCIGARSet.clear();
			CIGARSet3.clear();
			circIdSet3.clear();
			falseCIGARSet3.clear();
			for (String circRNA : circSet) {
				circRNAArr = circRNA.split("\t");
				if (circRNAArr[2].equals("1")) {
					if(circRNAArr[9].equals("3")) {							
						CIGARSet3.add(circRNAArr[1]);							
						circIdSet3.add(circRNAArr[0]);							
					}else {
						String[] CIGARArr = circRNAArr[1].split(";");
						if (circRNAArr[9].equals("1")) {
							tag++;
						}
						for (int i = 0; i < CIGARArr.length; i++) {
							CIGARSet.add(CIGARArr[i]);
						}
						circIdSet.add(circRNAArr[0]);
					}
				} else if (circRNAArr[2].equals("-1")) {						
					if(circRNAArr[9].equals("3")) {							
						falseCIGARSet3.add(circRNAArr[1]);							
						nonReads3++;							
					}else {
						String[] CIGARArr = circRNAArr[1].split(";");
						for (int i = 0; i < CIGARArr.length; i++) {
							falseCIGARSet.add(CIGARArr[i]);
						}
						nonReads++;
					}					
					
				} else if (circRNAArr[2].equals("-2")) {					
					if(circRNAArr[9].equals("3")) {
						falseCIGARSet3.add(circRNAArr[1]);
						FPReads3++;
					}else {
						String[] CIGARArr = circRNAArr[1].split(";");
						for (int i = 0; i < CIGARArr.length; i++) {
							falseCIGARSet.add(CIGARArr[i]);
						}
						FPReads++;							
					}
					
				}

			}				
			CIGARSet3.addAll(CIGARSet);
			circIdSet3.addAll(circIdSet);
			falseCIGARSet3.addAll(falseCIGARSet);
			
			for (String site : CIGARSet3) {
				String temCIGAR = site.replaceAll("[^(a-zA-Z)]", "");
				if (temCIGAR.equals("M")) {
				} else if (temCIGAR.indexOf("M") > 0 && temCIGAR.lastIndexOf("M") < temCIGAR.length() - 1) {
					CIGARCount3[2] += 1;
				} else if (temCIGAR.indexOf("M") == 0) {
					CIGARCount3[1] += 1;
				} else {
					CIGARCount3[0] += 1;
				}
			}				
			TPReads = circIdSet.size();
			TPReads3 = circIdSet3.size();
			FPReads3 += FPReads;
			nonReads3 += nonReads;
			if (strigency == 2) {						
				if (((TPReads > 19 * FPReads || FPReads <= 1) && TPReads > nonReads + FPReads && CIGARSet.size() >= 3 && TPReads >=2) ||
						(tag > 0 &&  falseCIGARSet3.size() == 0 && CIGARSet3.size() >= 3 && TPReads3 >=2)) {
					SummaryCircList.add(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2] + "\t" + chrStartEnd
							+ "\t" + TPReads3 + "\t" + CIGARCount3[0] + "_" + CIGARCount3[1] + "_" + CIGARCount3[2] + "\t"
							+ circFSJNewMap.get(chrStartEnd) + "\t"
							+ String.format("%.2f", (TPReads3 * 2 / (TPReads3 * 2 +circFSJNewMap.get(chrStartEnd) + 0.0)) ) + "\t" + circRNAArr[6]
							+ "\t" + String.join(",", circIdSet3)+"\t"+tag);
					circIdMap.put(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2], "");
				}
				} else if (strigency == 1) {						
						if (((TPReads > 19 * FPReads || falseCIGARSet.size() <= 2) && TPReads > nonReads + FPReads && TPReads >= 2) ||
								(tag > 0 &&  falseCIGARSet3.size() == 0 && TPReads3 >= 2)) {
							SummaryCircList.add(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2] + "\t" + chrStartEnd
									+ "\t" + TPReads3 + "\t" + CIGARCount3[0] + "_" + CIGARCount3[1] + "_" + CIGARCount3[2] + "\t"
									+ circFSJNewMap.get(chrStartEnd) + "\t"
									+ String.format("%.2f", (TPReads3 * 2 / (TPReads3 * 2 + circFSJNewMap.get(chrStartEnd) + 0.0))) + "\t" + circRNAArr[6]
									+ "\t" + String.join(",", circIdSet3)+"\t"+tag);
							circIdMap.put(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2], "");
						}
				} else if (strigency == 0 ) {					
						if (((TPReads > 19 * FPReads || falseCIGARSet.size() <= 2)&& TPReads > nonReads + FPReads) ||
								tag > 0 &&  falseCIGARSet3.size() == 0  && TPReads3 >= 2) {
							SummaryCircList.add(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2] + "\t" + chrStartEnd
									+ "\t" + TPReads3 + "\t" + CIGARCount3[0] + "_" + CIGARCount3[1] + "_" + CIGARCount3[2] + "\t"
									+ circFSJNewMap.get(chrStartEnd) + "\t"
									+ String.format("%.2f", (TPReads3 * 2 / (TPReads3 * 2 + circFSJNewMap.get(chrStartEnd) + 0.0))) + "\t" + circRNAArr[6]
									+ "\t" + String.join(",", circIdSet3)+"\t"+tag);
							circIdMap.put(chrStartEndArr[0] + ":" + chrStartEndArr[1] + "|" + chrStartEndArr[2], "");
						}
				}
		}			
      return SummaryCircList;
	}

}
