package com.zx.hg38;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import com.zx.findcircrna.SiteSort;

public class AnnotationIntron {
	public void annotation(ArrayList<String> SummaryCircList,HashMap<String, ArrayList<SiteSort>> geneExonMap,
			HashMap<String, ArrayList<String>> chrExonStartTranscriptMap,HashMap<String, ArrayList<String>> chrExonEndTranscriptMap,
			HashMap<String, ArrayList<Integer[]>> exonListMap,
			String outPut,int relExp) throws IOException {
		HashMap<String, ArrayList<SiteSort>> circChrMap = new HashMap<String, ArrayList<SiteSort>>();
		ArrayList<SiteSort> circList = new ArrayList<SiteSort>();
		for (String Circ : SummaryCircList) {
			String[] circArr = Circ.split("\t");
			if (Double.valueOf(circArr[7]) < relExp) {
				continue;
			}
			if(!circChrMap.containsKey(circArr[1])) {
				circList = new ArrayList<SiteSort>();
				circList.add(new SiteSort(Integer.valueOf(circArr[2]),circArr));
				circChrMap.put(circArr[1], circList);
			}else {
				circList = circChrMap.get(circArr[1]);
				circList.add(new SiteSort(Integer.valueOf(circArr[2]),circArr));
				circChrMap.put(circArr[1], circList);
			}	
		}
        ////////排序
		for (String chrKey : circChrMap.keySet()) {
			circList = circChrMap.get(chrKey);
			Collections.sort(circList);
			circChrMap.put(chrKey, circList);
		}
		ArrayList<SiteSort> geneSiteList = new ArrayList<SiteSort>();
		ArrayList<Integer[]> exonSiteList = new ArrayList<Integer[]>();
		BufferedWriter annotationCirc = new BufferedWriter(new FileWriter(new File(outPut)));
		annotationCirc.write("circRNA_ID"+ "\t"+"chr"+ "\t"+"circRNA_start"+ "\t"+"circRNA_end"+ "\t"+"#junction_reads"+ "\t"+"SM_MS_SMS"
				+ "\t"+"#non_junction_reads"+ "\t"+"junction_reads_ratio"+ "\t"+"circRNA_type"+ "\t"+"gene_id"+ "\t"+"strand"+ "\t"+ "junction_reads_ID"+"\n");
		for (String chrKey : circChrMap.keySet()) {
			circList = circChrMap.get(chrKey);
			if (geneExonMap.containsKey(chrKey)) {
				geneSiteList = geneExonMap.get(chrKey);
				int tag = 0;
				for (int i = 0; i < circList.size(); i++) {
					int tagType = 0;
					    String[] circTem = circList.get(i).getLength();
						String startKey = circTem[1] +"\t"+ circTem[2];
						String endKey = circTem[1] +"\t"+ circTem[3];
						//判断起点和终点是否是注释文件中的exon起点和终点
						if (chrExonStartTranscriptMap.containsKey(startKey) && chrExonEndTranscriptMap.containsKey(endKey)) {
							String[] geneStandStartArr = chrExonStartTranscriptMap.get(startKey).get(0).split("\t");
							String[] geneStandEndArr = chrExonEndTranscriptMap.get(endKey).get(0).split("\t");
							if (geneStandStartArr[0].equals(geneStandEndArr[0])) {
								annotationCirc.write(circTem[0]+ "\t"+circTem[1]+ "\t"+circTem[2]+ "\t"+circTem[3]+ "\t"+circTem[4]+ "\t"+circTem[5]+ "\t"+circTem[6]+ "\t"+circTem[7]+ "\t"+
										"exon"+ "\t"+geneStandStartArr[0]+ "\t"+circTem[8]+ "\t"+circTem[9]+ "\n");
								continue;
							}
						}
					
					for (int j = tag; j < geneSiteList.size(); j++) {
						if (circList.get(i).getSite()>Integer.valueOf(geneSiteList.get(j).getLength()[3])) {
							
						}else {
							tag =j;
							break;
						}
					}
					for (int j = tag; j < geneSiteList.size(); j++) {
						if (circList.get(i).getSite()>=geneSiteList.get(j).getSite()) {
							if (Integer.valueOf(circTem[3])<=Integer.valueOf(geneSiteList.get(j).getLength()[3])) {
								String[] geneTem =  geneSiteList.get(j).getLength();
								exonSiteList = exonListMap.get(geneTem[0]+ "\t"+geneTem[1]);
								int if_start_ok = 0,if_end_ok = 0;
								//System.out.println(geneTem[0]+ "\t"+geneTem[1]+"\t"+exonSiteList.size());
								for (int k = 0; k < exonSiteList.size(); k++) {
									Integer[] exonTem= exonSiteList.get(k);
									if (exonTem[0]<= Integer.valueOf(circTem[2]) && exonTem[1]>= Integer.valueOf(circTem[2])) {
										if_start_ok = 1;
									}
									if (exonTem[0]<= Integer.valueOf(circTem[3]) && exonTem[1]>= Integer.valueOf(circTem[3])) {
										if_end_ok = 1;
									}
								}
								if (if_start_ok ==1 && if_end_ok ==1) {
									annotationCirc.write(circTem[0]+ "\t"+circTem[1]+ "\t"+circTem[2]+ "\t"+circTem[3]+ "\t"+circTem[4]+ "\t"+circTem[5]+ "\t"+circTem[6]+ "\t"+circTem[7]+ "\t"+
											"exon"+ "\t"+geneTem[1]+ "\t"+circTem[8]+ "\t"+circTem[9]+ "\n");
									tagType = 1;
									break;
								}else {
									annotationCirc.write(circTem[0]+ "\t"+circTem[1]+ "\t"+circTem[2]+ "\t"+circTem[3]+ "\t"+circTem[4]+ "\t"+circTem[5]+ "\t"+circTem[6]+ "\t"+circTem[7]+ "\t"+
											"intron"+ "\t"+geneTem[1]+ "\t"+circTem[8]+ "\t"+circTem[9]+ "\n");
									tagType = 1;
									break;
								}
								
								
							}
						}else {
							
							break;
						}
					}
					if(tagType == 0) {
						circTem = circList.get(i).getLength();
						annotationCirc.write(circTem[0]+ "\t"+circTem[1]+ "\t"+circTem[2]+ "\t"+circTem[3]+ "\t"+circTem[4]+ "\t"+circTem[5]+ "\t"+circTem[6]+ "\t"+circTem[7]+ "\t"+
								"intergenic_region"+ "\t"+"NA"+ "\t"+circTem[8]+ "\t"+circTem[9]+ "\n");
					}
					
				}
			}else {
				for (int i = 0; i < circList.size(); i++) {
					String[] circTem = circList.get(i).getLength();
					annotationCirc.write(circTem[0]+ "\t"+circTem[1]+ "\t"+circTem[2]+ "\t"+circTem[3]+ "\t"+circTem[4]+ "\t"+circTem[5]+ "\t"+circTem[6]+ "\t"+circTem[7]+ "\t"+
							"intergenic_region"+ "\t"+"NA"+ "\t"+circTem[8]+ "\t"+circTem[9]+ "\n");
				}
				
			}

		}
		annotationCirc.close();
	}
}
