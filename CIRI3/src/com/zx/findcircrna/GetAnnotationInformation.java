package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GetAnnotationInformation {
	// chr-exonStart:gene-stand
	HashMap<String, String> chrExonStartMap = new HashMap<String, String>();	
	// chr-exonEnd:gene-stand
	HashMap<String, String> chrExonEndMap = new HashMap<String, String>();
	// chr:{geneStart:chr-geneId-geneStart-geneEnd}
	HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
	ArrayList<SiteSort> geneSiteList = new ArrayList<SiteSort>();
	// chr-geneId:{exonStart, exonEnd}
	HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
	ArrayList<Integer[]> exonSiteList = new ArrayList<Integer[]>();
	
    //chr-exonStart:list{gene-stand-transcript}
	HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>();
	//chr-exonEnd:list{gene-stand-transcript}
	HashMap<String, ArrayList<String>> chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();
	ArrayList<String> geneStandTranscriptList1 = new ArrayList<String>();
	ArrayList<String> geneStandTranscriptList2 = new ArrayList<String>();
	
	public void hand(String annotationFile, boolean lable) throws IOException {
		String geneId = "", chr = "";
		int geneStart = 0, geneEnd = 0;
		BufferedReader annotationRead = new BufferedReader(new FileReader(new File(annotationFile)));
		String annotation = annotationRead.readLine();
		if (lable) {
			// 包含intron circRNA
			Pattern r = null;
			if (Pattern.matches(".*gff$", annotationFile) || Pattern.matches(".*gff3$", annotationFile)) {
				r = Pattern.compile("gene_id=(\\S+?);transcript_id=(\\S+?);");
			} else if (Pattern.matches(".*gtf$", annotationFile)) {
				r = Pattern.compile("gene_id \\\"(\\S+)\\\";[\\s\\S]+?transcript_id \\\"(\\S+)\\\"");
			}
			while (annotation != null) {
				if (annotation.startsWith("#")) {
					annotation = annotationRead.readLine();
					continue;
				}
				String[] annotationArr = annotation.split("\t");
				if (annotationArr[2].equals("exon")) {
					Matcher m = r.matcher(annotationArr[8]);
					if (m.find()) {
						if (!chrExonStartTranscriptMap.containsKey(annotationArr[0] + "\t" + annotationArr[3])) {
							geneStandTranscriptList1 = new ArrayList<String>();
							geneStandTranscriptList1.add(m.group(1) + "\t" + annotationArr[6] + "\t" + m.group(2));
							chrExonStartTranscriptMap.put(annotationArr[0] + "\t" + annotationArr[3],
									geneStandTranscriptList1);
						} else {
							geneStandTranscriptList1 = chrExonStartTranscriptMap
									.get(annotationArr[0] + "\t" + annotationArr[3]);
							geneStandTranscriptList1.add(m.group(1) + "\t" + annotationArr[6] + "\t" + m.group(2));
							chrExonStartTranscriptMap.put(annotationArr[0] + "\t" + annotationArr[3],
									geneStandTranscriptList1);
						}

						if (!chrExonEndTranscriptMap.containsKey(annotationArr[0] + "\t" + annotationArr[4])) {
							geneStandTranscriptList2 = new ArrayList<String>();
							geneStandTranscriptList2.add(m.group(1) + "\t" + annotationArr[6] + "\t" + m.group(2));
							chrExonEndTranscriptMap.put(annotationArr[0] + "\t" + annotationArr[4],
									geneStandTranscriptList2);
						} else {
							geneStandTranscriptList2 = chrExonEndTranscriptMap
									.get(annotationArr[0] + "\t" + annotationArr[4]);
							geneStandTranscriptList2.add(m.group(1) + "\t" + annotationArr[6] + "\t" + m.group(2));
							chrExonEndTranscriptMap.put(annotationArr[0] + "\t" + annotationArr[4],
									geneStandTranscriptList2);
						}

						if (!geneId.equals(m.group(1))) {
							String exonListMapKey = chr + "\t" + geneId;
							exonListMap.put(exonListMapKey, exonSiteList);

							exonSiteList = new ArrayList<Integer[]>();
							String[] geneTem = { chr, geneId, geneStart + "", geneEnd + "" };

							if (!geneExonMap.containsKey(chr)) {
								geneSiteList = new ArrayList<SiteSort>();
								geneSiteList.add(new SiteSort(geneStart, geneTem));
								geneExonMap.put(chr, geneSiteList);
							} else {
								geneSiteList = geneExonMap.get(chr);
								geneSiteList.add(new SiteSort(geneStart, geneTem));
								geneExonMap.put(chr, geneSiteList);
							}
							geneId = m.group(1);
							chr = annotationArr[0];
							geneStart = Integer.valueOf(annotationArr[3]);
							geneEnd = Integer.valueOf(annotationArr[4]);
							Integer[] temArr = { geneStart, geneEnd };
							exonSiteList.add(temArr);
						} else {
							int exonStartTem = Integer.valueOf(annotationArr[3]);
							int exonEndTem = Integer.valueOf(annotationArr[4]);
							Integer[] temArr = { exonStartTem, exonEndTem };
							exonSiteList.add(temArr);
							if (geneStart > exonStartTem) {
								geneStart = exonStartTem;
							}
							if (geneEnd < exonEndTem) {
								geneEnd = exonEndTem;
							}
						}

					}
				}
				annotation = annotationRead.readLine();
			}
			String exonListMapKey = chr + "\t" + geneId;
			exonListMap.put(exonListMapKey, exonSiteList);
			String[] geneTem = { chr, geneId, geneStart + "", geneEnd + "" };
			if (!geneExonMap.containsKey(chr)) {
				geneSiteList = new ArrayList<SiteSort>();
				geneSiteList.add(new SiteSort(geneStart, geneTem));
				geneExonMap.put(chr, geneSiteList);
			} else {
				geneSiteList = geneExonMap.get(chr);
				geneSiteList.add(new SiteSort(geneStart, geneTem));
				geneExonMap.put(chr, geneSiteList);
			}
		} else {
			// 不包含intron circRNA
			Pattern r = null;
			if (Pattern.matches(".*gff$", annotationFile) || Pattern.matches(".*gff3$", annotationFile)) {
				r = Pattern.compile("gene_id=(\\S+?);");
			} else if (Pattern.matches(".*gtf$", annotationFile)) {
				r = Pattern.compile("gene_id \\\"(\\S+)\\\"");
			}
			while (annotation != null) {
				if (annotation.startsWith("#")) {
					annotation = annotationRead.readLine();
					continue;
				}
				String[] annotationArr = annotation.split("\t");
				if (annotationArr[2].equals("exon")) {
					Matcher m = r.matcher(annotationArr[8]);
					if (m.find()) {
						// chr exonStart:gene stand
						chrExonStartMap.put(annotationArr[0] + "\t" + annotationArr[3],
								m.group(1) + "\t" + annotationArr[6]);
						// chr exonEnd:gene stand
						chrExonEndMap.put(annotationArr[0] + "\t" + annotationArr[4],
								m.group(1) + "\t" + annotationArr[6]);
						if (!geneId.equals(m.group(1))) {
							String exonListMapKey = chr + "\t" + geneId;
							exonListMap.put(exonListMapKey, exonSiteList);
							exonSiteList = new ArrayList<Integer[]>();
							String[] geneTem = { chr, geneId, geneStart + "", geneEnd + "" };
							if (!geneExonMap.containsKey(chr)) {
								geneSiteList = new ArrayList<SiteSort>();
								geneSiteList.add(new SiteSort(geneStart, geneTem));
								geneExonMap.put(chr, geneSiteList);
							} else {
								geneSiteList = geneExonMap.get(chr);
								geneSiteList.add(new SiteSort(geneStart, geneTem));
								geneExonMap.put(chr, geneSiteList);
							}
							geneId = m.group(1);
							chr = annotationArr[0];
							geneStart = Integer.valueOf(annotationArr[3]);
							geneEnd = Integer.valueOf(annotationArr[4]);
							Integer[] temArr = { geneStart, geneEnd };
							exonSiteList.add(temArr);
						} else {
							int exonStartTem = Integer.valueOf(annotationArr[3]);
							int exonEndTem = Integer.valueOf(annotationArr[4]);
							Integer[] temArr = { exonStartTem, exonEndTem };
							exonSiteList.add(temArr);
							if (geneStart > exonStartTem) {
								geneStart = exonStartTem;
							}
							if (geneEnd < exonEndTem) {
								geneEnd = exonEndTem;
							}
						}

					}
				}
				annotation = annotationRead.readLine();
			}
			// 最后一个
			String exonListMapKey = chr + "\t" + geneId;
			exonListMap.put(exonListMapKey, exonSiteList);
			String[] geneTem = { chr, geneId, geneStart + "", geneEnd + "" };
			if (!geneExonMap.containsKey(chr)) {
				geneSiteList = new ArrayList<SiteSort>();
				geneSiteList.add(new SiteSort(geneStart, geneTem));
				geneExonMap.put(chr, geneSiteList);
			} else {
				geneSiteList = geneExonMap.get(chr);
				geneSiteList.add(new SiteSort(geneStart, geneTem));
				geneExonMap.put(chr, geneSiteList);
			}
		}

		//////// 排序
		for (String chrKey : geneExonMap.keySet()) {
			geneSiteList = geneExonMap.get(chrKey);
			Collections.sort(geneSiteList);
			geneExonMap.put(chrKey, geneSiteList);
		}
	}

	public HashMap<String, String> getChrExonStartMap() {
		return chrExonStartMap;
	}

	public HashMap<String, String> getChrExonEndMap() {
		return chrExonEndMap;
	}

	public HashMap<String, ArrayList<SiteSort>> getGeneExonMap() {
		return geneExonMap;
	}

	public HashMap<String, ArrayList<Integer[]>> getExonListMap() {
		return exonListMap;
	}

	public HashMap<String, ArrayList<String>> getChrExonStartTranscriptMap() {
		return chrExonStartTranscriptMap;
	}

	public HashMap<String, ArrayList<String>> getChrExonEndTranscriptMap() {
		return chrExonEndTranscriptMap;
	}
}
