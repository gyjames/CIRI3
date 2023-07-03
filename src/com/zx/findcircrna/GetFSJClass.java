package com.zx.findcircrna;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GetFSJClass {
	HashSet<String> FSJSet = new HashSet<String>();
	int seqLen;
	public GetFSJClass(int seqLen) {
		this.seqLen = seqLen;
	}
	public HashSet<String> getFSJ(int startTem,int endTem,String chr,int style,byte[] siteArray1,byte[] siteArray2,
			HashMap<Integer, ArrayList<SiteSort>> siteMap1,HashMap<Integer, ArrayList<SiteSort>> siteMap2) {
		FSJSet.clear();
		int num1 = startTem/seqLen;
		int num2 = endTem/seqLen;
		if (siteArray1[num1] == 1) {
			ArrayList<SiteSort> siteList1 = siteMap1.get(num1);
			for (int k = siteList1.size()-1; k >= 0; k--) {
				if (siteList1.get(k).getSite() >= startTem) {
					String[] temChrSite = siteList1.get(k).getLength();
					if (siteList1.get(k).getSite() <= endTem || style == 10 || style == 1) {
						String temKey = chr + "\t" +temChrSite[0] + "\t" + temChrSite[1];
						FSJSet.add(temKey);
					}
				}else {
					break;
				}
			}
		}
		if (siteArray2[num1] == 1) {
			ArrayList<SiteSort> siteList2 = siteMap2.get(num1);
			for (int k = siteList2.size()-1; k >= 0; k--) {
				if (siteList2.get(k).getSite() >= startTem) {
					String[] temChrSite = siteList2.get(k).getLength();
					if (siteList2.get(k).getSite() <= endTem || style == 10 || style == 1) {
						String temKey = chr + "\t" +temChrSite[0] + "\t" + temChrSite[1];
						FSJSet.add(temKey);
					}
				}else {
					break;
				}
			}
		}
		if (num2 != num1) {
			if(siteArray1[num2] == 1) {
				ArrayList<SiteSort> siteList1 = siteMap1.get(num2);
				for (int k = 0; k < siteList1.size(); k++) {
					if (siteList1.get(k).getSite() <= endTem) {
						String[] temChrSite = siteList1.get(k).getLength();
						if (siteList1.get(k).getSite() >= startTem || style == 10 || style == (-1)) {
							String temKey = chr + "\t" +temChrSite[0] + "\t" + temChrSite[1];
							FSJSet.add(temKey);
						}
					}else {
						break;
					}
				}
			}
			if(siteArray2[num2] == 1) {
				ArrayList<SiteSort> siteList2 = siteMap2.get(num2);
				for (int k = 0; k < siteList2.size(); k++) {
					if (siteList2.get(k).getSite() <= endTem) {
						String[] temChrSite = siteList2.get(k).getLength();
						if (siteList2.get(k).getSite() >= startTem || style == 10 || style == (-1)) {
							String temKey = chr + "\t" +temChrSite[0] + "\t" + temChrSite[1];
							FSJSet.add(temKey);
						}
					}else {
						break;
					}
				}
			}
		}		
		return FSJSet;
	}
}
