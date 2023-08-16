package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FileFindCircRNAScan2 {
	int readNum = 0;
	IsBSJScan2 isBSJScan2;
	IsStand stand = new IsStand();
	HashMap<String, Integer> circFSJMap;

	public FileFindCircRNAScan2(int minMapqUni, HashMap<String, Integer> circFSJMap, int linear_range_size_min,
			byte[] siteArray1, byte[] siteArray2, HashMap<Integer, ArrayList<SiteSort>> chrSiteMap1,
			HashMap<Integer, ArrayList<SiteSort>> chrSiteMap2, HashMap<String, String> chrTCGAMap, int seqLen,
			boolean intronLable) throws IOException {

		super();
		if (intronLable) {
			isBSJScan2 = new IsBSJScan2Intron(minMapqUni, circFSJMap, linear_range_size_min, siteArray1, siteArray2,
					chrSiteMap1, chrSiteMap2, chrTCGAMap, seqLen);
		} else {
			isBSJScan2 = new IsBSJScan2(minMapqUni, circFSJMap, linear_range_size_min, siteArray1, siteArray2,
					chrSiteMap1, chrSiteMap2, chrTCGAMap, seqLen);
		}
	}

	public ArrayList<String> getCircList() {
		return isBSJScan2.getCircList();
	}

	public HashMap<String, Integer> getCircFSJMap() throws IOException {
		return isBSJScan2.getCircFSJMap();
	}

	public void setCircScan2List() {
		isBSJScan2.setCircScan2List();
	}

	public void setFSJScan2List() {
		isBSJScan2.setFSJScan2List();
	}

	public long getReadNum() {
		return readNum;
	}

	public void findCircRNAScan2(String input, HashMap<String, String> scan1IdMap) throws IOException {
		readNum = 0;
		boolean lable = false;
		// 存放第一遍扫描的circRNA id

		BufferedReader samReadScan2 = new BufferedReader(new FileReader(new File(input)));
		//////////////////
		String id = "", oneRead = "";
		HashMap<Integer, ArrayList<String>> standMap = new HashMap<Integer, ArrayList<String>>();
		ArrayList<String> serveInforList = new ArrayList<String>();
		standMap.put(0, serveInforList);
		standMap.put(1, serveInforList);
		HashMap<Integer, String> readsMap = new HashMap<Integer, String>();
		String line = samReadScan2.readLine();
		while (line != null) {
			if (line.startsWith("@")) {
				line = samReadScan2.readLine();
				continue;
			}
			String lineArr[] = line.split("\t");
			int standKey = Integer.valueOf(stand.stand7(lineArr[1]));
			String serveInfor = lineArr[0] + "\t" + stand.stand5(lineArr[1]) + "\t" + lineArr[2] + "\t" + lineArr[3]
					+ "\t" + lineArr[4] + "\t" + lineArr[5];
			if (!id.equals(lineArr[0])) {
				// 判断是否匹配上
				if (lable) {
					readNum++;
				}
				if (scan1IdMap.containsKey(id)) {

				} else {
					isBSJScan2.isCandidate(standMap, readsMap);

				}
				standMap = new HashMap<Integer, ArrayList<String>>();
				serveInforList = new ArrayList<String>();
				serveInforList.add(serveInfor);
				standMap.put(standKey, serveInforList);
				readsMap.clear();
				readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
				id = lineArr[0];
				oneRead = standKey + "";
			} else {
				if (!stand.stand7(lineArr[1]).equalsIgnoreCase(oneRead)) {
					readsMap.put(standKey, stand.stand5(lineArr[1]) + lineArr[9].toUpperCase());
					serveInforList = new ArrayList<String>();
					serveInforList.add(serveInfor);
					standMap.put(standKey, serveInforList);
					oneRead = standKey + "";

				} else {
					serveInforList = standMap.get(standKey);
					serveInforList.add(serveInfor);
					standMap.put(standKey, serveInforList);
				}
			}

			if (!lineArr[5].equals("*")) {
				lable = true;
			}
			line = samReadScan2.readLine();

		}
		samReadScan2.close();

		// 判断是否匹配上
		if (lable) {
			readNum++;
		}
		// 判断是否含有BSJ
		if (scan1IdMap.containsKey(id)) {

		} else {
			isBSJScan2.isCandidate(standMap, readsMap);

		}

	}

}