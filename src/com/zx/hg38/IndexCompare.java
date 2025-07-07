package com.zx.hg38;

import java.util.HashMap;

public class IndexCompare {
	private String[][] bibasesMut = { { "AC", "AG", "GC", "AG", "AT", "AC", "AT", "AG"}, 
			                          { "CT", "GT", "CT", "GC", "GT", "AT", "CT", "AT"} };
	private String[][] bibases = { { "AC", "AG" }, { "CT", "GT" } };
	private String[] strandIndex = { "-", "+" };
	HashMap<Integer, String> indexStrandMap = new HashMap<Integer, String>();
	public HashMap<Integer, String> indexCompare(String end_string1, String end_string2) {
		indexStrandMap.clear();
		String[] upEndString = { end_string1.toUpperCase(), end_string2.toUpperCase() };
		for (int i = 0; i <= 1; i++) {
			int preIndex = -1;
			while (true) {
				int index = upEndString[0].indexOf(bibases[0][i], preIndex + 1);
				if (index == -1) {
					break;
				}
				if (upEndString[1].length() >= index + 2) {
					if (upEndString[1].substring(index, index + 2).equalsIgnoreCase(bibases[1][i])) {
						indexStrandMap.put(index, index + "\t" + strandIndex[i] + "\t" + bibases[0][i] + "\t" + bibases[1][i]);
					}
				} else {
					break;
				}
				preIndex = index;
			}

		}
		return indexStrandMap;
	}
	public HashMap<Integer, String> indexCompareChrM(String end_string1, String end_string2) {
		indexStrandMap.clear();
		String[] upEndString = { end_string1.toUpperCase(), end_string2.toUpperCase() };
		for (int i = 0; i <= 7; i++) {
			int preIndex = -1;
			int stand = i%2;
			while (true) {
				int index = upEndString[0].indexOf(bibasesMut[0][i], preIndex + 1);
				if (index == -1) {
					break;
				}
				if (upEndString[1].length() >= index + 2) {
					if (upEndString[1].substring(index, index + 2).equalsIgnoreCase(bibasesMut[1][i])) {
						indexStrandMap.put(index, index + "\t" + strandIndex[stand] + "\t" + bibasesMut[0][i] + "\t" + bibasesMut[1][i]);
					}
				} else {
					break;
				}
				preIndex = index;
			}

		}
		return indexStrandMap;
	}
}
