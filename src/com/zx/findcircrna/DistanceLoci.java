package com.zx.findcircrna;

import java.util.ArrayList;

public class DistanceLoci {
public Integer distanceLoci(ArrayList<Integer> locusList, ArrayList<Integer> locus2List,int window_step) {
	if (locusList.size() < 2 && locus2List.size() <2) {
		return 0;
	}
	int  locusSum = 0;
	int  locus2Sum = 0;
	for (int i = 1; i < locusList.size(); i++) {
		locusSum += Math.abs(locusList.get(i)-locusList.get(i-1));
	}
	for (int i = 1; i < locus2List.size(); i++) {
		locus2Sum += Math.abs(locus2List.get(i)-locus2List.get(i-1));
	}
	if (locusSum <= window_step * locusList.size() && locusSum *20 < locus2Sum) {
		return 1;
	}else {
		return 0;
	}	
}
}
