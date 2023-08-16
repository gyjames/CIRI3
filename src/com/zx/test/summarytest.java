package com.zx.test;

import java.util.ArrayList;

import com.zx.findcircrna.Summary;

public class summarytest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Summary summary = new Summary(2);
		ArrayList<String> SummaryCircList = summary.summary(samFile,1,circFSJNewMap,UserGivecircRNA);
	}

}
