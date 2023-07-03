package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class MutSamToBam {
	ArrayList<String> filePathList = new ArrayList<String>();
	ArrayList<String> fileNameList = new ArrayList<String>();
	public void getBam(String samFolder,int numClass) throws IOException {
		ArrayList<Long> fileSizeList = new ArrayList<Long>();
		
		String fileClass = "sam";
		BufferedReader samFileRead = new BufferedReader(new FileReader(new File(samFolder)));			
		String fileLine = samFileRead.readLine();
		while (fileLine != null) {
			if (fileLine.startsWith("#") || fileLine.equals("")) {
				fileLine = samFileRead.readLine();
				continue;
			}
			File file = new File(fileLine);
			fileSizeList.add(file.length());
			fileClass = fileLine.substring(fileLine.length()-3,fileLine.length());
			String[] arrTem = fileLine.split("/");
			fileNameList.add(arrTem[arrTem.length-1]);			
			filePathList.add(fileLine);
			fileLine = samFileRead.readLine();		
		}
		samFileRead.close();
	}
	
}
