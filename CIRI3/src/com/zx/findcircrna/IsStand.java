package com.zx.findcircrna;

public class IsStand {
	public String stand7(String infor) {
		int stand = Integer.parseInt(infor);
		// 转换成二进制
		String standNumber = Integer.toBinaryString(stand);
		if (standNumber.length()<7) {
			return "0";
		}else {
			return standNumber.substring(standNumber.length() - 7, standNumber.length() - 6);
		}
		
	}
	public String stand5(String infor) {
		int stand = Integer.parseInt(infor);
		// 转换成二进制
		String standNumber = Integer.toBinaryString(stand);
		if (standNumber.length()<5) {
			return "0";
		}else {
			return standNumber.substring(standNumber.length() - 5, standNumber.length() - 4);	
		}		
	}
}
