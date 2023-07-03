package com.zx.hg19;

import java.util.ArrayList;

import com.zx.findcircrna.DistanceLoci;

public class IsInCircRNA3 {
	int miss_count_min = 3,miss_count_max = 5;
	int window_step = 5,window_size = 10;
	ArrayList<Integer> locusListTem = new  ArrayList<Integer> ();
	ArrayList<Integer> locus2ListTem = new  ArrayList<Integer> ();
	DistanceLoci distanceLoci = new DistanceLoci();
	public int isInCircRNA3(String anoRead,String preJudge,String circ_range_seq,String pem_null_range_seq) {
		
		int anoReadLen = anoRead.length();
		int trial = (int) ((anoReadLen-window_size+0.0)/window_step);
		//$miss_count, $total_miss_count, $cont_miss_count
		int[] miss_count = {0,0,0};
		int[] miss_count2 = {0,0,0};
		locusListTem.clear();
	    locus2ListTem.clear();
	    for (int j = 0; j <= trial; j++) {
			String seq = anoRead.substring(j*window_step,window_size+j*window_step);
			int locus = circ_range_seq.indexOf(seq);														
			if (locus >= 0) {
				locusListTem.add(locus);
				miss_count[0] = 0;
			}else {
				miss_count[1]++;
				miss_count[0]++;
				if(miss_count[0] > miss_count[2]) {
					 miss_count[2] = miss_count[0];
				}
			}
			if (pem_null_range_seq.length()>0) {
				int locus2 = pem_null_range_seq.indexOf(seq);
				if (locus2>=0) {
					locus2ListTem.add(locus2);
				}else {
					miss_count2[1]++;
				}	
			}																	
		}
	   
		//0:id 1:stand 2:chr 3:"sm" 4:start 5:end 6:str 7:str2 8:str3 9:str2_ok 10:+	 11:AG	12:GT 13:cir
		if (pem_null_range_seq.length()>0) {
			if (miss_count2[1] == 0 && miss_count[1] == 0) {
				if (distanceLoci.distanceLoci(locusListTem, locus2ListTem, window_step) == 1) {
					
					return 1;						
				}else {											
					return -2;	
				}
			}else if (miss_count2[1] <= miss_count[1]) {
				if (locus2ListTem.size()>0) {
					return -2;
				}else {
					
					return -1;
				}
			}else if (miss_count[1]*4 > trial*3 && preJudge.equals("0") ) {
				return -2;
			}else if (miss_count[2] > miss_count_max ) {
				return -1;
			}else if (miss_count[1]*2 > trial ) {
				return -1;
			}else {
				
				return 1;
			}
		}else if (miss_count[1]*4 > trial*3 && preJudge.equals("0") ) {
			return -2;
		}else if (miss_count[2] > miss_count_max ) {
			return -1;
		}else if (miss_count[1]*2 > trial ) {
			return -1;
		}
		
		return 1;
		
	}
}
