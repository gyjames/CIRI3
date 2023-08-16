package com.zx.hg38;

import java.util.ArrayList;

import com.zx.findcircrna.DistanceLoci;

public class IsInCircRNA1_2 {
	int miss_count_min = 3,miss_count_max = 5;
	int[] window_unit = {9, 7, 5, 4, 3};
	ArrayList<Integer> locusListTem = new  ArrayList<Integer> ();
	ArrayList<Integer> locus2ListTem = new  ArrayList<Integer> ();
	DistanceLoci distanceLoci = new DistanceLoci();
	public int isInCircRNA1_2(int len_str,String str,String circ_range_seq,String linear_range) {
		for (int k = 0; k < window_unit.length; k++) {
			if (len_str < window_unit[k] * 2) {
				continue;
			}
			locusListTem.clear();
		    locus2ListTem.clear();
			int window_step = window_unit[k];
			int window_size = window_unit[k]*2;
			int trial = (int) ((len_str-window_size+0.0)/window_step);
			
			//$miss_count, $total_miss_count, $cont_miss_count
			int[] miss_count = {0,0,0};
			int[] miss_count2 = {0,0,0};													
			for (int j = 0; j <= trial; j++) {
				
				String seq = str.substring(j*window_step,j*window_step+window_size);
				
				int locus = circ_range_seq.indexOf(seq);
				int locus2 = linear_range.indexOf(seq);	
				
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
				if (locus2>=0) {
					locus2ListTem.add(locus2);
				}else {
					miss_count2[1]++;
				}				
			}
			if (len_str % window_unit[k] != 0) {
				trial++;
				String seq = str.substring(len_str-window_size,len_str);
				int locus = circ_range_seq.indexOf(seq);
				int locus2 = linear_range.indexOf(seq);
				
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
				if (locus2>=0) {
					locus2ListTem.add(locus2);
				}else {
					miss_count2[1]++;
				}								
			}	
			if (miss_count2[1] ==0 && miss_count[1] == 0) {
				if (distanceLoci.distanceLoci(locusListTem, locus2ListTem, window_step) == 1) {
					return 1; 	
				}else {
					
					return 0;
				}				
			}else if (miss_count2[1] <= miss_count[1]) {
				if (locus2ListTem.size()>0) {
					
					return 0;
				}else {
					
				}
			}else if (miss_count[2] > miss_count_max ) {
				
			}else if (miss_count[1]*2 > trial ) {
				
			}else {
				return 1; 	
			}
			
		} 
		 return 0;		
	}
}
