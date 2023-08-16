package com.zx.hg19;

public class IsInCircRNA2 {
	int miss_count_min = 3,miss_count_max = 5;
	int window_step = 5,window_size = 10;	
	public int isInCircRNA2(String unMapSeq,String circ_range_seq) {
		int SeqLen = unMapSeq.length();
		int[] miss_count3 = {0,0,0};												
		for (int j = 0; j <= (int)(SeqLen-window_size-0.0)/window_step; j++) {
			String seq;
			if(SeqLen<=10) {
				 seq = unMapSeq.substring(0,SeqLen);
			}else {
				seq = unMapSeq.substring(j*window_step,j*window_step+window_size);
			}
		   
			if (circ_range_seq.indexOf(seq) >= 0) {
				miss_count3[0] = 0;
				
			}else {
				miss_count3[1]++;
				miss_count3[0]++;
				if (miss_count3[0] > miss_count3[2]) {
					miss_count3[2] = miss_count3[0];
				}
			}	
			}												
		if (miss_count3[2] > miss_count_max) {
			return 0;
			
		}else if ((miss_count3[1]-1)*2 > (int)(SeqLen-window_size)/window_step ) {
			return 0;
			
		}else if (miss_count3[0] <= miss_count_max && miss_count3[0] >= miss_count_min) {
			
		}
		return 1;
	}
}
