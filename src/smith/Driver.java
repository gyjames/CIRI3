package smith;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;



/**
 *
 * @author DulanDias (AS2014347)
 */
public class Driver {

    /**
     * @param args the command line arguments
     * @throws IOException 
     */
    public static void main(String[] args) throws IOException {
    	BufferedReader hg = new BufferedReader(new FileReader(new File("D:/BioInformation/hg19/hg19.fa")));
		ArrayList<String> hgList = new ArrayList<String>();	
		HashMap<String, String> chrTCGAMap = new HashMap<String, String>();
		String chrName = "",chrTAGA;	
		String HgLine = hg.readLine();
		while (HgLine != null) {
			if (HgLine.startsWith(">chr")) {		
			    chrTAGA = String.join("", hgList).toUpperCase();
			    chrTCGAMap.put(chrName, chrTAGA);			    
			    chrName = HgLine.replace(">", "");
			    hgList.clear();
			    HgLine = hg.readLine();
				continue;
			}
			hgList.add(HgLine);
			HgLine = hg.readLine();
		}
		hg.close();
		chrTAGA = String.join("", hgList).toUpperCase();
		chrTCGAMap.put(chrName, chrTAGA);
		chrTAGA=null;
		System.out.println("准备文件导入完毕");  
		
		
		int match = 2;
        int mismatch = -4;
        int gap = -4;
        SmithWaterman aligner = new SmithWaterman(match, mismatch, gap);
		BufferedReader br = new BufferedReader(new FileReader(new File("D:/Ciri/MapAgian/SRR650317_scan2.txt")));	
		
		String line = br.readLine();
		while (line != null) {
			String[] arr = line.split("\t");
			String chr = arr[2];
			chrTAGA = chrTCGAMap.get(chr);
			int uncertainSeqLen = arr[6].length();
			int strat = Integer.valueOf(arr[4]);
			int end = Integer.valueOf(arr[5]);
			if(arr[3].equals("sm")) {
				String linear_range;
				if (end - strat +5 >= 50000) {
					if (2*strat >= end +6) {	
						linear_range = chrTAGA.substring(2*strat-end -6, strat-1);
					}else {
						linear_range = chrTAGA.substring(0, strat-1);
					}
				}else {
                    if(strat >= 50000) {
                    	linear_range = chrTAGA.substring(strat-50000-1, strat-1);
                    }else {
                    	linear_range = chrTAGA.substring(0, strat-1);
					}										
				}
				
				if(linear_range.indexOf(arr[6])>=0) {
					line = br.readLine();
					continue;
				}
				
				String mapRef = chrTAGA.substring(end-uncertainSeqLen-3,end);
				aligner.setSeq(arr[6],mapRef);
				String[] alignment= aligner.getAlignment();
				//alignment[0]表示read alignment[1]表示ref
				if (alignment[0].length()>=(uncertainSeqLen-2)) {
					mapRef = chrTAGA.substring(strat-1,end);
					if (!arr[8].equals("*")) {
						aligner.setSeq(arr[8],mapRef);
						alignment= aligner.getAlignment();
						if (alignment[0].length()>=(arr[8].length()-2)) {
							
							System.out.println(arr[2]+"\t"+arr[4]+"\t"+arr[5]+"\t"+"1"+"\t"+arr[13]+"\t"+arr[10]+"\t"+arr[0]);
							/*aligner.setSeq(arr[7],mapRef);
							alignment= aligner.getAlignment();*/
						}
					}else {
						System.out.println(arr[2]+"\t"+arr[4]+"\t"+arr[5]+"\t"+"1"+"\t"+arr[13]+"\t"+arr[10]+"\t"+arr[0]);
					}
				}
			}else {
				String linear_range;
				if (end - strat +5 >= 50000) {	
					if(2*end-strat+5>chrTAGA.length()) {
                    	linear_range = chrTAGA.substring(end,chrTAGA.length());
                    }else {
                    	linear_range = chrTAGA.substring(end,2*end-strat+5);
					}
					
				}else {
                    if(end+50000>chrTAGA.length()) {
                    	linear_range = chrTAGA.substring(end,chrTAGA.length());
                    }else {
                    	linear_range = chrTAGA.substring(end,end+50000);
					}									
				}
	
				if(linear_range.indexOf(arr[6])>=0) {
					line = br.readLine();
					continue;
				}
				String mapRef = chrTAGA.substring(strat-1,strat+uncertainSeqLen+2);
				aligner.setSeq(arr[6],mapRef);
				String[] alignment= aligner.getAlignment();
				//alignment[0]表示read alignment[1]表示ref
				if (alignment[0].length()>=(uncertainSeqLen-2)) {
					mapRef = chrTAGA.substring(strat-1,end);
					if (!arr[8].equals("*")) {
						aligner.setSeq(arr[8],mapRef);
						alignment= aligner.getAlignment();
						if (alignment[0].length()>=(arr[8].length()-2)) {
							System.out.println(arr[2]+"\t"+arr[4]+"\t"+arr[5]+"\t"+"1"+"\t"+arr[13]+"\t"+arr[10]+"\t"+arr[0]);
							/*aligner.setSeq(arr[7],mapRef);
							alignment= aligner.getAlignment();*/
						}
					}else {
						System.out.println(arr[2]+"\t"+arr[4]+"\t"+arr[5]+"\t"+"1"+"\t"+arr[13]+"\t"+arr[10]+"\t"+arr[0]);
					}
				}
			}
			
			line = br.readLine();
		}
		
		br.close();
		
		
		
		/*chrTAGA = String.join("", hgList).toUpperCase();
		String end = chrTAGA.substring(34924307,34924319);
    	
        String start = "GCAGTGTATACTGGTTATGAATGCGTTAGATAACAGAG";
        //String end = "ATTGATCTGGATACTATTGATGTAAGCAACCTCAACAGACAGTTTTTGTTTCAAAAGAAACATGTTGGAAGATCAAAGGCACAGGTAACTATATTTCTCATACCATTTCTATAACTTGATGGAGCTTCTATTTGTGATACCAGATTAATCCTGCCTGTCATATAGGAGGGAAAAAATGGGTGAAGCTCTCATTTCAGCTGAGAGATTGTGAAGGATAAAGATATATCAAGTAAAGGACTAAAGATTAAGCCCTCCTTTTTTGTTCTAAAAAGGAAATAATCCCGTGATGGTGTAACTTTCATGGCCTCTTTATATTTAAGTTTCTGATGCTTTGCTTAGGATGCACCTTGTGTTCCTACACCCCCTAGAGGTTGACTACTGCTATGTTTGGGGGCTAGAGGAGCGTTATGCTGCTGTCCTAGTTGAAACAAGACAGTTCCCCCGATCTTTGTTTAGTTTTTTTAGTGTTCCTTGGTGATTCAGAGTCTTAAAACATCTTTTCCCTTTGCCTGTCTCCCAACCTGTTTGGCAAACCCCCACTTGTTTCTCTTTTTCTTTCATATTTTGGCTTCCTCCTCCAGGAGTACTTTCCTAAACCTTTTAGACTAAGTCAGGTTGGTTTTTCTCTCATAACAACCTTTACTTTACTTTTCCTATCATTTTCTTTTCTTTTTTTTTTTTTTTTTTTTTTTTAGAGATAGGGTCTCACTCTGTCGCCCCAGGCTGGAGTGGAGTGATGCGGTCATGGCTTACTGCAGCCTAGAACTCCTGGCCTCAAGGGATCCTCCTACACTGGCCTCCCAAGTGCTGGGATTACAGGCATGAGCCACCATGCCTGGCCTTTTCCTGTCATTCTCATGGCATTTATCAGAGTCTGTAGTTGAACTTAATCTGGATTGTTATTTGATTAATGTCTCTCTTGCTTGACATGTACGTTGTATAGTTCCTGGTACATGCAAGATGCTCAGTAAATACTTGTTAAGTGAATGAATGAATGACACTTGGTTACTTCATGTTGTATTTTTGCTGAAGCTTGCCATTCATGTAGAAGTGTGTTTGGTGGTACGTGGGTCAGGGGGGGCTCTCCGGTTATTTTATAAGTATAAAAGTATAATTTTATAAGTGTAGCTTAACCCTAATGCAATTTTCTCACATTTTTGGACTTCTTTGGTGTTTAATATAAATAAGGTACTATGTAAGTAGATAACAGCATTGTGGCTTCCAGAAGCAGATTTCAGATAGGAACAGAAATATTGTAGTAATTCAGTGTTGTTTATTTTTCCAGGTTGCCAAGGAAAGTGTACTGCAGTTTTACCCGAAAGCTAATATCGTTGCCTACCATGACAGCATCATGAAGTATGCTATAGTGATTACATTGCAAAGTTGTATAAGGGTTTTGTAAGCCAAATATATAAGCTCAAAGTCATTCAGCTTTTTTAAAAAAAATGATTTTTCTAGAATTATTAAACAGTGGTGGTTTCTGCTTACAGAGTGGCATTCTTTTTACCAGTTCATAAATACATGGTTTGAATGATTTAGAATGGTGAGATAGTTTAATACTTTTAATAAATAAGTTACTTGGAAATTTATAAAGACCTTATATGCTTGCCTTAAGTAAGAAGATATTTTAAACTATGAAATATCCTTATTTATGATGATTCAGACTTTTCTAAAGTGTTTGCCTTTTTTTGTTTAAAAGTTTCCCTCCCAAAATAATCTGTAGTATGACTCTGACTTCTGTAATATAACACTTGGAATTGTGGTAGGGCAGTAGTATTGACAGTACATAAATATATTTATCTCAGTAAAGATAGAGTAAATATATAATTTTGGAGGAGTTGGAGCAAAGGAACTGTGGGCCGAAATTGTCCTTTATCTATGGTGGTTTAGCTTGCAAGTAAATAGTGTTAAAAATAGTCCTTGTTTGAGGTACTGGAATTGAGGTACTTAACACTTTCTCAGGGCATTTGTATCTCCGGATTGTGTTTTTTTTGTGTGTGAGATAGGGTTTTGCGCTTATTGCCCAGGCTGGAGTACAGTGGCATGATCTTGGCTCACTGCAACCTCCGCCTCCTGGGTTCAAGCGATTCTCCTGTCTCAGCCTCCTGAGTAGTTGGGATTACAGGCGCACGCCACCACACTCAGCTAATTTTTGTATTTTAGTAGGCGGGGTTTCATCATGTTGGCCAGGCTGGTCTCAAACTCCCGACCTCAGGTGATCCACCCGCCTTGGCCTCCCAAAGTGCAGGGATTACAGGCGTGAGCCACTGCACCCGGCCTGTATCTCCGGCTTGTTTTTTAGTTTTTTACTCTGTGAACCTGCCTGTAGCCCTATAACCTCCTTACCCTGTGTTGCTTTCCCAGTGCTGCTAAAAACAGAGACCCTAGGGGCCTGAACCCTGCCTGATCCTGATCATGTGAGGTGAAGTCAGAGGCCAGGCAGATCCTGGGAATGTGTTCTTTATGCGTTCCCAGGCACTTGCTGCCCAAAGAAGGCCTAGTGGGACTAGTGAAGTCATGCCATCTTTCTGCCAGTTTGAGATCTTTCTCTTGTGCCCTAGCTGGCCAAAAAGGATAGGAATTTGAAGGAAATAATTTCTTTAAAATCTGCTTATGTTTGGAAAGAGTCAAACCTTACTATCTATATGACATATAATCCTGATATCAGTTTTGTTTACAGTTACTATGTTTTTTTCAGAACTGTTTATTATACTGCAAAGATCATGTGGAAGGCTAGTTATTTTGGTGACCTTTTTTTATTTTGTTTTGTAGCCCTGACTATAATGTGGAATTTTTCCGACAGTTTATACTGGTTATGAATGCTTTAGATAACAGAG";
		
		aligner.setSeq(start,end);
        String[] alignment= aligner.getAlignment();
        System.out.println(alignment[0]);
        System.out.println(alignment[1]);*/
       
        
    }
    
}
