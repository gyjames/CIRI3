package com.zx.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import DE.BSJ_no;
import DE.BSJ_yes;
import DE.Ratio;
import DE.Relative;

public class TestParameters {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String[] Parameters = args;
		ArrayList<String> dieReasonList = new ArrayList<String>();	
		if (Parameters[0].equals("DE_BSJ")) {
			String inforPath = null,outPath = null,matrixPath = null,genePath = null;
			double pval;
			HashMap<String,String> parameterMap = new HashMap<String,String>();
			parameterMap.put("-I", "F");
			parameterMap.put("-O", "F");
			parameterMap.put("-M", "F");
			parameterMap.put("-G", "F");
			parameterMap.put("-P", "0.05");
			parameterMap.put("-H", "F");
			parameterMap.put("--in", "F");
			parameterMap.put("--out", "F");
			parameterMap.put("--matrix", "F");
			parameterMap.put("--gene", "F");
			parameterMap.put("--pvalue", "0.05");
			parameterMap.put("--help", "F");
			
			
			for (int i = 1; i < Parameters.length-1; i++) {
				if (Parameters[i].startsWith("-")) {
					parameterMap.put(Parameters[i], Parameters[i+1]);
				}
			}
			if (Parameters[Parameters.length-1].startsWith("-")) {
				parameterMap.put(Parameters[Parameters.length-1], "T");
			}
						
			if (!parameterMap.get("-H").equals("F")) {
				System.out.println(
						"without replicate"
						+ "Usage:    java CIRI3.jar DE_BSJ -I infor.tsv -O output.txt \r\n"
						+ "with replicate"
						+ "Usage:    java CIRI3.jar DE_BSJ -I infor.tsv -G gene_expression.txt -O output.txt \r\n"
						+ "\r\n"
						+ "Arguments:\r\n"
						+ "\r\n"
						+ "    -I, --in\r\n"
						+ "          the file of sample list(required)\r\n"
						+ "    -O, --out\r\n"
						+ "          output differential expression result(required)\r\n"
						+ "    -M, --matrix\r\n"
						+ "          circRNA BSJ expression matrix (optional)\r\n"
						+ "    -G, --gene\r\n"
						+ "          gene expression matrix \r\n"
						+ "    -P, --pvalue\r\n"
						+ "          p value threshold for DE score calculation (default: 0.05)\r\n"
						+ "    -H, --help\r\n"
						+ "          show this help information\r\n");
			}else if (parameterMap.get("-I").equals("F") && parameterMap.get("--in").equals("F")) {
				System.out.println("Please use --in or -I option to designate input the file of sample list!");
			}else if (parameterMap.get("-O").equals("F") && parameterMap.get("--out").equals("F")) {
				System.out.println("Please use --out or -O option to designate output file!");
			}else {
				//input
				if (!parameterMap.get("-I").equals("F") || !parameterMap.get("--in").equals("F")) {
					if(!parameterMap.get("-I").equals("F")) {
						inforPath = parameterMap.get("-I");
					}else {
						inforPath = parameterMap.get("--in");
					}
					File inputFileTem = new File(inforPath); 
					inforPath = inputFileTem.getCanonicalPath();
				}
				File infile = new File(inforPath);
				if (!infile.isFile()) {
					dieReasonList.add("No sample list file found at designated directory!");
				}
				//output
				if (!parameterMap.get("-O").equals("F") || !parameterMap.get("--out").equals("F")) {
					if(!parameterMap.get("-O").equals("F")) {
						outPath = parameterMap.get("-O");
					}else {
						outPath = parameterMap.get("--out");
					}	
					File outputFileTem = new File(outPath); 
					outPath = outputFileTem.getCanonicalPath();
				}
				File outfile = new File(outPath);
				if (!outfile.isFile()) {
					dieReasonList.add("Output file cannot be written in the directory "+outPath+"!");
				}
				//matrix
				if (!parameterMap.get("-M").equals("F") || !parameterMap.get("--matrix").equals("F")) {
					if(!parameterMap.get("-M").equals("F")) {
						matrixPath = parameterMap.get("-M");
					}else {
						matrixPath = parameterMap.get("--matrix");
					}	
					File matrixPathTem = new File(matrixPath); 
					matrixPath = matrixPathTem.getCanonicalPath();
					File matrixfile = new File(matrixPath);
					if (!matrixfile.isFile()) {
						dieReasonList.add("No circRNA BSJ expression matrix file found at designated directory!");
					}
				}		
				//gene
				if (!parameterMap.get("-G").equals("F") || !parameterMap.get("--gene").equals("F")) {
					if(!parameterMap.get("-G").equals("F")) {
						genePath = parameterMap.get("-G");
					}else {
						genePath = parameterMap.get("--gene");
					}	
					File genePathTem = new File(genePath); 
					genePath = genePathTem.getCanonicalPath();
					File genefile = new File(genePath);
					if (!genefile.isFile()) {
						dieReasonList.add("No gene expression matrix file found at designated directory!");
					}
				}		
				//Pvalue
				if(!parameterMap.get("-P").equals("0.05")) {
					pval = Integer.valueOf(parameterMap.get("-P"));		
				}else if(!parameterMap.get("--pvalue").equals("0.05")) {
					pval = Integer.valueOf(parameterMap.get("--pvalue"));
				}else {
					pval = 0.05;
				}
				if (pval <= 0 || pval >= 1) {
					dieReasonList.add("The p value threshold can only be between 0 and 1!");
				}
				//运行
				if(genePath == null) {
					BSJ_no bn = new BSJ_no();
					if (matrixPath == null) {
						bn.DEBSJNoFile(inforPath, outPath, pval);				
					}else {
						bn.DEBSJNoMatrix(inforPath, matrixPath, outPath, pval);
					}	
				}else {
					BSJ_yes by = new BSJ_yes();					
					if (matrixPath == null) {
						by.DEBSJNoFile(inforPath, genePath, outPath);
					}else {						
						by.DEBSJNoMatrix(inforPath, genePath, outPath, matrixPath);
					}	
				}
				    	
	        }
		}else if (Parameters[0].equals("DE_Ratio")) {
			String inforPath = null,outPath = null,BSJMatrixPath = null,FSJMatrixPath = null;
			int thread;
			HashMap<String,String> parameterMap = new HashMap<String,String>();
			parameterMap.put("-I", "F");
			parameterMap.put("-O", "F");
			parameterMap.put("-BM", "F");
			parameterMap.put("-FM", "F");
			parameterMap.put("-T", "1");
			parameterMap.put("-H", "F");
			parameterMap.put("--in", "F");
			parameterMap.put("--out", "F");
			parameterMap.put("--bsj_matrix", "F");
			parameterMap.put("--fsj_matrix", "F");
			parameterMap.put("--thread_num", "1");
			parameterMap.put("--help", "F");
			
			
			for (int i = 1; i < Parameters.length-1; i++) {
				if (Parameters[i].startsWith("-")) {
					parameterMap.put(Parameters[i], Parameters[i+1]);
				}
			}
			if (Parameters[Parameters.length-1].startsWith("-")) {
				parameterMap.put(Parameters[Parameters.length-1], "T");
			}
						
			if (!parameterMap.get("-H").equals("F")) {
				System.out.println(
						"Usage:    java CIRI3.jar DE_Ratio -I infor.tsv -BM BSJMatrix.txt -FM FSJMatrix.txt -O output.txt \r\n"
						+ "\r\n"
						+ "Arguments:\r\n"
						+ "\r\n"
						+ "    -I, --in\r\n"
						+ "          the file of sample list(required)\r\n"
						+ "    -O, --out\r\n"
						+ "          output differential expression result(required)\r\n"
						+ "    -BM, --bsj_matrix\r\n"
						+ "          circRNA BSJ expression matrix (required)\r\n"
						+ "    -FM, --fsj_matrix\r\n"
						+ "          circRNA FSJ expression matrix (required)\r\n"
						+ "    -T, --thread_num\r\n"
						+ "          set number of threads for parallel running (default: 1)\r\n"
						+ "    -H, --help\r\n"
						+ "          show this help information\r\n");
			}else if (parameterMap.get("-I").equals("F") && parameterMap.get("--in").equals("F")) {
				System.out.println("Please use --in or -I option to designate input the file of sample list!");
			}else if (parameterMap.get("-O").equals("F") && parameterMap.get("--out").equals("F")) {
				System.out.println("Please use --out or -O option to designate output file!");
			}else if (parameterMap.get("-BM").equals("F") && parameterMap.get("--bsj_matrix").equals("F")) {
				System.out.println("Please use --bsj_matrix or -BM option to designate input the file of circRNA BSJ expression!");
			}else if (parameterMap.get("-FM").equals("F") && parameterMap.get("--fsj_matrix").equals("F")) {
				System.out.println("Please use --fsj_matrix or -FM option to designate input the file of circRNA FSJ expression!");
			}else {
				//UserGivecircRNA	
				//input
				if (!parameterMap.get("-I").equals("F") || !parameterMap.get("--in").equals("F")) {
					if(!parameterMap.get("-I").equals("F")) {
						inforPath = parameterMap.get("-I");
					}else {
						inforPath = parameterMap.get("--in");
					}
					File inputFileTem = new File(inforPath); 
					inforPath = inputFileTem.getCanonicalPath();
				}
				File infile = new File(inforPath);
				if (!infile.isFile()) {
					dieReasonList.add("No sample list file found at designated directory!");
				}
				//output
				if (!parameterMap.get("-O").equals("F") || !parameterMap.get("--out").equals("F")) {
					if(!parameterMap.get("-O").equals("F")) {
						outPath = parameterMap.get("-O");
					}else {
						outPath = parameterMap.get("--out");
					}	
					File outputFileTem = new File(outPath); 
					outPath = outputFileTem.getCanonicalPath();
				}
				File outfile = new File(outPath);
				if (!outfile.isFile()) {
					dieReasonList.add("Output file cannot be written in the directory "+outPath+"!");
				}
				//bSJmatrix
				if (!parameterMap.get("-BM").equals("F") || !parameterMap.get("--bsj_matrix").equals("F")) {
					if(!parameterMap.get("-BM").equals("F")) {
						BSJMatrixPath = parameterMap.get("-BM");
					}else {
						BSJMatrixPath = parameterMap.get("--bsj_matrix");
					}	
					File BSJMatrixPathTem = new File(BSJMatrixPath); 
					BSJMatrixPath = BSJMatrixPathTem.getCanonicalPath();
					File BSJMatrixfile = new File(BSJMatrixPath);
					if (!BSJMatrixfile.isFile()) {
						dieReasonList.add("No circRNA BSJ expression matrix file found at designated directory!");
					}
				}				
				//FSJmatrix
				if (!parameterMap.get("-FM").equals("F") || !parameterMap.get("--fsj_matrix").equals("F")) {
					if(!parameterMap.get("-FM").equals("F")) {
						FSJMatrixPath = parameterMap.get("-FM");
					}else {
						FSJMatrixPath = parameterMap.get("--fsj_matrix");
					}	
					File FSJMatrixPathTem = new File(FSJMatrixPath); 
					FSJMatrixPath = FSJMatrixPathTem.getCanonicalPath();
					File FSJMatrixfile = new File(FSJMatrixPath);
					if (!FSJMatrixfile.isFile()) {
						dieReasonList.add("No circRNA BSJ expression matrix file found at designated directory!");
					}
				}	
				//threadNum
				if(!parameterMap.get("-T").equals("1")) {
					thread = Integer.valueOf(parameterMap.get("-T"));		
				}else if(!parameterMap.get("--thread_num").equals("1")) {
					thread = Integer.valueOf(parameterMap.get("--thread_num"));
				}else {
					thread = 1;
				}
				//运行
				Ratio r = new Ratio();
				r.DERatio(inforPath, BSJMatrixPath, FSJMatrixPath, outPath, thread);    	
	        }
		}else if (Parameters[0].equals("DE_Relative")) {
			String inforPath = null,outPath = null,matrixPath = null,genePath = null;
			int thread;
			HashMap<String,String> parameterMap = new HashMap<String,String>();
			parameterMap.put("-I", "F");
			parameterMap.put("-O", "F");
			parameterMap.put("-M", "F");
			parameterMap.put("-GC", "F");
			parameterMap.put("-T", "1");
			parameterMap.put("-H", "F");
			parameterMap.put("--in", "F");
			parameterMap.put("--out", "F");
			parameterMap.put("--matrix", "F");
			parameterMap.put("--circ_gene", "F");
			parameterMap.put("--thread_num", "1");
			parameterMap.put("--help", "F");
			
			
			for (int i = 1; i < Parameters.length-1; i++) {
				if (Parameters[i].startsWith("-")) {
					parameterMap.put(Parameters[i], Parameters[i+1]);
				}
			}
			if (Parameters[Parameters.length-1].startsWith("-")) {
				parameterMap.put(Parameters[Parameters.length-1], "T");
			}
						
			if (!parameterMap.get("-H").equals("F")) {
				System.out.println(
						"without replicate"
						+ "Usage:    java CIRI3.jar DE_Relative -I infor.tsv -O output.txt \r\n"
						+ "\r\n"
						+ "Arguments:\r\n"
						+ "\r\n"
						+ "    -I, --in\r\n"
						+ "          the file of sample list(required)\r\n"
						+ "    -O, --out\r\n"
						+ "          output differential expression result(required)\r\n"
						+ "    -M, --matrix\r\n"
						+ "          circRNA BSJ expression matrix (optional)\r\n"
						+ "    -GC, --circ_gene\r\n"
						+ "          the gene information file corresponding to circRNAs (optional)\r\n"
						+ "    -T, --thread_num\r\n"
						+ "          set number of threads for parallel running (default: 1)\r\n"
						+ "    -H, --help\r\n"
						+ "          show this help information\r\n");
			}else if (parameterMap.get("-I").equals("F") && parameterMap.get("--in").equals("F")) {
				System.out.println("Please use --in or -I option to designate input the file of sample list!");
			}else if (parameterMap.get("-O").equals("F") && parameterMap.get("--out").equals("F")) {
				System.out.println("Please use --out or -O option to designate output file!");
			}else {
				//input
				if (!parameterMap.get("-I").equals("F") || !parameterMap.get("--in").equals("F")) {
					if(!parameterMap.get("-I").equals("F")) {
						inforPath = parameterMap.get("-I");
					}else {
						inforPath = parameterMap.get("--in");
					}
					File inputFileTem = new File(inforPath); 
					inforPath = inputFileTem.getCanonicalPath();
				}
				/*File infile = new File(inforPath);
				if (!infile.isFile()) {
					dieReasonList.add("No sample list file found at designated directory!");
				}*/
				//output
				if (!parameterMap.get("-O").equals("F") || !parameterMap.get("--out").equals("F")) {
					if(!parameterMap.get("-O").equals("F")) {
						outPath = parameterMap.get("-O");
					}else {
						outPath = parameterMap.get("--out");
					}	
					File outputFileTem = new File(outPath); 
					outPath = outputFileTem.getCanonicalPath();
				}
				File outfile = new File(outPath);
				if (!outfile.isFile()) {
					dieReasonList.add("Output file cannot be written in the directory "+outPath+"!");
				}
				//matrix
				if (!parameterMap.get("-M").equals("F") || !parameterMap.get("--matrix").equals("F")) {
					if(!parameterMap.get("-M").equals("F")) {
						matrixPath = parameterMap.get("-M");
					}else {
						matrixPath = parameterMap.get("--matrix");
					}	
					File matrixPathTem = new File(matrixPath); 
					matrixPath = matrixPathTem.getCanonicalPath();
					File matrixfile = new File(matrixPath);
					if (!matrixfile.isFile()) {
						dieReasonList.add("No circRNA BSJ expression matrix file found at designated directory!");
					}
				}		
				//gene
				if (!parameterMap.get("-GC").equals("F") || !parameterMap.get("--circ_gene").equals("F")) {
					if(!parameterMap.get("-GC").equals("F")) {
						genePath = parameterMap.get("-GC");
					}else {
						genePath = parameterMap.get("--circ_gene");
					}	
					File genePathTem = new File(genePath); 
					genePath = genePathTem.getCanonicalPath();
					File genefile = new File(genePath);
					if (!genefile.isFile()) {
						dieReasonList.add("No gene information file file found at designated directory!");
					}
				}		
				//threadNum
				if(!parameterMap.get("-T").equals("1")) {
					thread = Integer.valueOf(parameterMap.get("-T"));		
				}else if(!parameterMap.get("--thread_num").equals("1")) {
					thread = Integer.valueOf(parameterMap.get("--thread_num"));
				}else {
					thread = 1;
				}
				//运行
				Relative r = new Relative();
				if(genePath == null &&  matrixPath == null) {					
					r.DERelativeFile(inforPath,outPath, thread);
				}else {
					r.DERelativeMatrix(inforPath,matrixPath,genePath, outPath, thread);
				}
				    	
	        }
	        
		}else {
			String inputFile = null,outputFile = null,UserGivecircRNA = "",annotationFile,faFile = null,mitochondrion;
			int minMapqUni,maxCircle,minCircle,linear_range_size_min = 50000,strigency,relExp,threadNum,way;
			boolean intronLable,mLable,spLable,maLable;
			HashMap<String,String> parameterMap = new HashMap<String,String>();
			parameterMap.put("-I", "F");
			parameterMap.put("-O", "F");
			parameterMap.put("-F", "F");
			parameterMap.put("-A", "F");
			parameterMap.put("-H", "F");
			parameterMap.put("-Max", "200000");
			parameterMap.put("-Min", "140");
			parameterMap.put("-S", "2");
			parameterMap.put("-U", "10");
			parameterMap.put("-E", "0");
			parameterMap.put("-Mc", "0");
			parameterMap.put("-M", "chrM");
			parameterMap.put("-T", "1");			
			parameterMap.put("-It", "0");	
			parameterMap.put("-W", "0");
			parameterMap.put("-UC", "");
			parameterMap.put("-Sp", "0");//剪切信号 splicing signals
			parameterMap.put("-Ma", "0");//比对 map
			parameterMap.put("--in", "F");
			parameterMap.put("--out", "F");
			parameterMap.put("--ref_file", "F");
			parameterMap.put("--anno", "F");
			parameterMap.put("--help", "F");
			parameterMap.put("--max_span", "200000");
			parameterMap.put("--min_span", "140");
			parameterMap.put("--strigency", "2");
			parameterMap.put("--mapq_uni", "10");
			parameterMap.put("--rel_exp", "0");
			parameterMap.put("--mitochondria", "0");
			parameterMap.put("--chrM", "chrM");
			parameterMap.put("--thread_num", "1");		
			parameterMap.put("--intron", "0");
			parameterMap.put("--way", "0");
			parameterMap.put("--user_circ", "");
			parameterMap.put("--splicing_signals", "0");
			parameterMap.put("--mapper", "0");
			for (int i = 0; i < Parameters.length-1; i++) {
				if (Parameters[i].startsWith("-")) {
					parameterMap.put(Parameters[i], Parameters[i+1]);
				}
			}
			if (Parameters[Parameters.length-1].startsWith("-")) {
				parameterMap.put(Parameters[Parameters.length-1], "T");
			}
			
			if (!parameterMap.get("-H").equals("F")) {
				System.out.println("Program:  CIRI3 (circRNA identifier3)\r\n"
						+ "Version:  $version\r\n"
						+ "Contact:  Xin Zheng <zhengxin@big.ac.cn>\r\n"
						+ "\r\n"
						+ "Usage:    java CIRI3.jar -I in.sam -O output.ciri -F ref.fa \r\n"
						+ "\r\n"
						+ "Arguments:\r\n"
						+ "\r\n"
						+ "    -I, --in\r\n"
						+ "          Input SAM file name or SAM files list(required; generated by BWA-MEM)\r\n"
						+ "    -O, --out\r\n"
						+ "          Output circRNA file name(required)\r\n"
						+ "    -F, --ref_file\r\n"
						+ "          FASTA file of all reference sequences. Please make sure this file is\r\n"
						+ "          the same one provided to BWA-MEM (required).\r\n"
						+ "    -A, --anno\r\n"
						+ "          input GTF/GFF3 formatted annotation file name (optional)\r\n"
						+ "    -G, --log\r\n"
						+ "          output log file name (optional)\r\n"
						+ "    -H, --help\r\n"
						+ "          show this help information\r\n"
						+ "    -Max, --max_span\r\n"
						+ "          max spanning distance of circRNAs (default: 200000)\r\n"
						+ "    -Min, --min_span\r\n"
						+ "          min spanning distance of circRNAs (default: 140)\r\n"
						+ "    -S, --strigency\r\n"
						+ "          2: only output circRNAs supported by more than 2 distinct PCC signals (default)\r\n"
						+ "          1: only output circRNAs supported by more than 2 junction reads\r\n"
						+ "          0: output all circRNAs regardless junction read or PCC signal counts\r\n"
						+ "    -U, --mapq_uni\r\n"
						+ "          set threshold for mappqing quality of each segment of junction reads\r\n"
						+ "          (default: 10; should be within [0,30])\r\n"
						+ "    -E, --rel_exp\r\n"
						+ "          set threshold for relative expression calculated based on counts of\r\n"
						+ "          junction reads and non-junction reads (optional: e.g. 0.1)\r\n"
						+ "    -Mc, --mitochondria\r\n"
						+ "          0: Skip the recognition of  mitochondria circRNA (default)\r\n"
						+ "          1: Perform the recognition of mitochondria circRNA (the ID of mitochondrion in reference file is required)\r\n"
						+ "    -M, --chrM\r\n"
						+ "          tell CIRI3 the ID of mitochondrion in reference file(s) (default:\r\n"
						+ "          chrM)\r\n"
						+ "    -T, --thread_num\r\n"
						+ "          set number of threads for parallel running (default: 1)\r\n"
						+ "    -It, --intron\r\n"
						+ "          0: Skip the recognition of  intronic self-ligated circRNA (default)\r\n"
						+ "          1: Perform the recognition of intronic self-ligated circRNA (formatted annotation file is required)\r\n"
						+ "    -Sp, --splicing_signals\r\n"
						+ "          0: Only canonical GT-AG splice signals were considered (default)\r\n"
						+ "          1: Both canonical GT-AG and non-canonical splice signals were considered.\r\n"
						+ "    -Ma, --mapper\r\n"
						+ "          0: The SAM file generated by the BWA-MEM (default)\r\n"
						+ "          1: the SAM file generated by STAR\r\n"
						+ "    -W, --way\r\n"
						+ "          0: Input a single SAM file to identify circRNA (default)\r\n"
						+ "          1: Input a SAM files list (including the absolute path to SAM files \r\n"
						+ "             and sample name) to identify circRNA\r\n"
						+ "          2: Input a SAM files list (including the absolute path to SAM files, \r\n"
						+ "             RNase information, and sample names) to identify circRNA, this can generate \r\n" 
						+ "             circRNA confidence score\r\n"
						+ "    -UC, --user_circ\r\n"
						+ "          file of circRNA collections of interest to users \r\n");
			}else if (parameterMap.get("-I").equals("F") && parameterMap.get("--in").equals("F")) {
				System.out.println("Please use --in or -I option to designate input SAM alignment file!");
			}else if (parameterMap.get("-O").equals("F") && parameterMap.get("--out").equals("F")) {
				System.out.println("Please use --out or -O option to designate output file!");
			}else if (parameterMap.get("-F").equals("F") && parameterMap.get("--ref_file").equals("F")) {
				System.out.println("Please use --ref-file or -F to designate one file with all references in!");
			}else {
				//UserGivecircRNA
				if (!parameterMap.get("-UC").equals("") || !parameterMap.get("--user_circ").equals("")) {
					if(!parameterMap.get("-UC").equals("F")) {
						UserGivecircRNA = parameterMap.get("-UC");
					}else {
						UserGivecircRNA = parameterMap.get("--user_circ");
					}
					File UcFile = new File(UserGivecircRNA); 
					UserGivecircRNA = UcFile.getCanonicalPath();
					UcFile = new File(UserGivecircRNA);
					if (!UcFile.isFile()) {
						dieReasonList.add("No file of circRNA collections of interest to users found in designated directory!");
					}
				}			
				//input
				if (!parameterMap.get("-I").equals("F") || !parameterMap.get("--in").equals("F")) {
					if(!parameterMap.get("-I").equals("F")) {
						inputFile = parameterMap.get("-I");
					}else {
						inputFile = parameterMap.get("--in");
					}
					File inputFileTem = new File(inputFile); 
					inputFile = inputFileTem.getCanonicalPath();
				}
				File infile = new File(inputFile);
				if (!infile.isFile()) {
					dieReasonList.add("No SAM alignment file found at designated directory!");
				}
				//output
				if (!parameterMap.get("-O").equals("F") || !parameterMap.get("--out").equals("F")) {
					if(!parameterMap.get("-O").equals("F")) {
						outputFile = parameterMap.get("-O");
					}else {
						outputFile = parameterMap.get("--out");
					}	
					File outputFileTem = new File(outputFile); 
					outputFile = outputFileTem.getCanonicalPath();
				}
				
				//fa
				if (!parameterMap.get("-F").equals("F") || !parameterMap.get("--ref_file").equals("F")) {
					if(!parameterMap.get("-F").equals("F")) {
						faFile = parameterMap.get("-F");
					}else {
						faFile = parameterMap.get("--ref_file");
					}	
					File faFileTem = new File(faFile); 
					faFile = faFileTem.getCanonicalPath();
				}	
				File fafile = new File(faFile);
				if (!fafile.isFile()) {
					dieReasonList.add("No FASTA file found at designated directory!");
				}
				//annotationFile
				if(!parameterMap.get("-A").equals("F")) {
					annotationFile = parameterMap.get("-A");	
					File annotationFileTem = new File(annotationFile); 
					annotationFile = annotationFileTem.getCanonicalPath();
					File annofile = new File(annotationFile);
					if (!annofile.isFile()) {
						dieReasonList.add("No formatted annotation file found at designated directory!");
					}
				}else if(!parameterMap.get("--anno").equals("F")) {
					annotationFile = parameterMap.get("--anno");
					File annotationFileTem = new File(annotationFile); 
					annotationFile = annotationFileTem.getCanonicalPath();
					File annofile = new File(annotationFile);
					if (!annofile.isFile()) {
						dieReasonList.add("No formatted annotation file found at designated directory!");
					}
				}else {
					annotationFile = "F";
				}
				//maxCircle
				if(!parameterMap.get("-Max").equals("200000")) {
					maxCircle = Integer.valueOf(parameterMap.get("-Max"));		
				}else if(!parameterMap.get("--max_span").equals("200000")) {
					maxCircle = Integer.valueOf(parameterMap.get("--max_span"));
				}else {
					maxCircle = 200000;
				}
				if (maxCircle < 10000) {
					dieReasonList.add("Max span size on reference should be larger than 10000!");
				}
				//minCircle
				if(!parameterMap.get("-Min").equals("140")) {
					minCircle = Integer.valueOf(parameterMap.get("-Min"));		
				}else if(!parameterMap.get("--min_span").equals("140")) {
					minCircle = Integer.valueOf(parameterMap.get("--min_span"));
				}else {
					minCircle = 140;
				}
				if (minCircle < 0) {
					dieReasonList.add("The min spanning distance of circRNA should be larger than 0!");
				}
				//strigency
				if(!parameterMap.get("-S").equals("2")) {
					strigency = Integer.valueOf(parameterMap.get("-S"));		
				}else if(!parameterMap.get("--strigency").equals("2")) {
					strigency = Integer.valueOf(parameterMap.get("--strigency"));
				}else {
					strigency = 2;
				}
				if (strigency != 0 && strigency != 1 && strigency != 2) {
					dieReasonList.add("--strigency or -S can only be 0, 1, 2!");
				}
				//minMapqUni
				if(!parameterMap.get("-U").equals("10")) {
					minMapqUni = Integer.valueOf(parameterMap.get("-U"));		
				}else if(!parameterMap.get("--mapq_uni").equals("10")) {
					minMapqUni = Integer.valueOf(parameterMap.get("--mapq_uni"));
				}else {
					minMapqUni = 10;
				}
				if (minMapqUni < 0 || minMapqUni > 30) {
					dieReasonList.add("--mapq_uni or -U can only be between 0 and 30!");
				}			
				//relExp
				if(!parameterMap.get("-E").equals("0")) {
					relExp = Integer.valueOf(parameterMap.get("-E"));		
				}else if(!parameterMap.get("--rel_exp").equals("0")) {
					relExp = Integer.valueOf(parameterMap.get("--rel_exp"));
				}else {
					relExp = 0;
				}
				if (relExp < 0 || relExp > 1) {
					dieReasonList.add("--rel_exp or -E can only be between 0 and 1!");
				}
				//mitochondrion
				if(!parameterMap.get("-Mc").equals("0")) {
					mLable = true;		
				}else if(!parameterMap.get("--mitochondria").equals("0")) {
					mLable = true;
				}else {
					mLable = false;
				}
				//chrM
				if(!parameterMap.get("-M").equals("chrM")) {
					mitochondrion = parameterMap.get("-M");		
				}else if(!parameterMap.get("--chrM").equals("chrM")) {
					mitochondrion = parameterMap.get("--chrM");
				}else {
					mitochondrion = "chrM";
				}
				//threadNum
				if(!parameterMap.get("-T").equals("1")) {
					threadNum = Integer.valueOf(parameterMap.get("-T"));		
				}else if(!parameterMap.get("--thread_num").equals("1")) {
					threadNum = Integer.valueOf(parameterMap.get("--thread_num"));
				}else {
					threadNum = 1;
				}
				//intronLable
				if(!parameterMap.get("-It").equals("0")) {//不等于
					intronLable = true;		
				}else if(!parameterMap.get("--intron").equals("0")) {
					intronLable = true;		
				}else {
					intronLable = false;
				}
				//splicing_signals
				if(!parameterMap.get("-Sp").equals("0")) {
					spLable = true;		
				}else if(!parameterMap.get("--splicing_signals").equals("0")) {
					spLable = true;		
				}else {
					spLable = false;
				}
				//mapper
				if(!parameterMap.get("-Ma").equals("0")) {
					maLable = false;		
				}else if(!parameterMap.get("--mapper").equals("0")) {
					maLable = false;		
				}else {
					maLable = true;
				}
				
				//way
				if(!parameterMap.get("-W").equals("0")) {
					way = Integer.valueOf(parameterMap.get("-W"));		
				}else if(!parameterMap.get("--way").equals("0")) {
					way = Integer.valueOf(parameterMap.get("--way"));
				}else {
					way = 0;
				}
				if (way != 0 && way != 1 && way != 2) {
					dieReasonList.add("--way or -W can only be 0, 1, 2!");
				}
				//other
				if(intronLable == true && annotationFile == "F") {
					dieReasonList.add("Identification of intronic self-ligated circRNA requires formatted annotation file!");
				}
				if(dieReasonList.size()>0) {
					for (String die : dieReasonList) {
						System.out.println(die);
					}
				}else {
					if (way == 0) {
						if (threadNum == 1) {
							if(maLable) {
								SingleTest ST = new SingleTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
								ST.CIRI3(inputFile, outputFile,annotationFile, faFile,UserGivecircRNA);
							}else {
								SingleSTARTest ST = new SingleSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
								ST.CIRI3(inputFile, outputFile,annotationFile, faFile,UserGivecircRNA);
							}
							
						}else {
							if(maLable) {
								MutTest mutTest = new MutTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
								mutTest.CIRI3(inputFile, outputFile,annotationFile, faFile,threadNum,UserGivecircRNA);
							}else{
								MutSTARTest mutTest = new MutSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
								mutTest.CIRI3(inputFile, outputFile,annotationFile, faFile,threadNum,UserGivecircRNA);
							}
							
						}
					}else if (way == 1) {
						if (threadNum == 1) {
                            if(maLable) {
                            	FileTest fileTest = new FileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
    							fileTest.CIRI3(inputFile, outputFile,annotationFile, faFile,UserGivecircRNA);
							}else{
								FileSTARTest fileTest = new FileSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
    							fileTest.CIRI3(inputFile, outputFile,annotationFile, faFile,UserGivecircRNA);
							}
							
						}else {
                            if(maLable) {
                            	MutFileTest mutFileTest = new MutFileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
    							mutFileTest.CIRI3(inputFile, outputFile,annotationFile, faFile, threadNum,UserGivecircRNA);
							}else{
								MutFileSTARTest mutFileTest = new MutFileSTARTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
								mutFileTest.CIRI3(inputFile, outputFile,annotationFile, faFile, threadNum,UserGivecircRNA);
							}	
						}
					}else if (way == 2) {
						if (threadNum == 1) {
							FileTsvTest fileTsvTest = new FileTsvTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
							fileTsvTest.CIRI3(inputFile, outputFile,annotationFile, faFile,UserGivecircRNA);
						}else {
							MutTsvFileTest mutTsvFileTest = new MutTsvFileTest(minMapqUni, maxCircle, minCircle,linear_range_size_min,intronLable,strigency,relExp,mitochondrion,mLable,spLable);
							mutTsvFileTest.CIRI3(inputFile, outputFile,annotationFile, faFile, threadNum,UserGivecircRNA);
						}
					}
					
					
				}
				
				
			}
		}
		
	}

}
