package DE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.test.test;



public class Relative {
	
	public void DERelativeFile(String inforPath,String outPath,int thread) throws IOException {
	//分组
	BufferedReader lab = new BufferedReader(new FileReader(new File(inforPath)));
	HashMap<String, ArrayList<String>> groupMap = new HashMap<String, ArrayList<String>>();	
	HashMap<String, String> samplePathMap = new HashMap<String, String>();
	ArrayList<String> sampleTemList = new ArrayList<String>();
	String line = lab.readLine();
	line = lab.readLine();
	while (line != null) {
		String[] arr = line.split("\t");
		if(groupMap.containsKey(arr[2])) {
			sampleTemList = groupMap.get(arr[2]);
			sampleTemList.add(arr[0]);
			groupMap.put(arr[2], sampleTemList);
		}else {
			sampleTemList = new ArrayList<String>();
			sampleTemList.add(arr[0]);
			groupMap.put(arr[2], sampleTemList);
		}
		samplePathMap.put(arr[0], arr[1]);
		line = lab.readLine();
	}
	lab.close();
	
	
	//提取信息
	HashMap<String, HashMap<String, int[]>> allMap = new HashMap<String, HashMap<String, int[]>>();
	HashMap<String, int[]> temBSJMap = new HashMap<String, int[]>();
	
	HashMap<String, HashSet<String>> geneCircMap = new HashMap<String, HashSet<String>>();
	HashSet<String> circSet = new HashSet<String>();
	
	ArrayList<String> classList = new ArrayList<>();
	classList.addAll(groupMap.keySet());
	for (int k = 0; k < classList.size()-1; k++) {
		for (int j = k+1; j < classList.size(); j++) {
			
			String classStr1 = classList.get(k);
			String classStr2 = classList.get(j);
			String outPathTem = outPath+"_"+classStr1+"_"+classStr2;
			ArrayList<String> caseSampleList = groupMap.get(classStr1);
			ArrayList<String> controlSampleList = groupMap.get(classStr2);
			//case		
			for (String sample : caseSampleList) {
				String path = samplePathMap.get(sample);
				temBSJMap = new HashMap<String, int[]>();
				BufferedReader br = new BufferedReader(new FileReader(new File(path)));						
				line = br.readLine();
				line = br.readLine();
				while (line != null) {
					String[] arr = line.split("\t");
					if(arr[9].equals("NA")) {
						line = br.readLine();
						continue;
					}
					String gene = arr[9];
					String idKey = arr[0]+";"+gene;
					if (!geneCircMap.containsKey(gene)) {
						circSet = new HashSet<String>();
						circSet.add(idKey);
						geneCircMap.put(gene, circSet);
					}else {					
						circSet = geneCircMap.get(gene);
						circSet.add(idKey);
						geneCircMap.put(gene, circSet);				
					}
					int[] temArr = {Integer.valueOf(arr[4]),(Integer.valueOf(arr[3])-Integer.valueOf(arr[2])+1)};
					temBSJMap.put(idKey, temArr);		
					line = br.readLine();
				}
				br.close();
				allMap.put(sample, temBSJMap);
			}
			//control	
			for (String sample : controlSampleList) {
				String path = samplePathMap.get(sample);
				temBSJMap = new HashMap<String, int[]>();
				BufferedReader br = new BufferedReader(new FileReader(new File(path)));						
				line = br.readLine();
				line = br.readLine();
				while (line != null) {
					String[] arr = line.split("\t");
					if(arr[9].equals("NA")) {
						line = br.readLine();
						continue;
					}
					String gene = arr[9];
					String idKey = arr[0]+";"+gene;
					if (!geneCircMap.containsKey(gene)) {
						circSet = new HashSet<String>();
						circSet.add(idKey);
						geneCircMap.put(gene, circSet);
					}else {					
						circSet = geneCircMap.get(gene);
						circSet.add(idKey);
						geneCircMap.put(gene, circSet);				
					}
					int[] temArr = {Integer.valueOf(arr[4]),(Integer.valueOf(arr[3])-Integer.valueOf(arr[2])+1)};
					temBSJMap.put(idKey, temArr);		
					line = br.readLine();
				}
				br.close();
				allMap.put(sample, temBSJMap);
			}
			
			//先计算每一个基因总和，再分别输出circRNA
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPathTem)));
			bw.write("circRNA_ID"+"\t"+"Case_BSJ"+"\t"+"Gene_Case_BSJ"+"\t"+"Control_BSJ"+"\t"+"Gene_Control_BSJ"+"\t"+"circRNALen1"+"\t"+"circRNALen2"+"\n");
			int caseSampleNum = caseSampleList.size();
			int controlSampleNum = controlSampleList.size();
			ArrayList<String> sampleList = new ArrayList<String>();
			sampleList.addAll(caseSampleList);
			sampleList.addAll(controlSampleList);		
			
			ArrayList<String> temList = new ArrayList<String>();
			ArrayList<String> temCaseSumList = new ArrayList<String>();
			ArrayList<String> temControlSumList = new ArrayList<String>();
			for (String gene : geneCircMap.keySet()) {
				circSet = geneCircMap.get(gene);
				int[] temSumBSJArr = new int[sampleList.size()];		
				for (String circ : circSet) {
					for (int i = 0; i < sampleList.size(); i++) {
						temBSJMap = allMap.get(sampleList.get(i));
						if (temBSJMap.containsKey(circ)) {
							temSumBSJArr[i] += temBSJMap.get(circ)[0];
						}
					}		
				}
				temCaseSumList.clear();
				for (int i = 0; i < caseSampleNum; i++) {
					temCaseSumList.add(temSumBSJArr[i]+"");
				}
				temControlSumList.clear();
				for (int i = caseSampleNum; i < caseSampleNum+controlSampleNum; i++) {
					temControlSumList.add(temSumBSJArr[i]+"");
				}
				String caseSum = String.join(",", temCaseSumList);
				String controlSum = String.join(",", temControlSumList);
				int circLen = 0;
				for (String circ : circSet) {
					bw.write(circ+"\t");
					int[] temBSJArr = new int[sampleList.size()];
					for (int i = 0; i < sampleList.size(); i++) {
						temBSJMap = allMap.get(sampleList.get(i));
						if (temBSJMap.containsKey(circ)) {
							temBSJArr[i] += temBSJMap.get(circ)[0];
							circLen = temBSJMap.get(circ)[1];
						}
					}
					temList.clear();
					for (int i = 0; i < caseSampleNum; i++) {
						temList.add(temBSJArr[i]+"");
					}
					bw.write(String.join(",", temList)+"\t"+caseSum+"\t");
					
					temList.clear();
					for (int i = caseSampleNum; i < caseSampleNum+controlSampleNum; i++) {
						temList.add(temBSJArr[i]+"");
					}
					bw.write(String.join(",", temList)+"\t"+controlSum+"\t"+circLen+"\t"+circLen+"\n");
				}
			}	
			bw.close();
			
			try {
				// 获取当前 JAR 文件的目录路径
		        String jarDirectory = getJarDirectory();		
		        // 设置命令及其参数
		        System.out.println(jarDirectory+"/scripts/rMATSexe"+" -i "+ outPathTem+ " -t "+ thread+ " -o "+ outPathTem+"P-V.txt"+ " -c "+ "0.001");
		        ProcessBuilder processBuilder = new ProcessBuilder(jarDirectory+"/scripts/rMATSexe", "-i", outPathTem, "-t", thread+"", "-o", outPathTem+"P-V.txt", "-c", "0.001");

		        // 启动进程
		        Process process = processBuilder.start();

		        // 获取进程的输出流（标准输出）
		        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
		        while ((line = reader.readLine()) != null) {
		            System.out.println(line);
		        }

		        // 等待命令执行完成
		        int exitCode = process.waitFor();
		        //System.out.println("Differential expression analysis completed.");

		        	
		        // 指定R脚本文件的路径
		        String rScriptPath = jarDirectory+"/scripts/FDR.R";
		        // 构建命令
		        String command = "Rscript " + rScriptPath +" "+outPathTem+"P-V.txt"+" "+outPathTem;
		        System.out.println(command);		
		        // 创建进程并执行命令
		        processBuilder = new ProcessBuilder(command.split(" "));
		        process = processBuilder.start();
		        // 等待进程执行完毕
		        exitCode = process.waitFor();
		        // 输出执行结果
		        if (exitCode == 0) {
		            System.out.println("Differential expression analysis completed.");
		        } else {
		            System.out.println("Error executing R script. Exit code: " + exitCode);
		        }
		        new File(outPathTem+"P-V.txt").delete();
		    } catch (IOException | InterruptedException e) {
		        e.printStackTrace();
		    }
		}
	}
	
 }
	public void DERelativeMatrix(String inforPath,String BSJMatrixPath,String circGnenPath,String outPath,int thread) throws IOException {
		// TODO Auto-generated method stub
		//分组
		BufferedReader lab = new BufferedReader(new FileReader(new File(inforPath)));
		HashMap<String, ArrayList<String>> groupMap = new HashMap<String, ArrayList<String>>();	
		ArrayList<String> sampleTemList = new ArrayList<String>();
		String line = lab.readLine();
		line = lab.readLine();
		while (line != null) {
			String[] arr = line.split("\t");
			if(groupMap.containsKey(arr[2])) {
				sampleTemList = groupMap.get(arr[2]);
				sampleTemList.add(arr[0]);
				groupMap.put(arr[2], sampleTemList);
			}else {
				sampleTemList = new ArrayList<String>();
				sampleTemList.add(arr[0]);
				groupMap.put(arr[2], sampleTemList);
			}
			line = lab.readLine();
		}
		lab.close();
		
		
		//提取gene信息
		HashMap<String, HashSet<String>> geneCircMap = new HashMap<String, HashSet<String>>();
		HashSet<String> circSet = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(circGnenPath)));
		line = br.readLine();
		line = br.readLine();
		while (line != null) {
			String[] arr = line.split("\t");
			if (!geneCircMap.containsKey(arr[1])) {
				circSet = new HashSet<String>();
				circSet.add(arr[0]);
				geneCircMap.put(arr[1], circSet);
			}else {					
				circSet = geneCircMap.get(arr[1]);
				circSet.add(arr[0]);
				geneCircMap.put(arr[1], circSet);				
			}
			line = br.readLine();
		}
		br.close();
		
		//提取BSJ信息
		HashMap<String, Integer> sampleSiteMap = new HashMap<String, Integer>();
		HashMap<String, String[]> circBSJMap = new HashMap<String, String[]>();
		
		br = new BufferedReader(new FileReader(new File(BSJMatrixPath)));
		line = br.readLine();
		String[] arr = line.split("\t");
		for (int i = 1; i < arr.length; i++) {
			sampleSiteMap.put(arr[i], i-1);
		}
		line = br.readLine();
		while (line != null) {
			arr = line.split("\t",2);
			String[] BSJArr = arr[1].split("\t");
			circBSJMap.put(arr[0], BSJArr);
			line = br.readLine();
		}
		br.close();
				
		
		ArrayList<String> classList = new ArrayList<>();
		classList.addAll(groupMap.keySet());
		for (int k = 0; k < classList.size()-1; k++) {
			for (int j = k+1; j < classList.size(); j++) {
				String classStr1 = classList.get(k);
				String classStr2 = classList.get(j);
				String outPathTem = outPath+"_"+classStr1+"_"+classStr2;
				ArrayList<String> caseSampleList = groupMap.get(classStr1);
				ArrayList<String> controlSampleList = groupMap.get(classStr2);
				//先计算每一个基因总和，再分别输出circRNA
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPathTem)));
				bw.write("circRNA_ID"+"\t"+"Case_BSJ"+"\t"+"Gene_Case_BSJ"+"\t"+"Control_BSJ"+"\t"+"Gene_Control_BSJ"+"\t"+"circRNALen1"+"\t"+"circRNALen2"+"\n");
				ArrayList<String> sampleList = new ArrayList<String>();
				sampleList.addAll(caseSampleList);
				sampleList.addAll(controlSampleList);		
				
				ArrayList<String> temList = new ArrayList<String>();
				ArrayList<String> temCaseSumList = new ArrayList<String>();
				ArrayList<String> temControlSumList = new ArrayList<String>();
				for (String gene : geneCircMap.keySet()) {
					circSet = geneCircMap.get(gene);
					int[] temSumBSJArr = new int[sampleList.size()];		
					for (String circ : circSet) {
						String[] BSJArr = circBSJMap.get(circ);
						for (int i = 0; i < sampleList.size(); i++) {					
							temSumBSJArr[i] += Integer.valueOf(BSJArr[i]);					
						}		
					}
					temCaseSumList.clear();
					for (String sample : caseSampleList) {
						temCaseSumList.add(temSumBSJArr[sampleSiteMap.get(sample)]+"");
					}
					temControlSumList.clear();
					for (String sample : controlSampleList) {
						temControlSumList.add(temSumBSJArr[sampleSiteMap.get(sample)]+"");
					}
					
					String caseSum = String.join(",", temCaseSumList);
					String controlSum = String.join(",", temControlSumList);
					int circLen = 0;
					for (String circ : circSet) {
					    String[] arr1 = circ.split(":");
					    String[] arr2 = arr1[1].split("\\|");
					    circLen = Integer.valueOf(arr2[1])-Integer.valueOf(arr2[0])+1;
						bw.write(circ+";"+gene+"\t");
						String[] temBSJArr = circBSJMap.get(circ);
						temList.clear();
						for (String sample : caseSampleList) {
							temList.add(temBSJArr[sampleSiteMap.get(sample)]+"");
						}
						bw.write(String.join(",", temList)+"\t"+caseSum+"\t");
						temList.clear();
						for (String sample : controlSampleList) {
							temList.add(temBSJArr[sampleSiteMap.get(sample)]+"");
						}
					    bw.write(String.join(",", temList)+"\t"+controlSum+"\t"+circLen+"\t"+circLen+"\n");
					}
				}	
				bw.close();
				try {
					// 获取当前 JAR 文件的目录路径
		            String jarDirectory = getJarDirectory();		
		            // 设置命令及其参数
		            System.out.println(jarDirectory+"/scripts/rMATSexe"+" -i "+ outPathTem+ " -t "+ thread+ " -o "+ outPathTem+"P-V.txt"+ " -c "+ "0.001");
		            ProcessBuilder processBuilder = new ProcessBuilder(jarDirectory+"/scripts/rMATSexe", "-i", outPathTem, "-t", thread+"", "-o", outPathTem+"P-V.txt", "-c", "0.001");

		            // 启动进程
		            Process process = processBuilder.start();

		            // 获取进程的输出流（标准输出）
		            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
		            while ((line = reader.readLine()) != null) {
		                System.out.println(line);
		            }

		            // 等待命令执行完成
		            int exitCode = process.waitFor();
		            //System.out.println("Differential expression analysis completed.");

		            	
		            // 指定R脚本文件的路径
		            String rScriptPath = jarDirectory+"/scripts/FDR.R";
		            // 构建命令
		            String command = "Rscript " + rScriptPath +" "+outPathTem+"P-V.txt"+" "+outPathTem;
		            System.out.println(command);		
		            // 创建进程并执行命令
		            processBuilder = new ProcessBuilder(command.split(" "));
		            process = processBuilder.start();
		            // 等待进程执行完毕
		            exitCode = process.waitFor();
		            // 输出执行结果
		            if (exitCode == 0) {
		                System.out.println("Differential expression analysis completed.");
		            } else {
		                System.out.println("Error executing R script. Exit code: " + exitCode);
		            }
		            new File(outPathTem+"P-V.txt").delete();
		        } catch (IOException | InterruptedException e) {
		            e.printStackTrace();
		        }
			}
		}
		
	 }
	private static String getJarDirectory() {
        try {
            String jarFilePath = test.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath();
            File jarFile = new File(jarFilePath);
            String jarDirectory = jarFile.getParentFile().getAbsolutePath();
            return jarDirectory;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
}
