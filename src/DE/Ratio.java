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

import com.zx.test.test;



public class Ratio {

	
	public void DERatio(String inforPath,String BSJMatrixPath,String FSJMatrixPath,String outPath,int thread) throws IOException {
		//分组
		BufferedReader infor = new BufferedReader(new FileReader(new File(inforPath)));
		HashMap<String, ArrayList<String>> groupMap = new HashMap<String, ArrayList<String>>();	
		ArrayList<String> sampleList = new ArrayList<String>();
		groupMap.put("Case", new ArrayList<String>());
		groupMap.put("Control", new ArrayList<String>());
		String line = infor.readLine();
		while (line != null) {
			if (line.startsWith("#")) {
				line = infor.readLine();
				continue;
			}
			String[] arr = line.split("\t");
			if (arr[1].equals("Case")) {
				sampleList = groupMap.get(arr[1]);
				sampleList.add(arr[0]);
				groupMap.put(arr[1], sampleList);
			}else {					
				sampleList = groupMap.get(arr[1]);
				sampleList.add(arr[0]);
				groupMap.put(arr[1], sampleList);				
			}
			line = infor.readLine();
		}
		infor.close();
		
		//存放circRNA
		ArrayList<String> circRNAList = new ArrayList<String>();

		ArrayList<String> group1List = groupMap.get("Case");
		ArrayList<String> group2List = groupMap.get("Control");
		//标记那些列是case，那些列是control
		ArrayList<Integer> group1SiteList = new ArrayList<Integer>();
		ArrayList<Integer> group2SiteList = new ArrayList<Integer>();
		BufferedReader br = new BufferedReader(new FileReader(new File(BSJMatrixPath)));	
		line = br.readLine();
		String[] arr = line.split("\t");
		for (int i = 1; i < arr.length; i++) {
			if (group1List.contains(arr[i])) {
				group1SiteList.add(i);
			}else if (group2List.contains(arr[i])) {
				group2SiteList.add(i);
			}
		}
		//整理成字典
		HashMap<String, String> group1BSJMap = new HashMap<String, String>();	
		HashMap<String, String> group2BSJMap = new HashMap<String, String>();	
		ArrayList<String> temList = new ArrayList<String>();
		line = br.readLine();
		while (line != null) {
			arr = line.split("\t");
			circRNAList.add(arr[0]);
			temList.clear();
			for (Integer site : group1SiteList) {
				temList.add((Integer.valueOf(arr[site])*2)+"");
			}
			group1BSJMap.put(arr[0], String.join(",", temList));
			temList.clear();
			for (Integer site : group2SiteList) {
				temList.add((Integer.valueOf(arr[site])*2)+"");
			}
			group2BSJMap.put(arr[0], String.join(",", temList));
			line = br.readLine();
		}
		br.close();
		
		br = new BufferedReader(new FileReader(new File(FSJMatrixPath)));	
		line = br.readLine();
		HashMap<String, String> group1FSJMap = new HashMap<String, String>();	
		HashMap<String, String> group2FSJMap = new HashMap<String, String>();	
		temList = new ArrayList<String>();
		line = br.readLine();
		while (line != null) {
			arr = line.split("\t");
			temList.clear();
			for (Integer site : group1SiteList) {
				temList.add(arr[site]);
			}
			group1FSJMap.put(arr[0], String.join(",", temList));
			temList.clear();
			for (Integer site : group2SiteList) {
				temList.add(arr[site]);
			}
			group2FSJMap.put(arr[0], String.join(",", temList));
			line = br.readLine();
		}
		br.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPath)));		
		bw.write("circRNA_ID"+"\t"+"Case_BSJ"+"\t"+"Case_FSJ"+"\t"+"Control_BSJ"+"\t"+"Control_FSJ"+"\t"+"circRNALen1"+"\t"+"circRNALen2"+"\n");
		for (String id : circRNAList) {
			String Group1_BSJ = group1BSJMap.get(id);
			String Group2_BSJ = group2BSJMap.get(id);
			String Group1_FSJ = group1FSJMap.get(id);
			String Group2_FSJ = group2FSJMap.get(id);
			arr = id.split(":");
			String[] arr1 = arr[1].split("\\|");
			int circRNALen = Integer.valueOf(arr1[1])-Integer.valueOf(arr1[0])+1;
			bw.write(id+"\t"+Group1_BSJ+"\t"+Group1_FSJ+"\t"+Group2_BSJ+"\t"+Group2_FSJ+"\t"+circRNALen+"\t"+circRNALen+"\n");
		}	
		bw.close();	
		try {
			// 获取当前 JAR 文件的目录路径
            String jarDirectory = getJarDirectory();		
            // 设置命令及其参数
            System.out.println(jarDirectory+"/scripts/rMATSexe"+" -i "+ outPath+ " -t "+ thread+ " -o "+ outPath+"P-V.txt"+ " -c "+ "0.001");
            ProcessBuilder processBuilder = new ProcessBuilder(jarDirectory+"/scripts/rMATSexe", "-i", outPath, "-t", thread+"", "-o", outPath+"P-V.txt", "-c", "0.001");

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
            String command = "Rscript " + rScriptPath +" "+outPath+"P-V.txt"+" "+outPath;
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
            new File(outPath+"P-V.txt").delete();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
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
