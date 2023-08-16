package DE;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.test.test;


public class BSJ_yes {
	
	public boolean DEBSJNoFile(String inforPath,String genePath,String outPath) throws IOException {		
		ArrayList<String> pathList = new ArrayList<String>();
		ArrayList<String> sampleList = new ArrayList<String>();
		
		BufferedReader infor = new BufferedReader(new FileReader(new File(inforPath)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPath)));	
		String line = infor.readLine();
		while (line != null) {
			if (line.startsWith("#")) {
				line = infor.readLine();
				continue;
			}
			String[] arr = line.split("\t");
			pathList.add(arr[1]);
			sampleList.add(arr[0]);			
			line = infor.readLine();
		}
		infor.close();
				
		HashSet<String> circSet = new HashSet<String>();
		//导入case信息
		BufferedReader br;
		for (String path : pathList) {
			br = new BufferedReader(new FileReader(new File(path)));
			line = br.readLine();
			line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				circSet.add(arr[0]);
				line = br.readLine();
			}
			br.close();	
		}
		
				
		//构建BSJ矩阵
		int[][] BSJMatrix = new int[circSet.size()][sampleList.size()];
		HashMap<String, Integer> circMap = new  HashMap<String, Integer>();
		int num = 0;
		for (String circ : circSet) {
			circMap.put(circ, num);
			num++;
		}
		
		num = 0;
		for (String path : pathList) {
			br = new BufferedReader(new FileReader(new File(path)));
			line = br.readLine();
			line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				BSJMatrix[circMap.get(arr[0])][num] = Integer.valueOf(arr[4]);
				line = br.readLine();
			}
			br.close();	
			num++;
		}
			
		//输出
		bw.write("circRNA_ID");
		for (String sample : sampleList) {
			bw.write("\t"+sample);
		}
		bw.write("\n");
		
		
			
		for (String circ : circSet) {
			bw.write(circ);
			for (int i = 0; i < sampleList.size(); i++) {
				bw.write("\t"+BSJMatrix[circMap.get(circ)][i]);			
			}	
			bw.write("\n");
		}
		bw.close();
		
		try {
			// 获取当前 JAR 文件的目录路径
            String jarDirectory = getJarDirectory();	
            // 指定R脚本文件的路径
            String rScriptPath = jarDirectory+"/scripts/BSJ_yes.R";
            // 构建命令
            String command = "Rscript " + rScriptPath +" "+inforPath+" "+outPath+" "+genePath+" "+ outPath;
            System.out.println(command);		
            // 创建进程并执行命令
            ProcessBuilder processBuilder = new ProcessBuilder(command.split(" "));
            Process process = processBuilder.start();
            // 等待进程执行完毕
            int exitCode = process.waitFor();
            // 输出执行结果
            if (exitCode == 0) {
                System.out.println("Differential expression analysis completed.");
            } else {
                System.out.println("Error executing R script. Exit code: " + exitCode);
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
		return false;
	}
	
	
	public boolean DEBSJNoMatrix(String inforPath,String genePath,String outPath,String BSJPath) throws IOException {		
		try {
			// 获取当前 JAR 文件的目录路径
            String jarDirectory = getJarDirectory();	
            // 指定R脚本文件的路径
            String rScriptPath = jarDirectory+"/scripts/BSJ_yes.R";
            // 构建命令
            String command = "Rscript " + rScriptPath +" "+inforPath+" "+BSJPath+" "+genePath+" "+ outPath;
            System.out.println(command);		
            // 创建进程并执行命令
            ProcessBuilder processBuilder = new ProcessBuilder(command.split(" "));
            Process process = processBuilder.start();
            // 等待进程执行完毕
            int exitCode = process.waitFor();
            // 输出执行结果
            if (exitCode == 0) {
                System.out.println("Differential expression analysis completed.");
            } else {
                System.out.println("Error executing R script. Exit code: " + exitCode);
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
		return false;
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
