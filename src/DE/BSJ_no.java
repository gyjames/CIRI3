package DE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.test.test;


public class BSJ_no {

	public boolean DEBSJNoFile(String inforPath,String outPath,double pval) throws IOException {		
		String casePath = "",controlPath = "",caseName = null,controlName = null;
		double factor = 1;
		long caseMapRead = 1,conrtolMapRead = 1;
		BufferedReader infor = new BufferedReader(new FileReader(new File(inforPath)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPath)));	
		String line = infor.readLine();
		while (line != null) {
			if (line.startsWith("#")) {
				line = infor.readLine();
				continue;
			}
			String[] arr = line.split("\t");
			if(arr[2].equals("Case")) {
				casePath = arr[1];
				caseMapRead = Long.valueOf(arr[3]);
				caseName = arr[0];
			}else {
				controlPath = arr[1];
				conrtolMapRead = Long.valueOf(arr[3]);
				controlName = arr[0];
			}		
			line = infor.readLine();
		}
		infor.close();
		factor = caseMapRead / (conrtolMapRead+0.0);
		
		HashSet<String> circSet = new HashSet<String>();
		//导入case信息
		HashMap<String, String> circCaseMap = new HashMap<String, String>();
		BufferedReader caseFile = new BufferedReader(new FileReader(new File(casePath)));
		line = caseFile.readLine();
		line = caseFile.readLine();
		while (line != null) {
			String[] arr = line.split("\t");
			circCaseMap.put(arr[0], arr[4]);
			circSet.add(arr[0]);
			line = caseFile.readLine();
		}
		caseFile.close();		
		//导入control信息
		HashMap<String, String> circControlMap = new HashMap<String, String>();
		BufferedReader controlFile = new BufferedReader(new FileReader(new File(controlPath)));
		line = controlFile.readLine();
		line = controlFile.readLine();
		while (line != null) {
			String[] arr = line.split("\t");
			circControlMap.put(arr[0], arr[4]);
			circSet.add(arr[0]);
			line = controlFile.readLine();
		}
		controlFile.close();
		//构建BSJ矩阵
		bw.write("circRNA_ID"+"\t"+caseName+"\t"+controlName+"\n");
		for (String circ : circSet) {
			if (circCaseMap.containsKey(circ) && circControlMap.containsKey(circ)) {
				bw.write(circ+"\t"+circCaseMap.get(circ)+"\t"+circControlMap.get(circ)+"\n");				
			}else if (circCaseMap.containsKey(circ)) {			
				bw.write(circ+"\t"+circCaseMap.get(circ)+"\t"+"0"+"\n");
			}else {
				bw.write(circ+"\t"+"0"+"\t"+circControlMap.get(circ)+"\n");
			}
		}
		bw.close();
		
		try {
			// 获取当前 JAR 文件的目录路径
            String jarDirectory = getJarDirectory();	
            // 指定R脚本文件的路径
            String rScriptPath = jarDirectory+"/scripts/DE_Score.R";
            // 构建命令
            String command = "Rscript " + rScriptPath +" "+outPath+" "+inforPath+" "+outPath+" "+ factor+" "+ pval;
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
	
	
	public boolean DEBSJNoMatrix(String inforPath,String matrixPath,String outPath,double pval) throws IOException {		
		double factor = 1;
		long caseMapRead = 1,conrtolMapRead = 1;
		BufferedReader infor = new BufferedReader(new FileReader(new File(inforPath)));
		String line = infor.readLine();
		while (line != null) {
			if (line.startsWith("#")) {
				line = infor.readLine();
				continue;
			}
			String[] arr = line.split("\t");
			if(arr[2].equals("Case")) {
				caseMapRead = Long.valueOf(arr[3]);
			}else {
				conrtolMapRead = Long.valueOf(arr[3]);
			}		
			line = infor.readLine();
		}
		infor.close();
		factor = caseMapRead / (conrtolMapRead+0.0);
				
		try {
			// 获取当前 JAR 文件的目录路径
            String jarDirectory = getJarDirectory();		
            // 指定R脚本文件的路径
            String rScriptPath = jarDirectory+"/scripts/DE_Score.R";
            // 构建命令
            String command = "Rscript " + rScriptPath +" "+matrixPath+" "+inforPath+" "+outPath+" "+ factor+" "+ pval;
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
