package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.zx.findcircrna.BamToSam;
import com.zx.findcircrna.CircRNAUniverseIO;
import com.zx.findcircrna.MutFindCircRNAScan2;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class Scan2Test {
    private int minMapqUni, linear_range_size_min, strigency, relExp, seqLen = 0, AllFileSplitNum = 10;
    private boolean intronLable;
    private String mitochondrion;

    public Scan2Test(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
            boolean intronLable, int strigency, int relExp, String mitochondrion,
            boolean mlable, boolean spLable) {
        this.minMapqUni = minMapqUni;
        this.linear_range_size_min = linear_range_size_min;
        this.intronLable = intronLable;
        this.strigency = strigency;
        this.relExp = relExp;
        this.mitochondrion = mitochondrion;
    }

    public static String samFile;
    private static Lock lock = new ReentrantLock();
    private static HashMap<String, Integer> circFSJMap;
    private static HashMap<String, String> chrTCGAMap;
    private static HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1, chrSiteMap2;
    public static HashMap<String, byte[]> siteArrayMap1, siteArrayMap2;

    /**
     * inputFile:    SAM/BAM file (same as used in SCAN1; BSJ files must exist at {samFile}BSJ{i})
     * outputFile:   output prefix; writes {outputFile}.fsj_counts
     * universeFile: path to .universe file from BUILD_UNIVERSE
     * faFile:       FASTA reference genome
     * annotationFile: annotation file path or "F" (not used for scan2 itself, kept for API symmetry)
     * threads:      thread count (should match the count used in SCAN1)
     */
    public boolean CIRI3(String inputFile, String outputFile, String universeFile,
            String faFile, String annotationFile, int threads) throws IOException {
        long startTime = System.currentTimeMillis();
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String outputFileLog = outputFile + ".log";
        BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
        System.out.println(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN2 start");
        fileLog.write(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN2 start\n");

        // BAM->SAM conversion
        if (inputFile.substring(inputFile.length() - 3).equals("sam")) {
            samFile = inputFile;
        } else if (inputFile.substring(inputFile.length() - 3).equals("bam")) {
            samFile = inputFile.substring(0, inputFile.length() - 3) + "sam";
            BamToSam bts = new BamToSam();
            bts.bamToBam(inputFile, samFile);
        } else {
            System.out.println("Please enter the file that ends with sam or bam");
            return false;
        }

        // Step 1: Read universe file
        int[] seqLenOut = new int[1];
        HashMap<String, String[]> universeDataMap = CircRNAUniverseIO.readUniverse(universeFile, seqLenOut);
        seqLen = seqLenOut[0];
        System.out.println(df.format(System.currentTimeMillis()) + " :Universe loaded: " + universeDataMap.size() + " circRNAs, seqLen=" + seqLen);
        fileLog.write(df.format(System.currentTimeMillis()) + " :Universe loaded: " + universeDataMap.size() + " circRNAs, seqLen=" + seqLen + "\n");

        // Step 2: Load FA file
        ReadFaFile RF = new ReadFaFile();
        RF.readFa(faFile);
        HashMap<String, Integer> chrLenMap = RF.getChrLenMap();
        chrTCGAMap = RF.getChrTCGAMap();
        RF = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");
        fileLog.write(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files\n");

        // Step 3: Reconstruct circFSJMap and site-index structures from universe
        // (identical logic to MutTest lines 282-358, but reading from universeDataMap)
        circFSJMap = new HashMap<String, Integer>();
        chrSiteMap1 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        chrSiteMap2 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        siteArrayMap1 = new HashMap<String, byte[]>();
        siteArrayMap2 = new HashMap<String, byte[]>();

        // Group universe entries by chromosome
        HashMap<String, HashMap<String, String[]>> chrUniverseMap = new HashMap<String, HashMap<String, String[]>>();
        for (String circKey : universeDataMap.keySet()) {
            String chrKey = circKey.split("\t")[0];
            if (!chrUniverseMap.containsKey(chrKey)) {
                chrUniverseMap.put(chrKey, new HashMap<String, String[]>());
            }
            chrUniverseMap.get(chrKey).put(circKey, universeDataMap.get(circKey));
        }

        for (String chrKey : chrUniverseMap.keySet()) {
            HashMap<String, String[]> circMap = chrUniverseMap.get(chrKey);
            byte[] siteArray1 = new byte[(chrLenMap.get(chrKey) / seqLen) + 1];
            byte[] siteArray2 = new byte[(chrLenMap.get(chrKey) / seqLen) + 1];
            HashMap<Integer, ArrayList<SiteSort>> SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
            HashMap<Integer, ArrayList<SiteSort>> SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
            for (String circKey : circMap.keySet()) {
                String[] arr = circMap.get(circKey);
                // arr[0]=start, arr[1]=end (from universeDataMap — matches siteInfor.split("\t"))
                circFSJMap.put(circKey, 0);
                int site1 = Integer.parseInt(arr[0]) / seqLen;
                int site2 = Integer.parseInt(arr[1]) / seqLen;
                siteArray1[site1] = 1;
                siteArray2[site2] = 1;
                if (!SiteMap1.containsKey(site1)) {
                    SiteMap1.put(site1, new ArrayList<SiteSort>());
                }
                SiteMap1.get(site1).add(new SiteSort(Integer.parseInt(arr[0]), arr));
                if (!SiteMap2.containsKey(site2)) {
                    SiteMap2.put(site2, new ArrayList<SiteSort>());
                }
                SiteMap2.get(site2).add(new SiteSort(Integer.parseInt(arr[1]), arr));
            }
            // Sort site lists
            for (Integer site : SiteMap1.keySet()) {
                Collections.sort(SiteMap1.get(site));
            }
            for (Integer site : SiteMap2.keySet()) {
                Collections.sort(SiteMap2.get(site));
            }
            chrSiteMap1.put(chrKey, SiteMap1);
            chrSiteMap2.put(chrKey, SiteMap2);
            siteArrayMap1.put(chrKey, siteArray1);
            siteArrayMap2.put(chrKey, siteArray2);
        }
        chrUniverseMap = null;
        universeDataMap = null;

        // Determine AllFileSplitNum from BSJ file count
        AllFileSplitNum = 0;
        while (new File(samFile + "BSJ" + (AllFileSplitNum + 1)).exists()) {
            AllFileSplitNum++;
        }
        if (AllFileSplitNum == 0) {
            System.out.println("ERROR: No BSJ files found at " + samFile + "BSJ1. Run SCAN1 first.");
            fileLog.close();
            return false;
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :AllFileSplitNum=" + AllFileSplitNum);
        fileLog.write(df.format(System.currentTimeMillis()) + " :AllFileSplitNum=" + AllFileSplitNum + "\n");

        // Step 4: Launch thread pool and run scan2 (identical to MutTest lines 118-154)
        ExecutorService poolExe = Executors.newFixedThreadPool(threads);
        final CyclicBarrier threadSub = new CyclicBarrier(threads + 1);
        final CyclicBarrier threadMain = new CyclicBarrier(threads + 1);
        AtomicInteger incr = new AtomicInteger(1);

        for (int i = 0; i < threads; i++) {
            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        threadMain.await();
                        MutFindCircRNAScan2 scan2 = new MutFindCircRNAScan2(minMapqUni, circFSJMap,
                                linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1, chrSiteMap2,
                                chrTCGAMap, seqLen, intronLable);
                        HashMap<String, String> scan1IdMap = new HashMap<String, String>();
                        while (true) {
                            int threadNum = incr.getAndIncrement();
                            if (threadNum > AllFileSplitNum) {
                                break;
                            } else {
                                scan1IdMap.clear();
                                BufferedReader BSJbr = new BufferedReader(
                                        new FileReader(new File(samFile + "BSJ" + threadNum)));
                                String line = BSJbr.readLine();
                                while (line != null) {
                                    String[] BSJArr = line.split("\t", 2);
                                    scan1IdMap.put(BSJArr[0], "");
                                    line = BSJbr.readLine();
                                }
                                BSJbr.close();
                                scan2.findCircRNAScan2(samFile, scan1IdMap, AllFileSplitNum, threadNum);
                                System.out.println(df.format(System.currentTimeMillis()) + " :Second scan completed " + threadNum);
                                fileLog.write(df.format(System.currentTimeMillis()) + " :Second scan completed " + threadNum + "\n");
                            }
                        }
                        HashMap<String, Integer> circFSJMapTem = scan2.getCircFSJMap();
                        scan2 = null;
                        lock.lock();
                        for (String circKey : circFSJMap.keySet()) {
                            int num = circFSJMap.get(circKey);
                            int numNew = circFSJMapTem.get(circKey);
                            circFSJMap.put(circKey, num + numNew);
                        }
                        lock.unlock();
                        scan1IdMap = null;
                        circFSJMapTem = null;
                        threadSub.await();
                        threadMain.await();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            };
            poolExe.execute(runnable);
        }

        try {
            threadMain.await();
            threadMain.reset();
            threadSub.await(); // wait for scan2 done

            // Step 5: Cleanup site maps
            chrSiteMap1 = null;
            chrSiteMap2 = null;
            siteArrayMap1 = null;
            siteArrayMap2 = null;

            System.out.println(df.format(System.currentTimeMillis()) + " :Second scan completed");
            fileLog.write(df.format(System.currentTimeMillis()) + " :Second scan completed\n");

            // Step 6: Write FSJ counts
            String fsjCountsPath = outputFile + ".fsj_counts";
            CircRNAUniverseIO.writeFSJCounts(fsjCountsPath, circFSJMap);
            circFSJMap = null;
            System.out.println(df.format(System.currentTimeMillis()) + " :FSJ counts written to " + fsjCountsPath);
            fileLog.write(df.format(System.currentTimeMillis()) + " :FSJ counts written to " + fsjCountsPath + "\n");

            threadMain.await(); // release threads to exit

        } catch (Exception e) {
            e.printStackTrace();
        }

        poolExe.shutdown();
        // Do NOT delete BSJ files — FINALIZE needs them for BSJ matrix and Summary

        long endTime = System.currentTimeMillis();
        System.out.println("Program run time: " + (endTime - startTime) + "ms");
        fileLog.write("Program run time: " + (endTime - startTime) + "ms\n");
        fileLog.close();
        return true;
    }
}
