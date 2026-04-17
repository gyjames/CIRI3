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
import com.zx.findcircrna.MutFindCircRNASTARScan2;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class Scan2STARTest {
    private int minMapqUni, linear_range_size_min, strigency, relExp, seqLen = 0, AllFileSplitNum = 10;
    private long matchNum = 0;
    private boolean intronLable;

    public Scan2STARTest(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
            boolean intronLable, int strigency, int relExp, String mitochondrion,
            boolean mlable, boolean spLable) {
        this.minMapqUni = minMapqUni;
        this.linear_range_size_min = linear_range_size_min;
        this.intronLable = intronLable;
        this.strigency = strigency;
        this.relExp = relExp;
    }

    public static String bwaSamFile, starSamFile;
    private static Lock lock = new ReentrantLock();
    private static HashMap<String, Integer> circFSJMap;
    private static HashMap<String, String> chrTCGAMap;
    private static HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1, chrSiteMap2;
    public static HashMap<String, byte[]> siteArrayMap1, siteArrayMap2;

    /**
     * inputFile:    comma-separated "chimericPath,starSamPath,unmappedSamPath" (same as SCAN1 STAR input)
     * outputFile:   output prefix; writes {outputFile}.fsj_counts
     * universeFile: path to .universe file from BUILD_UNIVERSE
     * faFile:       FASTA reference genome
     * annotationFile: "F" or annotation path (not used directly, kept for symmetry)
     * threads:      thread count
     */
    public boolean CIRI3(String inputFile, String outputFile, String universeFile,
            String faFile, String annotationFile, int threads) throws IOException {
        long startTime = System.currentTimeMillis();
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String outputFileLog = outputFile + ".log";
        BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
        System.out.println(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN2 (STAR) start");
        fileLog.write(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN2 (STAR) start\n");

        // Parse STAR input triple
        String[] samFileArr = inputFile.split(",");
        String starSamPath = samFileArr[1];
        String unmappedSamPath = samFileArr[2];

        if (unmappedSamPath.substring(unmappedSamPath.length() - 3).equals("sam")) {
            bwaSamFile = unmappedSamPath;
        } else if (unmappedSamPath.substring(unmappedSamPath.length() - 3).equals("bam")) {
            bwaSamFile = unmappedSamPath.substring(0, unmappedSamPath.length() - 3) + "sam";
            BamToSam bts = new BamToSam();
            bts.bamToBam(unmappedSamPath, bwaSamFile);
        } else {
            System.out.println("Please enter the file that ends with sam or bam");
            return false;
        }
        if (starSamPath.substring(starSamPath.length() - 3).equals("sam")) {
            starSamFile = starSamPath;
        } else if (starSamPath.substring(starSamPath.length() - 3).equals("bam")) {
            starSamFile = starSamPath.substring(0, starSamPath.length() - 3) + "sam";
            BamToSam bts = new BamToSam();
            bts.bamToBam(starSamPath, starSamFile);
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
        circFSJMap = new HashMap<String, Integer>();
        chrSiteMap1 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        chrSiteMap2 = new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        siteArrayMap1 = new HashMap<String, byte[]>();
        siteArrayMap2 = new HashMap<String, byte[]>();

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

        // Step 4: Determine AllFileSplitNum for bwaSamFile from BSJ file count
        AllFileSplitNum = 0;
        while (new File(bwaSamFile + "BSJ" + (AllFileSplitNum + 1)).exists()) {
            AllFileSplitNum++;
        }
        if (AllFileSplitNum == 0) {
            System.out.println("ERROR: No BSJ files found at " + bwaSamFile + "BSJ1. Run SCAN1 first.");
            fileLog.close();
            return false;
        }
        final int bwaFileSplitNum = AllFileSplitNum;

        // Build all-upfront idCircMap from all SCAN1 BSJ files (matches original -Ma 1 pipeline)
        HashMap<String, String> idCircMap = new HashMap<String, String>();
        long bsjTotalAfterScan1 = 0;
        long bsjTagOneAfterScan1 = 0;
        for (int i = 1; i <= bwaFileSplitNum; i++) {
            long fileLines = 0, fileTagOne = 0;
            BufferedReader BSJbr = new BufferedReader(new FileReader(new File(bwaSamFile + "BSJ" + i)));
            String bsjLine = BSJbr.readLine();
            while (bsjLine != null) {
                String[] BSJArr = bsjLine.split("\t", 2);
                idCircMap.put(BSJArr[0], "");
                fileLines++;
                String[] cols = bsjLine.split("\t", 7);
                if (cols.length > 2 && cols[2].equals("1")) fileTagOne++;
                bsjLine = BSJbr.readLine();
            }
            BSJbr.close();
            bsjTotalAfterScan1 += fileLines;
            bsjTagOneAfterScan1 += fileTagOne;
            System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_scan1 BSJ" + i + ": total=" + fileLines + " tag1=" + fileTagOne);
            fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_scan1 BSJ" + i + ": total=" + fileLines + " tag1=" + fileTagOne + "\n");
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_scan1 BSJ_TOTAL: total=" + bsjTotalAfterScan1 + " tag1=" + bsjTagOneAfterScan1);
        fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_scan1 BSJ_TOTAL: total=" + bsjTotalAfterScan1 + " tag1=" + bsjTagOneAfterScan1 + "\n");
        System.out.println(df.format(System.currentTimeMillis()) + " :idCircMap loaded: " + idCircMap.size() + " reads");
        fileLog.write(df.format(System.currentTimeMillis()) + " :idCircMap loaded: " + idCircMap.size() + " reads\n");

        ExecutorService poolExe = Executors.newFixedThreadPool(threads);
        final CyclicBarrier threadSub = new CyclicBarrier(threads + 1);
        final CyclicBarrier threadMain = new CyclicBarrier(threads + 1);
        AtomicInteger incr = new AtomicInteger(1);

        for (int i = 0; i < threads; i++) {
            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        // Phase A: BWA scan2 — appends new BSJ reads to BSJ files
                        threadMain.await();
                        MutFindCircRNAScan2 scan2 = new MutFindCircRNAScan2(minMapqUni, circFSJMap,
                                linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1, chrSiteMap2,
                                chrTCGAMap, seqLen, intronLable);
                        while (true) {
                            int threadNum = incr.getAndIncrement();
                            if (threadNum > bwaFileSplitNum) {
                                break;
                            } else {
                                scan2.findCircRNAScan2(bwaSamFile, idCircMap, bwaFileSplitNum, threadNum);
                                System.out.println(df.format(System.currentTimeMillis()) + " :unmapSam Second scan completed " + threadNum);
                                fileLog.write(df.format(System.currentTimeMillis()) + " :unmapSam Second scan completed " + threadNum + "\n");
                            }
                        }
                        HashMap<String, Integer> circFSJMapTem = scan2.getCircFSJMap();
                        lock.lock();
                        for (String circKey : circFSJMap.keySet()) {
                            circFSJMap.put(circKey, circFSJMap.get(circKey) + circFSJMapTem.get(circKey));
                        }
                        lock.unlock();
                        scan2.setFSJScan2List();
                        scan2 = null;
                        circFSJMapTem = null;
                        threadSub.await(); // [B] bwa scan2 done
                        threadMain.await(); // [C] wait for star scan2 start

                        // Phase B: STAR scan2 — appends new BSJ reads to BSJ files
                        MutFindCircRNASTARScan2 starScan2 = new MutFindCircRNASTARScan2(minMapqUni, circFSJMap,
                                linear_range_size_min, siteArrayMap1, siteArrayMap2, chrSiteMap1, chrSiteMap2,
                                chrTCGAMap, seqLen, intronLable);
                        while (true) {
                            int threadNum = incr.getAndIncrement();
                            if (threadNum > AllFileSplitNum) {
                                break;
                            } else {
                                starScan2.findCircRNAScan2(starSamFile, idCircMap, AllFileSplitNum, threadNum, bwaSamFile);
                                System.out.println(df.format(System.currentTimeMillis()) + " :starsam Second scan completed " + threadNum);
                                fileLog.write(df.format(System.currentTimeMillis()) + " :starsam Second scan completed " + threadNum + "\n");
                            }
                        }
                        circFSJMapTem = starScan2.getCircFSJMap();
                        long matchNumTem = starScan2.getReadNum();
                        lock.lock();
                        for (String circKey : circFSJMap.keySet()) {
                            circFSJMap.put(circKey, circFSJMap.get(circKey) + circFSJMapTem.get(circKey));
                        }
                        matchNum = matchNum + matchNumTem;
                        lock.unlock();
                        starScan2.setFSJScan2List();
                        starScan2 = null;
                        circFSJMapTem = null;
                        threadSub.await(); // [D] star scan2 done
                        threadMain.await(); // [E] wait for final release, then exit
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            };
            poolExe.execute(runnable);
        }

        try {
            // Release threads for BWA scan2 [A]
            threadMain.await();
            threadMain.reset();
            threadSub.await(); // [B] wait for BWA scan2 done
            System.out.println(df.format(System.currentTimeMillis()) + " :BWA unmapped scan2 completed");
            fileLog.write(df.format(System.currentTimeMillis()) + " :BWA unmapped scan2 completed\n");

            // Log BSJ counts after BWA scan2
            long bsjTotalAfterBwa = 0, bsjTagOneAfterBwa = 0;
            for (int i = 1; i <= bwaFileSplitNum; i++) {
                long fl = 0, ft1 = 0;
                BufferedReader br = new BufferedReader(new FileReader(new File(bwaSamFile + "BSJ" + i)));
                String ln = br.readLine();
                while (ln != null) { fl++; String[] c = ln.split("\t", 7); if (c.length > 2 && c[2].equals("1")) ft1++; ln = br.readLine(); }
                br.close();
                bsjTotalAfterBwa += fl; bsjTagOneAfterBwa += ft1;
                System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_bwa_scan2 BSJ" + i + ": total=" + fl + " tag1=" + ft1);
                fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_bwa_scan2 BSJ" + i + ": total=" + fl + " tag1=" + ft1 + "\n");
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_bwa_scan2 BSJ_TOTAL: total=" + bsjTotalAfterBwa + " tag1=" + bsjTagOneAfterBwa);
            fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_bwa_scan2 BSJ_TOTAL: total=" + bsjTotalAfterBwa + " tag1=" + bsjTagOneAfterBwa + "\n");

            // Compute AllFileSplitNum for STAR SAM
            File starFile = new File(starSamFile);
            long fileSizeGB = starFile.length() / 1024 / 1024 / 1024;
            if (fileSizeGB > 200 && fileSizeGB > threads * 10) {
                AllFileSplitNum = (int) fileSizeGB / 10;
            } else {
                AllFileSplitNum = threads;
            }
            incr.set(1);

            // Release threads for STAR scan2 [C]
            threadSub.reset();
            threadMain.await();
            threadMain.reset();
            threadSub.await(); // [D] wait for STAR scan2 done

            System.out.println(df.format(System.currentTimeMillis()) + " :STAR scan2 completed");
            System.out.println(df.format(System.currentTimeMillis()) + " :Mapped_Reads " + matchNum);
            fileLog.write(df.format(System.currentTimeMillis()) + " :STAR scan2 completed\n");
            fileLog.write(df.format(System.currentTimeMillis()) + " :Mapped_Reads " + matchNum + "\n");

            // Log BSJ counts after STAR scan2
            long bsjTotalAfterStar = 0, bsjTagOneAfterStar = 0;
            for (int i = 1; i <= AllFileSplitNum; i++) {
                long fl = 0, ft1 = 0;
                BufferedReader br = new BufferedReader(new FileReader(new File(bwaSamFile + "BSJ" + i)));
                String ln = br.readLine();
                while (ln != null) { fl++; String[] c = ln.split("\t", 7); if (c.length > 2 && c[2].equals("1")) ft1++; ln = br.readLine(); }
                br.close();
                bsjTotalAfterStar += fl; bsjTagOneAfterStar += ft1;
                System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_star_scan2 BSJ" + i + ": total=" + fl + " tag1=" + ft1);
                fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_star_scan2 BSJ" + i + ": total=" + fl + " tag1=" + ft1 + "\n");
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :DIAG after_star_scan2 BSJ_TOTAL: total=" + bsjTotalAfterStar + " tag1=" + bsjTagOneAfterStar);
            fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG after_star_scan2 BSJ_TOTAL: total=" + bsjTotalAfterStar + " tag1=" + bsjTagOneAfterStar + "\n");

            // Cleanup site maps
            chrSiteMap1 = null;
            chrSiteMap2 = null;
            siteArrayMap1 = null;
            siteArrayMap2 = null;

            // Write FSJ counts
            String fsjCountsPath = outputFile + ".fsj_counts";
            CircRNAUniverseIO.writeFSJCounts(fsjCountsPath, circFSJMap);
            circFSJMap = null;
            System.out.println(df.format(System.currentTimeMillis()) + " :FSJ counts written to " + fsjCountsPath);
            fileLog.write(df.format(System.currentTimeMillis()) + " :FSJ counts written to " + fsjCountsPath + "\n");

            threadMain.await(); // [E] release threads to exit

        } catch (Exception e) {
            e.printStackTrace();
        }

        poolExe.shutdown();

        long endTime = System.currentTimeMillis();
        System.out.println("Program run time: " + (endTime - startTime) + "ms");
        fileLog.write("Program run time: " + (endTime - startTime) + "ms\n");
        fileLog.close();
        return true;
    }
}
