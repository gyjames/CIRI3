package com.zx.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.zx.findcircrna.BamToSam;
import com.zx.findcircrna.GetAnnotationInformation;
import com.zx.findcircrna.GetChimericOut;
import com.zx.findcircrna.MutFindCircRNASTARScan1;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class Scan1STARTest {
    private int minMapqUni, maxCircle, minCircle, linear_range_size_min, strigency, relExp, seqLen = 0, AllFileSplitNum = 10;
    private boolean intronLable, mlable, spLable;
    private String mitochondrion;

    public Scan1STARTest(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
            boolean intronLable, int strigency, int relExp, String mitochondrion,
            boolean mlable, boolean spLable) {
        this.minMapqUni = minMapqUni;
        this.maxCircle = maxCircle;
        this.minCircle = minCircle;
        this.linear_range_size_min = linear_range_size_min;
        this.intronLable = intronLable;
        this.strigency = strigency;
        this.relExp = relExp;
        this.mitochondrion = mitochondrion;
        this.mlable = mlable;
        this.spLable = spLable;
    }

    public static String bwaSamFile, starSamFile;
    private static Lock lock = new ReentrantLock();
    private static HashMap<String, String> chrTCGAMap;
    private static HashMap<String, String> idCircMap;
    HashMap<String, String> chrExonStartMap = new HashMap<String, String>(),
            chrExonEndMap = new HashMap<String, String>();
    HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>(),
            chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();

    /**
     * inputFile: comma-separated "chimericPath,starSamPath,unmappedSamPath"
     */
    public boolean CIRI3(String inputFile, String outputFile, String annotationFile,
            String faFile, int threads, String UserGivecircRNA) throws IOException {
        long startTime = System.currentTimeMillis();
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String outputFileLog = outputFile + ".log";
        BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
        System.out.println(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN1 (STAR) start");
        fileLog.write(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN1 (STAR) start\n");

        // Parse STAR input triple
        String[] samFileArr = inputFile.split(",");
        String chimericPath = samFileArr[0];
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

        ExecutorService poolExe = Executors.newFixedThreadPool(threads);
        final CyclicBarrier threadSub = new CyclicBarrier(threads + 1);
        final CyclicBarrier threadMain = new CyclicBarrier(threads + 1);
        AtomicInteger incr = new AtomicInteger(1);

        for (int i = 0; i < threads; i++) {
            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        threadMain.await();
                        MutFindCircRNASTARScan1 scan1 = new MutFindCircRNASTARScan1(minMapqUni, maxCircle, minCircle,
                                linear_range_size_min, intronLable, chrExonStartMap, chrExonEndMap, chrTCGAMap,
                                chrExonStartTranscriptMap, chrExonEndTranscriptMap, mitochondrion, mlable, spLable);
                        while (true) {
                            int threadNum = incr.getAndIncrement();
                            if (threadNum > AllFileSplitNum) {
                                break;
                            } else {
                                scan1.findCircRNAScan1(bwaSamFile, AllFileSplitNum, threadNum, idCircMap);
                                System.out.println(df.format(System.currentTimeMillis()) + " :First scan completed " + threadNum);
                                fileLog.write(df.format(System.currentTimeMillis()) + " :First scan completed " + threadNum + "\n");
                            }
                        }
                        int seqLenTem = scan1.getReadLen();
                        lock.lock();
                        if (seqLenTem > seqLen) {
                            seqLen = seqLenTem;
                        }
                        lock.unlock();
                        scan1 = null;
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
            // Load annotation
            if (!annotationFile.equals("F")) {
                GetAnnotationInformation GAI = new GetAnnotationInformation();
                GAI.hand(annotationFile, intronLable);
                if (intronLable) {
                    chrExonStartTranscriptMap = GAI.getChrExonStartTranscriptMap();
                    chrExonEndTranscriptMap = GAI.getChrExonEndTranscriptMap();
                } else {
                    chrExonStartMap = GAI.getChrExonStartMap();
                    chrExonEndMap = GAI.getChrExonEndMap();
                }
                GAI = null;
                System.out.println(df.format(System.currentTimeMillis()) + " :Successfully imported comment files");
                fileLog.write(df.format(System.currentTimeMillis()) + " :Successfully imported comment files\n");
            }

            // Load FA file
            ReadFaFile RF = new ReadFaFile();
            RF.readFa(faFile);
            chrTCGAMap = RF.getChrTCGAMap();
            RF = null;
            System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");
            fileLog.write(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files\n");

            // Build idCircMap from chimeric file and write to BSJ1
            GetChimericOut getChiCirc = new GetChimericOut(maxCircle, minCircle, linear_range_size_min, intronLable,
                    chrExonStartMap, chrExonEndMap, chrTCGAMap, chrExonStartTranscriptMap, chrExonEndTranscriptMap,
                    mitochondrion, mlable, spLable);
            getChiCirc.getBSJ(chimericPath);
            idCircMap = getChiCirc.getIdCircMap();
            getChiCirc = null;
            String outBSJPath = bwaSamFile + "BSJ1";
            BufferedWriter BSJOut = new BufferedWriter(new FileWriter(new File(outBSJPath)));
            for (String idKey : idCircMap.keySet()) {
                BSJOut.write(idCircMap.get(idKey) + "\n");
            }
            BSJOut.close();

            // AllFileSplitNum = threads for STAR (bwaSamFile)
            AllFileSplitNum = threads;

            // Delete stale BSJ2-BSJn from any previous SCAN2 run to prevent contamination
            for (int k = 2; k <= AllFileSplitNum; k++) {
                new File(bwaSamFile + "BSJ" + k).delete();
            }

            // Release threads for scan1
            threadMain.await();
            threadMain.reset();
            threadSub.await(); // wait for scan1 done

            // Log BSJ file counts after scan1
            long bsjTotal = 0, bsjTagOne = 0;
            for (int k = 1; k <= AllFileSplitNum; k++) {
                java.io.File bsjFile = new java.io.File(bwaSamFile + "BSJ" + k);
                long fl = 0, ft1 = 0;
                if (bsjFile.exists()) {
                    java.io.BufferedReader bsjBr = new java.io.BufferedReader(new java.io.FileReader(bsjFile));
                    String bsjLn = bsjBr.readLine();
                    while (bsjLn != null) { fl++; String[] c = bsjLn.split("\t", 7); if (c.length > 2 && c[2].equals("1")) ft1++; bsjLn = bsjBr.readLine(); }
                    bsjBr.close();
                }
                bsjTotal += fl; bsjTagOne += ft1;
                System.out.println(df.format(System.currentTimeMillis()) + " :DIAG scan1 BSJ" + k + ": total=" + fl + " tag1=" + ft1);
                fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG scan1 BSJ" + k + ": total=" + fl + " tag1=" + ft1 + "\n");
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :DIAG scan1 BSJ_TOTAL: total=" + bsjTotal + " tag1=" + bsjTagOne + " readLen=" + seqLen);
            fileLog.write(df.format(System.currentTimeMillis()) + " :DIAG scan1 BSJ_TOTAL: total=" + bsjTotal + " tag1=" + bsjTagOne + " readLen=" + seqLen + "\n");

            // Write scan1_meta
            String metaPath = outputFile + ".scan1_meta";
            BufferedWriter metaBw = new BufferedWriter(new FileWriter(new File(metaPath)));
            metaBw.write("samFile=" + bwaSamFile + "\n");
            metaBw.write("starSamFile=" + starSamFile + "\n");
            metaBw.write("readLen=" + seqLen + "\n");
            metaBw.write("readNum=0\n");
            metaBw.write("fileSplitNum=" + AllFileSplitNum + "\n");
            metaBw.close();
            System.out.println(df.format(System.currentTimeMillis()) + " :Scan1 metadata written to " + metaPath);
            fileLog.write(df.format(System.currentTimeMillis()) + " :Scan1 metadata written to " + metaPath + "\n");

            // Release threads to exit
            threadSub.reset();
            threadMain.await();

        } catch (Exception e) {
            e.printStackTrace();
        }

        poolExe.shutdown();
        // Do NOT delete BSJ files or SAM files — needed for subsequent stages

        long endTime = System.currentTimeMillis();
        System.out.println("Program run time: " + (endTime - startTime) + "ms");
        fileLog.write("Program run time: " + (endTime - startTime) + "ms\n");
        fileLog.close();
        return true;
    }
}
