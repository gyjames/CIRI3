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
import com.zx.findcircrna.MutFindCircRNAScan1;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class Scan1Test {
    private int minMapqUni, maxCircle, minCircle, linear_range_size_min, strigency, relExp, seqLen = 0, AllFileSplitNum = 10;
    private long matchNum = 0;
    private boolean intronLable, mlable, spLable;
    private String mitochondrion;

    public Scan1Test(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
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

    public static String samFile;
    private static Lock lock = new ReentrantLock();
    private static HashMap<String, String> chrTCGAMap;
    HashMap<String, String> chrExonStartMap = new HashMap<String, String>(),
            chrExonEndMap = new HashMap<String, String>();
    HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>(),
            chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();

    public boolean CIRI3(String inputFile, String outputFile, String annotationFile,
            String faFile, int threads, String UserGivecircRNA) throws IOException {
        long startTime = System.currentTimeMillis();
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String outputFileLog = outputFile + ".log";
        BufferedWriter fileLog = new BufferedWriter(new FileWriter(new File(outputFileLog)));
        System.out.println(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN1 start");
        fileLog.write(df.format(System.currentTimeMillis()) + " :CIRI3 SCAN1 start\n");

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

        ExecutorService poolExe = Executors.newFixedThreadPool(threads);
        final CyclicBarrier threadSub = new CyclicBarrier(threads + 1);
        final CyclicBarrier threadMain = new CyclicBarrier(threads + 1);
        AtomicInteger incr = new AtomicInteger(1);

        for (int i = 0; i < threads; i++) {
            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        threadMain.await();
                        MutFindCircRNAScan1 scan1 = new MutFindCircRNAScan1(minMapqUni, maxCircle, minCircle,
                                linear_range_size_min, intronLable, chrExonStartMap, chrExonEndMap, chrTCGAMap,
                                chrExonStartTranscriptMap, chrExonEndTranscriptMap, mitochondrion, mlable, spLable);
                        while (true) {
                            int threadNum = incr.getAndIncrement();
                            if (threadNum > AllFileSplitNum) {
                                break;
                            } else {
                                scan1.findCircRNAScan1(samFile, AllFileSplitNum, threadNum);
                                System.out.println(df.format(System.currentTimeMillis()) + " :First scan completed " + threadNum);
                                fileLog.write(df.format(System.currentTimeMillis()) + " :First scan completed " + threadNum + "\n");
                            }
                        }
                        int seqLenTem = scan1.getReadLen();
                        long matchNumTem = scan1.getReadNum();
                        lock.lock();
                        if (seqLenTem > seqLen) {
                            seqLen = seqLenTem;
                        }
                        matchNum = matchNum + matchNumTem;
                        lock.unlock();
                        scan1 = null;
                        threadSub.await();
                        threadMain.await(); // wait for final release, then exit
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            };
            poolExe.execute(runnable);
        }

        try {
            // Load annotation
            HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
            HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
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

            // Compute AllFileSplitNum
            File file = new File(samFile);
            long fileSizeGB = file.length() / 1024 / 1024 / 1024;
            if (fileSizeGB > 200 && fileSizeGB > threads * 10) {
                AllFileSplitNum = (int) fileSizeGB / 10;
            } else {
                AllFileSplitNum = threads;
            }

            // Release threads for scan1
            threadMain.await();
            threadMain.reset();
            threadSub.await(); // wait for scan1 done

            System.out.println(df.format(System.currentTimeMillis()) + " :Mapped_Reads " + matchNum);
            fileLog.write(df.format(System.currentTimeMillis()) + " :Mapped_Reads " + matchNum + "\n");

            // Write scan1_meta
            String metaPath = outputFile + ".scan1_meta";
            BufferedWriter metaBw = new BufferedWriter(new FileWriter(new File(metaPath)));
            metaBw.write("samFile=" + samFile + "\n");
            metaBw.write("readLen=" + seqLen + "\n");
            metaBw.write("readNum=" + matchNum + "\n");
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
        // Do NOT delete BSJ files — they are needed for BUILD_UNIVERSE and SCAN2
        // Do NOT delete converted SAM — it is needed for SCAN2

        long endTime = System.currentTimeMillis();
        System.out.println("Program run time: " + (endTime - startTime) + "ms");
        fileLog.write("Program run time: " + (endTime - startTime) + "ms\n");
        fileLog.close();
        return true;
    }
}
