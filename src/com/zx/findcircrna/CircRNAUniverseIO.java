package com.zx.findcircrna;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

public class CircRNAUniverseIO {

    /**
     * Writes the universe file.
     * universeDataMap: circFSJMap-key ("chr\tstart\tend") -> siteInfor string ("start\tend[\textra...]")
     * File format: first line "seqLen=N", then one line per circRNA: "chr\tstart\tend[\textra...]"
     */
    public static void writeUniverse(String universePath, int seqLen,
            HashMap<String, String> universeDataMap) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(universePath)));
        bw.write("seqLen=" + seqLen + "\n");
        for (String circKey : universeDataMap.keySet()) {
            // circKey = "chr\tstart\tend"; value = siteInfor = "start\tend[\textra...]"
            String chrKey = circKey.split("\t")[0];
            bw.write(chrKey + "\t" + universeDataMap.get(circKey) + "\n");
        }
        bw.close();
    }

    /**
     * Reads the universe file.
     * Returns map of circFSJMap-key ("chr\tstart\tend") -> arr[] where arr[0]=start, arr[1]=end, arr[2+]=extra.
     * seqLenOut[0] is set to the seqLen from the header line.
     */
    public static HashMap<String, String[]> readUniverse(String universePath,
            int[] seqLenOut) throws IOException {
        HashMap<String, String[]> result = new HashMap<String, String[]>();
        BufferedReader br = new BufferedReader(new FileReader(new File(universePath)));
        String line = br.readLine();
        // first line: seqLen=N
        seqLenOut[0] = Integer.parseInt(line.split("=")[1].trim());
        line = br.readLine();
        while (line != null) {
            if (line.equals("")) {
                line = br.readLine();
                continue;
            }
            String[] parts = line.split("\t");
            // parts[0]=chr, parts[1]=start, parts[2]=end, parts[3+]=extra
            String circKey = parts[0] + "\t" + parts[1] + "\t" + parts[2];
            // arr mirrors siteInfor.split("\t"): arr[0]=start, arr[1]=end, arr[2+]=extra
            String[] arr = Arrays.copyOfRange(parts, 1, parts.length);
            result.put(circKey, arr);
            line = br.readLine();
        }
        br.close();
        return result;
    }

    /**
     * Writes per-sample FSJ counts.
     * circFSJMap: "chr\tstart\tend" -> count
     * File format: one line per circRNA: "chr\tstart\tend\tcount"
     */
    public static void writeFSJCounts(String fsjCountsPath,
            HashMap<String, Integer> circFSJMap) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fsjCountsPath)));
        for (String circKey : circFSJMap.keySet()) {
            bw.write(circKey + "\t" + circFSJMap.get(circKey) + "\n");
        }
        bw.close();
    }
}
