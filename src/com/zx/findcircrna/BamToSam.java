package com.zx.findcircrna;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
public class BamToSam {
	public void bamToBam(String BamFile,String SamFile) throws IOException {
        // Create SamReaderFactory for reading BAM file
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
        SamReader samReader = samReaderFactory.open(new File(BamFile));

        // Create SAMFileWriter for writing SAM file
        SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        SAMFileWriter samFileWriter = samFileWriterFactory.makeSAMWriter(samReader.getFileHeader(),
                false, new File(SamFile));

        // Convert BAM to SAM and extract only the first ten columns
        for (SAMRecord samRecord : samReader) {
            // Create a new SAMRecord with only the first ten columns
            SAMRecord trimmedRecord = createTrimmedSAMRecord(samRecord);

            // Write the trimmed SAMRecord to the output file
            samFileWriter.addAlignment(trimmedRecord);
        }

        // Close SamReader and SAMFileWriter
        samReader.close();
        samFileWriter.close();
	}
    // Create a new SAMRecord with only the first ten columns
    private static SAMRecord createTrimmedSAMRecord(SAMRecord samRecord) {
        SAMRecord trimmedRecord = new SAMRecord(samRecord.getHeader());
        trimmedRecord.setReadName(samRecord.getReadName());
        trimmedRecord.setFlags(samRecord.getFlags());
        trimmedRecord.setReferenceIndex(samRecord.getReferenceIndex());
        trimmedRecord.setAlignmentStart(samRecord.getAlignmentStart());
        trimmedRecord.setMappingQuality(samRecord.getMappingQuality());
        trimmedRecord.setCigar(samRecord.getCigar());
        trimmedRecord.setMateReferenceIndex(samRecord.getMateReferenceIndex());
        trimmedRecord.setMateAlignmentStart(samRecord.getMateAlignmentStart());
        trimmedRecord.setInferredInsertSize(samRecord.getInferredInsertSize());
        trimmedRecord.setReadString(samRecord.getReadString());

        return trimmedRecord;
    }
}
