package com.zx.MultiThreaded;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

public class FileReader {  
		 
		private FileChannel fileChanne;
	 
	 
		private ByteBuffer byteBuffer;
	 
		private int bufferSize;
		private long start;
		public FileReader(FileChannel fileChannel, int bufferSize,long start) throws IOException {
			this.fileChanne = fileChannel;

			this.bufferSize = bufferSize;
			this.start = start;
			//this.end = end;
			fileChannel.position(start);
			
		}
	 
		public String readline() throws IOException {
	 
			if (byteBuffer == null) {
				byteBuffer = ByteBuffer.allocate(bufferSize);
	 
				int len = fileChanne.read(byteBuffer);
	 
				if (len == -1) {
					return null;
				}
					
	 
				byteBuffer.flip();
			}
	 
			byte[] bb = new byte[bufferSize];
	 
			int i = 0;
	 
			while (true) {
	 
				while (byteBuffer.hasRemaining()) {
	 
					byte b = byteBuffer.get();
	 
					if ('\n' == b || '\r' == b) {
	 
						if (byteBuffer.hasRemaining()) {
							byte n = byteBuffer.get();
	 
							if ('\n' != n) {
								byteBuffer.position(byteBuffer.position() - 1);
							}
	 
						} else {
	 
							byteBuffer.clear();
	 
							int len = fileChanne.read(byteBuffer);
	 
							byteBuffer.flip();
	 
							if (len != -1) {
								byte n = byteBuffer.get();
	 
								if ('\n' != n) {
									byteBuffer.position(byteBuffer.position() - 1);
								}
							}
	 
						}
	 
						return new String(bb, 0, i);
	 
					} else {
	 
						if (i >= bb.length) {
	 
							bb = Arrays.copyOf(bb, bb.length + bufferSize + 1);
						}
	 
						bb[i++] = b;
					}
	 
				}
	 
				byteBuffer.clear();
				int len = fileChanne.read(byteBuffer);
				byteBuffer.flip();
	 
				if (len == -1 && i == 0) {
					return null;
				}
	 
			}
	 
		}
	 
		public void close() throws IOException {
			this.fileChanne.close();
		}
	 
	
}
