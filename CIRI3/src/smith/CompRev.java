package smith;

public class CompRev {

	public String compRev(String seq) {
		StringBuffer newSeq = new StringBuffer(seq);
		String ReSeq = newSeq.reverse().toString();
		String NewReSeq = ReSeq.replace("A", "D")
				.replace("T", "A")
				.replace("G", "E")
				.replace("C", "G")
		        .replace("D", "T")
		        .replace("E", "C");
		return NewReSeq;
		
	}
}
