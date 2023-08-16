package smith;

/**
 *
 * @author DulanDias  (AS2014347)
 */
public abstract class SequenceAlignment extends DynamicProgramming {

   protected int match;
   protected int mismatch;
   protected int space;
   protected String[] alignments;

   /*public SequenceAlignment(String sequence1, String sequence2) {
      this(sequence1, sequence2, 1, -1, -1);
   }*/

   public SequenceAlignment( int match,int mismatch, int gap) {
      super();

      this.match = match;
      this.mismatch = mismatch;
      this.space = gap;
   }

   protected Object getTraceback() {
      StringBuffer align1Buf = new StringBuffer();
      StringBuffer align2Buf = new StringBuffer();
      Cell currentCell = getTracebackStartingCell();
      while (traceBackIsNotDone(currentCell)) {
         if (currentCell.getRow() - currentCell.getPrevCell().getRow() == 1) {
            align2Buf.insert(0, sequence2.charAt(currentCell.getRow() - 1));
         } else {
            align2Buf.insert(0, '-');
         }
         if (currentCell.getCol() - currentCell.getPrevCell().getCol() == 1) {
            align1Buf.insert(0, sequence1.charAt(currentCell.getCol() - 1));
         } else {
            align1Buf.insert(0, '-');
         }
         currentCell = currentCell.getPrevCell();
      }

      String[] alignments = new String[] { align1Buf.toString(),
            align2Buf.toString() };

      return alignments;
   }

   protected abstract boolean traceBackIsNotDone(Cell currentCell);

   public int getAlignmentScore() {
      if (alignments == null) {
         getAlignment();
      }

      int score = 0;
      for (int i = 0; i < alignments[0].length(); i++) {
         char c1 = alignments[0].charAt(i);
         char c2 = alignments[1].charAt(i);
         if (c1 == '-' || c2 == '-') {
            score += space;
         } else if (c1 == c2) {
            score += match;
         } else {
            score += mismatch;
         }
      }

      return score;
   }

   public String[] getAlignment() {
      ensureTableIsFilledIn();
      alignments = (String[]) getTraceback();
      return alignments;
   }

   protected abstract Cell getTracebackStartingCell();
}
