package smith;

/**
 *
 * @author DulanDias  (AS2014347)
 */
public class SmithWaterman extends SequenceAlignment {

   private Cell highScoreCell;

   /*public SmithWaterman(String sequence1, String sequence2) {
      super(sequence1, sequence2);
   }*/

   public SmithWaterman(int match,
         int mismatch, int gap) {
      super(match, mismatch, gap);
   }

   protected void initialize() {
      super.initialize();

      highScoreCell = scoreTable[0][0];
   }

   protected void fillInCell(Cell currentCell, Cell cellAbove, Cell cellToLeft,
         Cell cellAboveLeft) {
      int rowSpaceScore = cellAbove.getScore() + space;
      int colSpaceScore = cellToLeft.getScore() + space;
      int matchOrMismatchScore = cellAboveLeft.getScore();
      if (sequence2.charAt(currentCell.getRow() - 1) == sequence1
            .charAt(currentCell.getCol() - 1)) {
         matchOrMismatchScore += match;
      } else {
         matchOrMismatchScore += mismatch;
      }
      if (rowSpaceScore >= colSpaceScore) {
         if (matchOrMismatchScore >= rowSpaceScore) {
            if (matchOrMismatchScore > 0) {
               currentCell.setScore(matchOrMismatchScore);
               currentCell.setPrevCell(cellAboveLeft);
            }
         } else {
            if (rowSpaceScore > 0) {
               currentCell.setScore(rowSpaceScore);
               currentCell.setPrevCell(cellAbove);
            }
         }
      } else {
         if (matchOrMismatchScore >= colSpaceScore) {
            if (matchOrMismatchScore > 0) {
               currentCell.setScore(matchOrMismatchScore);
               currentCell.setPrevCell(cellAboveLeft);
            }
         } else {
            if (colSpaceScore > 0) {
               currentCell.setScore(colSpaceScore);
               currentCell.setPrevCell(cellToLeft);
            }
         }
      }
      if (currentCell.getScore() > highScoreCell.getScore()) {
         highScoreCell = currentCell;
      }
   }

   /*
    * (non-Javadoc)
    * 
    * @see java.lang.Object#toString()
    */
   @Override
   public String toString() {
      return "[NeedlemanWunsch: sequence1=" + sequence1 + ", sequence2="
            + sequence2 + "]";
   }

   @Override
   protected boolean traceBackIsNotDone(Cell currentCell) {
      return currentCell.getScore() != 0;
   }

   @Override
   protected Cell getTracebackStartingCell() {
      return highScoreCell;
   }

   @Override
   protected Cell getInitialPointer(int row, int col) {
      return null;
   }

   @Override
   protected int getInitialScore(int row, int col) {
      return 0;
   }
}