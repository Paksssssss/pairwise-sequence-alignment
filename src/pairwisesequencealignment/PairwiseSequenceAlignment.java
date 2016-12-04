/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pairwisesequencealignment;

/**
 *
 * @author paks
 */
public class PairwiseSequenceAlignment {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        PairwiseSequenceAligner psa = new PairwiseSequenceAligner(true,false);
        psa.matchScore = 3;
        psa.gapScore = -1;
        psa.mismatchScore = 0;
        
        Sequence seq1 = new Sequence("seq1", "ACTGACTAGA", false);
        Sequence seq2 = new Sequence("seq2", "ACTAAGATA", false);
        psa.seq1 = seq1;
        psa.seq2 = seq2;
        psa.initializeMatrix();
        psa.fillMatrix();
        psa.printMatrix();
    }
    
}
