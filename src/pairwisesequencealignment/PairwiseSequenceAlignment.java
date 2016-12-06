/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pairwisesequencealignment;

import java.io.IOException;

/**
 *
 * @author paks
 */
public class PairwiseSequenceAlignment {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        PairwiseSequenceAligner psa = new PairwiseSequenceAligner(true);
        psa.matchScore = 5;
        psa.gapScore = -4;
        psa.mismatchScore = -3;
        
        Sequence seq1 = new Sequence("seq1", "VLSPADKFLTNV", false);
        Sequence seq2 = new Sequence("seq2", "VFTELSPAKTV", false);
        psa.seq1 = seq1;
        psa.seq2 = seq2;
        psa.isProtein = true;
        psa.global = true;
        psa.setScoringMatrix();
        psa.initializeMatrix();
        psa.fillMatrix();
        psa.printMatrix();
        psa.solve();
        for (String s: psa.alignments) {
            System.out.println(s);
        }
    }
    
}
