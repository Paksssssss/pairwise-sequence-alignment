/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pairwisesequencealignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;

/**
 *
 * @author paks
 */
public class PairwiseSequenceAligner {

    Cell matrix[][];
    boolean global;
    boolean isProtein;
    int matchScore, mismatchScore, gapScore;
    Sequence seq1, seq2;
    ArrayList<String> alignments;
    ScoringMatrix scoringMatrix;
    int alignmentScore;

    public PairwiseSequenceAligner(boolean isGlobal) {
        this.global = isGlobal;
    }

    PairwiseSequenceAligner() {
        alignments = new ArrayList();
        seq1 = new Sequence();
        seq2 = new Sequence();
    }

    public boolean parseInput(ArrayList<String> userInput) {
        ArrayList<Sequence> sequenceSet = new ArrayList();
        Sequence newSequence = new Sequence();
        for (String s : userInput) {
            if (s.startsWith(">")) {
                if (!newSequence.name.isEmpty()&&!newSequence.sequence.isEmpty()) {
                    sequenceSet.add(newSequence);
                    newSequence = new Sequence();
                }
                newSequence.name = s.substring(1);
            } else if(s.isEmpty()){
                if (!newSequence.name.isEmpty()&&!newSequence.sequence.isEmpty()) {
                    sequenceSet.add(newSequence);
                    newSequence = new Sequence();
                }
                newSequence.name = s;
            }  else if (s.equals(userInput.get(userInput.size()-1))) {
                s=s.toUpperCase();
                newSequence.sequence += s;
                sequenceSet.add(newSequence);
            } else {
                s=s.toUpperCase();
                newSequence.sequence += s;
            }
        }
        if (sequenceSet.size()>2 || sequenceSet.size() == 0) {
            return false;
        } else {
            seq1 = sequenceSet.get(0);
            seq2 = sequenceSet.get(1);
            seq1.sequence = seq1.sequence.toUpperCase();
            seq2.sequence = seq2.sequence.toUpperCase();
            return true;
        }
    }
    
    public boolean checkInput(){
        if (seq1.sequence.isEmpty()||seq2.sequence.isEmpty()){
            return false;
        }
        if (isProtein) {
            if (seq1.sequence.contains("B")||seq1.sequence.contains("U")
                    ||seq1.sequence.contains("O")||seq1.sequence.contains("Z")
                            ||seq1.sequence.contains("X")||seq1.sequence.contains("J")) {
                return false;
            }
            if (seq2.sequence.contains("B")||seq2.sequence.contains("U")
                    ||seq2.sequence.contains("O")||seq2.sequence.contains("Z")
                            ||seq2.sequence.contains("X")||seq2.sequence.contains("J")) {
                return false;
            }
        } else {
            for (int i = 0; i < seq1.sequence.length(); i++) {
                if (seq1.sequence.charAt(i)!='A'&&seq1.sequence.charAt(i)!='T'&&seq1.sequence.charAt(i)!='C'&&seq1.sequence.charAt(i)!='G') {
                    return false;
                }
            } for (int i = 0; i < seq2.sequence.length(); i++) {
                if (seq2.sequence.charAt(i)!='A'&&seq2.sequence.charAt(i)!='T'&&seq2.sequence.charAt(i)!='C'&&seq2.sequence.charAt(i)!='G') {
                    return false;
                }
            }
        }
        return true;
    }

    public void setScoringMatrix() throws IOException {
        if (global) {
            scoringMatrix = new ScoringMatrix("./pam120.txt");
        } else {
            scoringMatrix = new ScoringMatrix("./blosum62.txt");
        }
        gapScore = scoringMatrix.getScore('*', 'A');
    }

    public boolean scoringSchemeChecker(ArrayList<String> scheme) {
        for (String s : scheme) {
            if (s.matches("[^-\\d]")) {
                return false;
            }
        }
        return true;
    }

    public void fillMatrix() {
        int top, left, diag;
        boolean mismatch;
        System.out.println(seq2.sequence.length());
        System.out.println(seq1.sequence.length());
        int score = 0;
        for (int j = 1; j < seq2.sequence.length(); j++) {
            for (int i = 1; i < seq1.sequence.length(); i++) {
                mismatch = false;
                if (isProtein) {
                    if (seq1.sequence.charAt(i) == seq2.sequence.charAt(j)) {
                        diag = scoringMatrix.getScore(seq1.sequence.charAt(i), seq2.sequence.charAt(j)) + matrix[i - 1][j - 1].value;
                    } else {
                        diag = scoringMatrix.getScore(seq1.sequence.charAt(i), seq2.sequence.charAt(j)) + matrix[i - 1][j - 1].value;
                        mismatch = true;
                    }
                    left = gapScore + matrix[i - 1][j].value;
                    top = gapScore + matrix[i][j - 1].value;
                } else {
                    if (seq1.sequence.charAt(i) == seq2.sequence.charAt(j)) {
                        diag = matchScore + matrix[i - 1][j - 1].value;
                    } else {
                        diag = mismatchScore + matrix[i - 1][j - 1].value;
                        mismatch = true;
                    }
                    left = gapScore + matrix[i - 1][j].value;
                    top = gapScore + matrix[i][j - 1].value;
                }
                int choice = getMax(diag, top, left);
                System.out.print(choice+" ");
                if (!global && diag <= 0 && top <= 0 && left <= 0) {
                    matrix[i][j] = new Cell(0, false, false, false, mismatch);
                } else {
                    matrix[i][j] = new Cell(choice, mismatch);
                    matrix[i][j].diag = choice == diag;
                    matrix[i][j].up = choice == top;
                    matrix[i][j].left = choice == left;
                }
            }
            System.out.println("");
        }
    }

    public void solve(){
        if (global) {
            traceback("",seq1.sequence.length()-1,seq2.sequence.length()-1);
            alignmentScore = matrix[seq1.sequence.length()-1][seq2.sequence.length()-1].value;
        } else {
            int max, i = 1, j = 1; 
            max = matrix[0][0].value;
            for (; j < seq2.sequence.length(); j++) {
                for (i = 1; i < seq1.sequence.length(); i++) {
                    if (matrix[i][j].value > max) {
                        max = matrix[i][j].value;
                    }
                }
            }
            System.out.println(max);
            for (j=1; j < seq2.sequence.length(); j++) {
                for (i = 1; i < seq1.sequence.length(); i++) {
                    if (matrix[i][j].value == max) {
                        traceback("",i,j);
                    }
                }
            }
            alignmentScore = max;
        }
        System.out.println(alignments.size());
        alignments = new ArrayList<String>(new LinkedHashSet<String>(alignments));
    }

    public String traceback(String alignment, int i, int j) {
        System.out.println(alignment + "");

        if (matrix[i][j].diag) {
            if (matrix[i][j].mismatch) {
                alignment = "." + alignment;
            } else {
                alignment = "*" + alignment;
            }
            traceback(alignment, --i, --j);
        }
        if (matrix[i][j].left) {
            alignment = "-" + alignment;
            traceback(alignment, --i, j);
        }
        if (matrix[i][j].up) {
            alignment = "^" + alignment;
            traceback(alignment, i, --j);
        }
        if (global && (i == 0 || j == 0)) {
            alignments.add(alignment);
            return alignment;
        } else if (!global && (i == 1 || j ==1)) {
            alignments.add(alignment);
            return alignment;
        }
        return alignment;
    }

    private int getMax(int diag, int top, int left) {
        if (diag >= top && diag >= left) {
            return diag;
        } else if (top >= diag && top >= left) {
            return top;
        } else {
            return left;
        }
    }

    public void initializeMatrix() {
        seq1.sequence = "-" + seq1.sequence;
        seq2.sequence = "-" + seq2.sequence;
        matrix = new Cell[seq1.sequence.length()][seq2.sequence.length()];
        System.out.println(matrix[0][0]);
        if (global) {// im setting initial values of global matrix
            for (int i = 0, j = 0; i < seq1.sequence.length(); i++, j += gapScore) {
                matrix[i][0] = new Cell(j);
            }
            for (int i = 0, j = 0; i < seq2.sequence.length(); i++, j += gapScore) {
                matrix[0][i] = new Cell(j);
            }
        } else {
            for (int i = 0; i < seq1.sequence.length(); i++) {
                matrix[i][0] = new Cell(0);
            }
            for (int i = 0; i < seq2.sequence.length(); i++) {
                matrix[0][i] = new Cell(0);
            }
        }
    }

    public void printMatrix() {
        for (int j = 0; j < seq2.sequence.length(); j++) {
            for (int i = 0; i < seq1.sequence.length(); i++) {
                System.out.print(matrix[i][j].value + "\t");
            }
            System.out.println("");
        }
    }
    
    public HashMap countOccurence(String sequence){
        HashMap tempMap = new HashMap();
        String a;
        char let;
        int occ;
        if (isProtein) {
            a = "ARNDCQEGHILKMFPSTWYVBZX";
        } else {
            a = "ACTG";
        }
        
        for (int i = 0; i < a.length(); i++) {
            let = a.charAt(i);
            occ = a.replaceAll("[^"+let+"]", sequence).length();
        }
        return tempMap;
    }
}

class Sequence {

    String name;
    String sequence;
    boolean isProtorNuc;

    public Sequence() {
        this.name = "";
        this.sequence = "";
    }

    public Sequence(String name, String sequence, boolean isProtorNuc) {
        this.name = name;
        this.sequence = sequence;
        this.isProtorNuc = isProtorNuc;
    }
}

class Cell {

    int value;
    boolean left;
    boolean up;
    boolean diag;
    boolean mismatch;

    public Cell() {
        int value = 0;
    }

    public Cell(int value) {
        this.value = value;
    }

    public Cell(int value, boolean mismatch){
        this.value = value;
        this.mismatch = mismatch;
    }

    public Cell(int value, boolean left, boolean up, boolean diag, boolean mismatch) {
        this.value = value;
        this.left = left;
        this.up = up;
        this.diag = diag;
        this.mismatch = mismatch;
    }
}
