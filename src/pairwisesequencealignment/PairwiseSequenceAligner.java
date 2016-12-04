/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pairwisesequencealignment;

import java.util.ArrayList;

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
    
    public PairwiseSequenceAligner(boolean isGlobal, boolean isProtein){
        this.global = isGlobal;
        this.isProtein = isProtein;
    }
    
    public boolean scoringSchemeChecker(ArrayList<String> scheme){
        for (String s:scheme) {
            if (s.matches("[^-\\d]")) {
                return false;
            }
        }
        return true;
    }
    
    public void fillMatrix(){
        int top, left, diag;
        for (int j = 1; j < seq2.sequence.length(); j++) {
            for (int i = 1; i < seq1.sequence.length(); i++) {
                if (seq1.sequence.charAt(i) == seq2.sequence.charAt(j)) {
                    diag = matchScore + matrix[i-1][j-1].value;
                } else {
                    diag = mismatchScore + matrix[i-1][j-1].value;
                }
                left = gapScore + matrix[i-1][j].value;
                top =  gapScore + matrix[i][j-1].value;
                int choice = getMax(diag,top,left);
                switch (choice) {
                    case 0:
                        matrix[i][j] = new Cell(diag,false,false,true);
                        break;
                    case 1:
                        matrix[i][j] = new Cell(top,false,true,false);
                        break;
                    case 2:
                        matrix[i][j] = new Cell(left,true,false,false);
                        break;
                    default:
                        break;
                }
            }
        }
    }
    
    private int getMax(int diag, int top, int left){
        if (diag>=top&& diag>=left) {
            return 0;
        } else if (top >= diag&& top>=left) {
            return 1;
        } else {
            return 2;
        }
    }
    
    public void initializeMatrix(){
        seq1.sequence = "-"+ seq1.sequence;
        seq2.sequence = "-"+ seq2.sequence;
        matrix = new Cell[seq1.sequence.length()][seq2.sequence.length()];
        System.out.println(matrix[0][0]);
        if (global) {// im setting initial values of global matrix
            for (int i = 0 , j=0; i < seq1.sequence.length(); i++, j+=-4) {
                matrix[i][0] = new Cell(j);
            }
            for (int i = 0 , j=0; i < seq2.sequence.length(); i++, j+=-4) {
                matrix[0][i] = new Cell(j);
            }
        } else {
            for (int i = 0 ; i < seq1.sequence.length(); i++) {
                matrix[i][0] = new Cell(0);
            }
            for (int i = 0 ; i < seq2.sequence.length(); i++) {
                matrix[0][i] = new Cell(0);
            }
        }
    }
    
    public void printMatrix(){
        for (int i = 0; i < seq1.sequence.length(); i++) {
            for (int j = 0; j < seq2.sequence.length(); j++) {
                System.out.print(matrix[i][j].value + " ");
            }
            System.out.println("");
        }
    }
}

class Sequence {
    String name;
    String sequence;
    boolean isProtorNuc;
    
    public Sequence(){
        this.name = "";
        this.sequence = "";
    }
    
    public Sequence(String name, String sequence, boolean isProtorNuc){
        this.name = name;
        this.sequence = sequence;
        this.isProtorNuc = isProtorNuc;
    }
}

class Cell{
    int value;
    boolean left;
    boolean up;
    boolean diag;
    public Cell(){
        int value = 0;
    }
    public Cell (int value){
        this.value = value;
    }
    public Cell(int value, boolean left, boolean up, boolean diag){
        this.value = value;
        this.left = left;
        this.up = up;
        this.diag= diag;
    }
}
