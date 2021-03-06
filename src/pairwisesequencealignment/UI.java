/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pairwisesequencealignment;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 *
 * @author paks
 */
public class UI extends javax.swing.JFrame {

    /**
     * Creates new form UI
     */
    public UI() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        bindingGroup = new org.jdesktop.beansbinding.BindingGroup();

        uploadFileChooser = new javax.swing.JFileChooser();
        outputFrame = new javax.swing.JFrame();
        jPanel1 = new javax.swing.JPanel();
        jScrollPane2 = new javax.swing.JScrollPane();
        outputTextArea = new javax.swing.JTextArea();
        saveOutput = new javax.swing.JButton();
        jScrollPane1 = new javax.swing.JScrollPane();
        userInput = new javax.swing.JTextArea();
        inputLabel = new javax.swing.JLabel();
        resetButton = new javax.swing.JButton();
        uploadFileButton = new javax.swing.JButton();
        submitButton = new javax.swing.JButton();
        inputOptionsPanel = new javax.swing.JPanel();
        nuclPanel = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        matchScore = new javax.swing.JTextField();
        mismatchScore = new javax.swing.JTextField();
        gapScore = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        nucleotideRadButton = new javax.swing.JRadioButton();
        glocalButton = new javax.swing.JToggleButton();
        proteinPanel = new javax.swing.JPanel();
        jLabel5 = new javax.swing.JLabel();
        scoringMatrixCBox = new javax.swing.JComboBox<>();
        proteinRadButton = new javax.swing.JRadioButton();

        outputTextArea.setEditable(false);
        outputTextArea.setColumns(20);
        outputTextArea.setFont(new java.awt.Font("Monospaced", 0, 13)); // NOI18N
        outputTextArea.setRows(5);
        jScrollPane2.setViewportView(outputTextArea);

        saveOutput.setText("Save");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane2)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(151, 151, 151)
                .addComponent(saveOutput, javax.swing.GroupLayout.PREFERRED_SIZE, 233, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(153, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 419, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(saveOutput, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                .addContainerGap())
        );

        javax.swing.GroupLayout outputFrameLayout = new javax.swing.GroupLayout(outputFrame.getContentPane());
        outputFrame.getContentPane().setLayout(outputFrameLayout);
        outputFrameLayout.setHorizontalGroup(
            outputFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        outputFrameLayout.setVerticalGroup(
            outputFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Pairwise Sequence Alignment");

        userInput.setColumns(20);
        userInput.setRows(5);
        jScrollPane1.setViewportView(userInput);

        inputLabel.setFont(new java.awt.Font("Ubuntu", 1, 15)); // NOI18N
        inputLabel.setText("Input");

        resetButton.setText("Reset");
        resetButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                resetButtonActionPerformed(evt);
            }
        });

        uploadFileButton.setText("Upload File");
        uploadFileButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                uploadFileButtonActionPerformed(evt);
            }
        });

        submitButton.setText("Submit");
        submitButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                submitButtonActionPerformed(evt);
            }
        });

        jLabel1.setFont(new java.awt.Font("Ubuntu", 1, 17)); // NOI18N
        jLabel1.setText("Scoring Scheme");

        matchScore.setText("1");

        mismatchScore.setText("0");
        mismatchScore.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mismatchScoreActionPerformed(evt);
            }
        });

        gapScore.setText("-1");
        gapScore.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                gapScoreActionPerformed(evt);
            }
        });

        jLabel2.setText("Match       :");

        org.jdesktop.beansbinding.Binding binding = org.jdesktop.beansbinding.Bindings.createAutoBinding(org.jdesktop.beansbinding.AutoBinding.UpdateStrategy.READ_WRITE, matchScore, org.jdesktop.beansbinding.ObjectProperty.create(), jLabel2, org.jdesktop.beansbinding.BeanProperty.create("labelFor"));
        bindingGroup.addBinding(binding);

        jLabel3.setText("Mismatch :");

        binding = org.jdesktop.beansbinding.Bindings.createAutoBinding(org.jdesktop.beansbinding.AutoBinding.UpdateStrategy.READ_WRITE, mismatchScore, org.jdesktop.beansbinding.ObjectProperty.create(), jLabel3, org.jdesktop.beansbinding.BeanProperty.create("labelFor"));
        bindingGroup.addBinding(binding);

        jLabel4.setText("Gap          :");

        binding = org.jdesktop.beansbinding.Bindings.createAutoBinding(org.jdesktop.beansbinding.AutoBinding.UpdateStrategy.READ_WRITE, gapScore, org.jdesktop.beansbinding.ObjectProperty.create(), jLabel4, org.jdesktop.beansbinding.BeanProperty.create("labelFor"));
        bindingGroup.addBinding(binding);

        nucleotideRadButton.setFont(new java.awt.Font("Ubuntu", 0, 17)); // NOI18N
        nucleotideRadButton.setSelected(true);
        nucleotideRadButton.setText("Nucleotide Sequence");
        nucleotideRadButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nucleotideRadButtonActionPerformed(evt);
            }
        });

        glocalButton.setText("Local");
        glocalButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                glocalButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout nuclPanelLayout = new javax.swing.GroupLayout(nuclPanel);
        nuclPanel.setLayout(nuclPanelLayout);
        nuclPanelLayout.setHorizontalGroup(
            nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(nuclPanelLayout.createSequentialGroup()
                .addGap(43, 43, 43)
                .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(nucleotideRadButton)
                    .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                        .addGroup(nuclPanelLayout.createSequentialGroup()
                            .addComponent(jLabel1)
                            .addGap(28, 28, 28))
                        .addGroup(nuclPanelLayout.createSequentialGroup()
                            .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(jLabel4, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(jLabel2, javax.swing.GroupLayout.DEFAULT_SIZE, 76, Short.MAX_VALUE)
                                .addComponent(jLabel3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                            .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(gapScore, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(matchScore, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(mismatchScore, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addGroup(nuclPanelLayout.createSequentialGroup()
                        .addGap(44, 44, 44)
                        .addComponent(glocalButton, javax.swing.GroupLayout.PREFERRED_SIZE, 86, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(70, Short.MAX_VALUE))
        );
        nuclPanelLayout.setVerticalGroup(
            nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(nuclPanelLayout.createSequentialGroup()
                .addComponent(nucleotideRadButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(matchScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(mismatchScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(nuclPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(gapScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(glocalButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jLabel5.setFont(new java.awt.Font("Ubuntu", 1, 17)); // NOI18N
        jLabel5.setText("Scoring Matrix");

        scoringMatrixCBox.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "PAM120 (Global)", "BIOSUM62 (Local)" }));
        scoringMatrixCBox.setEnabled(false);

        proteinRadButton.setFont(new java.awt.Font("Ubuntu", 0, 17)); // NOI18N
        proteinRadButton.setText("Protein Sequence");
        proteinRadButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                proteinRadButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout proteinPanelLayout = new javax.swing.GroupLayout(proteinPanel);
        proteinPanel.setLayout(proteinPanelLayout);
        proteinPanelLayout.setHorizontalGroup(
            proteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinPanelLayout.createSequentialGroup()
                .addGroup(proteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(proteinPanelLayout.createSequentialGroup()
                        .addGap(33, 33, 33)
                        .addComponent(proteinRadButton))
                    .addGroup(proteinPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(proteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(proteinPanelLayout.createSequentialGroup()
                                .addGap(12, 12, 12)
                                .addComponent(scoringMatrixCBox, javax.swing.GroupLayout.PREFERRED_SIZE, 201, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(jLabel5))))
                .addContainerGap(40, Short.MAX_VALUE))
        );
        proteinPanelLayout.setVerticalGroup(
            proteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinPanelLayout.createSequentialGroup()
                .addComponent(proteinRadButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel5, javax.swing.GroupLayout.PREFERRED_SIZE, 21, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(scoringMatrixCBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(86, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout inputOptionsPanelLayout = new javax.swing.GroupLayout(inputOptionsPanel);
        inputOptionsPanel.setLayout(inputOptionsPanelLayout);
        inputOptionsPanelLayout.setHorizontalGroup(
            inputOptionsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(inputOptionsPanelLayout.createSequentialGroup()
                .addContainerGap(29, Short.MAX_VALUE)
                .addComponent(nuclPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(proteinPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(37, 37, 37))
        );
        inputOptionsPanelLayout.setVerticalGroup(
            inputOptionsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(inputOptionsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(inputOptionsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(inputOptionsPanelLayout.createSequentialGroup()
                        .addComponent(proteinPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(nuclPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane1))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(23, 23, 23)
                        .addComponent(inputLabel)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
            .addGroup(layout.createSequentialGroup()
                .addGap(49, 49, 49)
                .addComponent(resetButton, javax.swing.GroupLayout.PREFERRED_SIZE, 119, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(97, 97, 97)
                .addComponent(uploadFileButton, javax.swing.GroupLayout.PREFERRED_SIZE, 145, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(submitButton, javax.swing.GroupLayout.PREFERRED_SIZE, 115, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(38, 38, 38))
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(inputOptionsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(inputLabel)
                .addGap(5, 5, 5)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 199, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(resetButton)
                    .addComponent(uploadFileButton)
                    .addComponent(submitButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(inputOptionsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        bindingGroup.bind();

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void submitButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_submitButtonActionPerformed
        if (!this.proteinRadButton.isSelected() && !this.nucleotideRadButton.isSelected()) {
            JOptionPane.showMessageDialog(this, "Please choose which kind of sequence!", "Error", JOptionPane.ERROR_MESSAGE);
        } else {
            if (!this.userInput.getText().startsWith(">")) {
                JOptionPane.showMessageDialog(this, "Input not in FASTA Format", "Input Error", JOptionPane.ERROR_MESSAGE);
            } else {
                if (proteinRadButton.isSelected()) {
                    psa.isProtein = true;
                    if (scoringMatrixCBox.getSelectedIndex() == 0) {
                        psa.global = true;
                    } else {
                        psa.global = false;
                    }
                    try {
                        psa.setScoringMatrix();
                    } catch (IOException ex) {
                        Logger.getLogger(UI.class.getName()).log(Level.SEVERE, null, ex);
                    }
                } else if (nucleotideRadButton.isSelected()) {
                    psa.isProtein = false;
                    psa.global = glocalButton.isSelected();
                    System.out.println();
                    ArrayList<String> scoreScheme = new ArrayList();
                    scoreScheme.add(matchScore.getText());
                    scoreScheme.add(mismatchScore.getText());
                    scoreScheme.add(gapScore.getText());
                    if (!psa.scoringSchemeChecker(scoreScheme)) {
                        JOptionPane.showMessageDialog(this, "Invalid Scoring Scheme!", "Error", JOptionPane.ERROR_MESSAGE);
                    } else {
                        psa.matchScore = Integer.parseInt(this.matchScore.getText());
                        psa.mismatchScore = Integer.parseInt(this.mismatchScore.getText());
                        psa.gapScore = Integer.parseInt(this.gapScore.getText());
                    }
                }
                ArrayList<String> parseThisInput = new ArrayList(Arrays.asList(this.userInput.getText().split("\n")));
                if (psa.parseInput(parseThisInput)) {
                    if (psa.checkInput()) {
                        psa.initializeMatrix();
                        psa.fillMatrix();
                        psa.solve();
                        generateOutput();

                    } else {
                        JOptionPane.showMessageDialog(this, "Invalid Input!", "Error", JOptionPane.ERROR_MESSAGE);
                    }
                } else {
                    JOptionPane.showMessageDialog(this, "Invalid Input!", "Error", JOptionPane.ERROR_MESSAGE);
                }
            }
        }
    }//GEN-LAST:event_submitButtonActionPerformed

    private void nucleotideRadButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nucleotideRadButtonActionPerformed
        this.proteinRadButton.setSelected(false);
        this.glocalButton.setEnabled(true);
        this.matchScore.setEnabled(true);
        this.mismatchScore.setEnabled(true);
        this.gapScore.setEnabled(true);
        this.scoringMatrixCBox.setEnabled(false);
    }//GEN-LAST:event_nucleotideRadButtonActionPerformed

    private void mismatchScoreActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mismatchScoreActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_mismatchScoreActionPerformed

    private void gapScoreActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_gapScoreActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_gapScoreActionPerformed

    private void proteinRadButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_proteinRadButtonActionPerformed
        this.nucleotideRadButton.setSelected(false);
        this.glocalButton.setEnabled(false);
        this.matchScore.setEnabled(false);
        this.mismatchScore.setEnabled(false);
        this.gapScore.setEnabled(false);
        this.scoringMatrixCBox.setEnabled(true);
    }//GEN-LAST:event_proteinRadButtonActionPerformed

    private void resetButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_resetButtonActionPerformed
        this.userInput.setText("");
        this.nucleotideRadButton.setSelected(true);
        this.proteinRadButton.setSelected(false);
        this.matchScore.setText("1");
        this.mismatchScore.setText("0");
        this.gapScore.setText("-1");
        this.scoringMatrixCBox.setSelectedIndex(0);
        this.glocalButton.setEnabled(true);
        this.matchScore.setEnabled(true);
        this.mismatchScore.setEnabled(true);
        this.gapScore.setEnabled(true);
        this.scoringMatrixCBox.setEnabled(false);
    }//GEN-LAST:event_resetButtonActionPerformed

    private void glocalButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_glocalButtonActionPerformed
        if (glocalButton.isSelected()) {
            glocalButton.setText("Global");
        } else {
            glocalButton.setText("Local");
        }
    }//GEN-LAST:event_glocalButtonActionPerformed

    private void uploadFileButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_uploadFileButtonActionPerformed
        uploadFileChooser.setMultiSelectionEnabled(false);
        uploadFileChooser.setFileFilter(new CustomFileFilter());
        int returnVal = uploadFileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = uploadFileChooser.getSelectedFile();
            try {
                // What to do with the file, e.g. display it in a TextArea
                this.userInput.read(new FileReader(file.getAbsolutePath()), null);
            } catch (FileNotFoundException ex) {
                System.out.println("problem accessing file" + file.getAbsolutePath());
            } catch (IOException ex) {
                System.out.println("problem accessing file" + file.getAbsolutePath());
            }
        } else {
            System.out.println("File access cancelled by user.");
        }
    }//GEN-LAST:event_uploadFileButtonActionPerformed

    public void showResults() {

    }

    public void generateOutput() {
        outputTextArea.setText("");
        Date date = new Date();
        String a;
        char let;
        int occ1, occ2;
        if (psa.isProtein) {
            a = "ARNDCQEGHILKMFPSTWYVBZX";
        } else {
            a = "ACTG";
        }
        String output = "Pairwise Sequence Alignment ver. 1.0 by John Vincent N. Pakson (2013-54677)\n";
        output += "Run date:"+ date.toString()    +"\n\n Submitted Sequences:\n";
        output += ">"+psa.seq1.name +"\n"+ psa.seq1.sequence +"\n\n";
        output += ">"+psa.seq2.name +"\n"+ psa.seq2.sequence +"\n\n";
        output += ">"+psa.seq1.name + " length: "+ (psa.seq1.sequence.length()-1)+"\n";
        output += ">"+psa.seq2.name + " length: "+ (psa.seq2.sequence.length()-1)+"\n";
        output += "Frequence Occurence:\nSQ:\tS1\ts2\tTotal\n";
        for (int i = 0; i < a.length(); i++) {
            let = a.charAt(i);
            occ1 = psa.seq1.sequence.replaceAll("[^"+let+"]", "").length();
            occ2 = psa.seq2.sequence.replaceAll("[^"+let+"]", "").length();
            output += let+":\t"+occ1+"\t"+occ2+"\t"+(occ1+occ2)+"\n";
        }
        output += printAlignment();
        output += "Score: "+ psa.alignmentScore;
    
        this.saveOutput.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                writeToFile();
            }
        });
        this.outputTextArea.setText(output);
        
        this.outputFrame.pack();
        this.outputFrame.setVisible(true);
    }
    
    public String printAlignment(){
        String seq1,seq2,alignment, finString = "";
        int a = 1;
        psa.seq1.sequence = psa.seq1.sequence.substring(1);
        psa.seq2.sequence = psa.seq2.sequence.substring(1);
        for (String al: psa.alignments) {
            seq1 = seq2 = alignment = "";
            finString += "\nAlignment " + a + " \n";
            for (int i = 0,j =0, k=0; i < al.length(); i++) {
                if (i%10==0 && i!=0) {
                    seq1 += " ";
                    alignment += " ";
                    seq2 +=" ";
                } else if (al.charAt(i) == '*' || al.charAt(i) == '.') {
                    seq1 += psa.seq1.sequence.charAt(j);
                    alignment += al.charAt(i);
                    seq2 +=psa.seq2.sequence.charAt(k);
                    j++;
                    k++;
                } else if (al.charAt(i) == '-') {
                    seq2 += '-';
                    alignment += " ";
                    seq1 += psa.seq1.sequence.charAt(j);
                    j++;
                } else if (al.charAt(i) == '^'){
                    seq1 += '-';
                    alignment += " ";
                    seq2 += psa.seq2.sequence.charAt(k);
                    k++;
                }
            }
            a++;
            finString +=seq1+"\n"+alignment+"\n"+seq2+"\n\n";
        }
        return finString;
    }
    
    public void writeToFile(){
        JFileChooser saveFile = new JFileChooser();
        saveFile.setMultiSelectionEnabled(false);
        saveFile.setDialogTitle("Export Report");
        // Demonstrate "Save" dialog:
        int returnVal = saveFile.showSaveDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            try {
                String dir = saveFile.getSelectedFile().getCanonicalPath();
                //Path out = Paths.get(dir);
                FileWriter fw = new FileWriter(dir);
                this.outputTextArea.write(fw);
            } catch (IOException ex) {
                System.out.println("NO SUCH PATH");
            }
        } else {
            System.out.println("File access cancelled by user.");
        }
    }

    static class CustomFileFilter extends javax.swing.filechooser.FileFilter {

        @Override
        public boolean accept(File f) {
            // Allow only directories, or files with ".txt" extension
            return f.isDirectory() || f.getAbsolutePath().endsWith(".FASTA") || f.getAbsolutePath().endsWith(".fasta") || f.getAbsolutePath().endsWith(".txt");
        }

        @Override
        public String getDescription() {
            return "Fasta Files (*.FASTA) or Text Files (*.txt)";
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Windows".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(UI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(UI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(UI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(UI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new UI().setVisible(true);
            }
        });
    }

    public HashMap occurences1;
    public HashMap occurences2;
    public PairwiseSequenceAligner psa = new PairwiseSequenceAligner();
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextField gapScore;
    private javax.swing.JToggleButton glocalButton;
    private javax.swing.JLabel inputLabel;
    private javax.swing.JPanel inputOptionsPanel;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JTextField matchScore;
    private javax.swing.JTextField mismatchScore;
    private javax.swing.JPanel nuclPanel;
    private javax.swing.JRadioButton nucleotideRadButton;
    private javax.swing.JFrame outputFrame;
    private javax.swing.JTextArea outputTextArea;
    private javax.swing.JPanel proteinPanel;
    private javax.swing.JRadioButton proteinRadButton;
    private javax.swing.JButton resetButton;
    private javax.swing.JButton saveOutput;
    private javax.swing.JComboBox<String> scoringMatrixCBox;
    private javax.swing.JButton submitButton;
    private javax.swing.JButton uploadFileButton;
    private javax.swing.JFileChooser uploadFileChooser;
    private javax.swing.JTextArea userInput;
    private org.jdesktop.beansbinding.BindingGroup bindingGroup;
    // End of variables declaration//GEN-END:variables
}
