__author__ = ['titrian','aydin']


from PyQt4 import QtCore, QtGui

from GUI.mainUI import Ui_MainWindow as Dlg

from tools.gibbs_plot import plot2phase, plot3phase, plotYoungsModulus

from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

from GUI.setConcUI import setConcentration_Form
from GUI.setLabelUI import setLabel_Form
from GUI.showTable import showTable_Form

import tools.utilities as u

import sys

class MainWindowGui(QtGui.QMainWindow, Dlg):

    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.crystalList = []
        self.concentrationList = []
        self.labelList = []
        self.kList = []
        self.myks = []
        self.valList = []
        self.youngsList = []
        self.bulkList = []
        self.shearList = []

        self.label_C13.setText("C44")
        self.groupBox_C22.setVisible(False)
        self.groupBox_C23.setVisible(False)
        self.groupBox_C33.setVisible(False)
        self.groupBox_C44.setVisible(False)
        self.groupBox_C55.setVisible(False)
        self.groupBox_C66.setVisible(False)

        self.lineEdit_C11.setText("106.75")
        self.lineEdit_C12.setText("60.41")
        self.lineEdit_C13.setText("28.34")

        self.connect(self.pushButton_Cancel, QtCore.SIGNAL("clicked(bool)"), self.cancelAction)
        self.connect(self.pushButton_Run, QtCore.SIGNAL("clicked(bool)"), self.runAction)
        self.connect(self.pushButton_Add, QtCore.SIGNAL("clicked(bool)"), self.addAction)
        self.connect(self.comboBox, QtCore.SIGNAL("currentIndexChanged(int)"), self.indexChanged)
        self.connect(self.pushButton_Remove, QtCore.SIGNAL("clicked(bool)"), self.removeAction)

#        self.connect(self.pushButton_plot2phase, QtCore.SIGNAL("clicked(bool)"), self.plot2phaseAction)
#        self.connect(self.pushButton_plot3phase, QtCore.SIGNAL("clicked(bool)"), self.plot3phaseAction)

        self.connect(self.comboBox_plot2phase, QtCore.SIGNAL("currentIndexChanged(int)"), self.plot2phaseAction)
        self.connect(self.comboBox_plot3phase, QtCore.SIGNAL("currentIndexChanged(int)"), self.plot3phaseAction)

        self.connect(self.showTable, QtCore.SIGNAL("clicked(bool)"), self.plotTable)
        self.connect(self.pushButton_Clear, QtCore.SIGNAL("clicked(bool)"), self.clearListWidget)
        self.connect(self.pushButton_DeleteAll, QtCore.SIGNAL("clicked(bool)"), self.deleteCrystalList)
        self.connect(self.pushButton, QtCore.SIGNAL("clicked(bool)"), self.runRangeAction)
        self.connect(self.pushButton_Youngs, QtCore.SIGNAL("clicked(bool)"), self.plotYoungsMod)

        self.connect(self.pushButton_CONC, QtCore.SIGNAL("clicked(bool)"), self.setConcentration)

        self.connect(self.listWidget_2, QtCore.SIGNAL("itemClicked(QListWidgetItem*)"), self.setConcentration)

        self.connect(self.listWidget_3, QtCore.SIGNAL("itemClicked(QListWidgetItem*)"), self.setLabel)

        #self.connect(self.pushButton_Clear, QtCore.SIGNAL("clicked(bool)"), self.readOutput)

    def plotYoungsMod(self):
        vals = {}
        type = self.comboBox.currentText()
        if not self.groupBox_C11.isHidden(): vals.update( {str(self.label_C11.text()) : float(self.lineEdit_C11.text())} )
        if not self.groupBox_C12.isHidden(): vals.update( {str(self.label_C12.text()) : float(self.lineEdit_C12.text())} )
        if not self.groupBox_C13.isHidden(): vals.update( {str(self.label_C13.text()) : float(self.lineEdit_C13.text())} )
        if not self.groupBox_C22.isHidden(): vals.update( {str(self.label_C22.text()) : float(self.lineEdit_C22.text())} )
        if not self.groupBox_C23.isHidden(): vals.update( {str(self.label_C23.text()) : float(self.lineEdit_C23.text())} )
        if not self.groupBox_C33.isHidden(): vals.update( {str(self.label_C33.text()) : float(self.lineEdit_C33.text())} )
        if not self.groupBox_C44.isHidden(): vals.update( {str(self.label_C44.text()) : float(self.lineEdit_C44.text())} )
        if not self.groupBox_C55.isHidden(): vals.update( {str(self.label_C55.text()) : float(self.lineEdit_C55.text())} )
        if not self.groupBox_C66.isHidden(): vals.update( {str(self.label_C66.text()) : float(self.lineEdit_C66.text())} )
        plotYoungsModulus(type, vals)
        print(vals)


    def setLabel(self):
        val = self.listWidget_3.currentItem()
        txt = val.text()
        row = self.listWidget_3.row(val)

        self.qdialog = setLabel_Form(txt)
        self.qdialog.exec_()

        self.labelList[row] = self.qdialog.label

        self.refreshListWidget()

    def setConcentration(self):
        self.qdialog = setConcentration_Form(self.crystalList, self.concentrationList)
#        self.qdialog.setFields(self.crystalList, self.concentrationList)
        self.qdialog.exec_()
        self.concentrationList = self.qdialog.newConcList
        self.refreshListWidget()

    def deleteCrystalList(self):
        self.crystalList[:] = []
        self.concentrationList[:] = []
        self.labelList[:] = []
        self.refreshListWidget()

    def clearListWidget(self):
        self.textBrowser.clear()

    def printTxt(self, txt):
        self.textBrowser.insertPlainText(str(txt) + str("\n"))

    def plotTable(self):
        # TODO: check if bulk, young, .. has been calculated
        # TODO: it's hard coded for 2-phase polycrystal. This must be coded for single-phase, 3-phase and 4-phase
        if len(self.youngsList) == 0:
            u.inputError("No data found. Check your calculation: (Range) Run")
        elif len(self.labelList) == 4:
            u.inputError("Table plot for 4-phase polycrystal not implemented!\n"
                         "Implemented for 2-phase or 3-phase polycrystal")
        else:
            self.tableHeaderList = []
            if len(self.labelList) >= 2:
                header0 = "concentration: " + self.labelList[0]
                header1 = " concentration: " + self.labelList[1]
                self.tableHeaderList.append(header0)
                self.tableHeaderList.append(header1)
                if len(self.labelList) == 3:
                    header2 = "concentration: " + self.labelList[2]
                    self.tableHeaderList.append(header2)

                header3 = " Youngs modulus (GPa)"
                header4 = " Bulk modulus (GPa)"
                header5 = " Shear modulus (GPa)"
                self.tableHeaderList.append(header3)
                self.tableHeaderList.append(header4)
                self.tableHeaderList.append(header5)

                concvals_1 = []
                concvals_2 = []
                concvals_3 = []
                for k in range(len(self.kList)):
                    concvals_1.append(self.kList[k][0])
                    concvals_2.append(self.kList[k][1])
                    if len(self.labelList) == 3:
                        concvals_3.append(self.kList[k][2])
                if len(self.labelList) == 2:
                    self.qdialog = showTable_Form(self.tableHeaderList, [concvals_1, concvals_2, self.youngsList, self.bulkList, self.shearList] )
                elif len(self.labelList) == 3:
                    self.qdialog = showTable_Form(self.tableHeaderList, [concvals_1, concvals_2, concvals_3, self.youngsList, self.bulkList, self.shearList] )
            self.qdialog.show()

#    def plot2phaseAction(self):
#        #TODO: move calculation to crystalStack.py
#        if len(self.listWidget) != 2:
#            u.inputError("Exact two phases needed.. nothing to do")
#        elif len(self.youngsList) == 0:
#            u.inputError( "Nothing calculated. Please run calculation first: (Range) Run")
#        else:
#            plot2phase(self.myks, self.youngsList, self.phaseA, self.phaseB, 'Youngs Modulus (GPa)')

    def plot2phaseAction(self):
        plot_name = self.comboBox_plot2phase.currentText()

        if plot_name == "Bulk Modulus":
            if len(self.listWidget) != 2:
                u.inputError("Exact two phases needed.. nothing to do")
            elif len(self.shearList) == 0:
                u.inputError( "Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot2phase(self.myks, self.bulkList, self.phaseA, self.phaseB, 'Bulk Modulus (GPa)')

        elif plot_name == "Shear Modulus":
            if len(self.listWidget) != 2:
                u.inputError("Exact two phases needed.. nothing to do")
            elif len(self.shearList) == 0:
                u.inputError( "Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot2phase(self.myks, self.shearList, self.phaseA, self.phaseB, 'Shear Modulus (GPa)')

        elif plot_name == "Youngs Modulus":
            if len(self.listWidget) != 2:
                u.inputError("Exact two phases needed.. nothing to do")
            elif len(self.shearList) == 0:
                u.inputError( "Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot2phase(self.myks, self.youngsList, self.phaseA, self.phaseB, 'Youngs Modulus (GPa)')

    def plot3phaseAction(self):
        plot_name = self.comboBox_plot3phase.currentText()

        if plot_name == "Bulk Modulus":
            if len(self.listWidget) != 3:
                u.inputError("Exact 3 phases needed.. nothing to do")
            elif len(self.youngsList) == 0:
                u.inputError("Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot3phase(self.myks, self.bulkList, self.phaseA, self.phaseB, self.phaseC,'Bulk Modulus (GPa)',type = "bulkMod")

        elif plot_name == "Shear Modulus":
            if len(self.listWidget) != 3:
                u.inputError("Exact 3 phases needed.. nothing to do")
            elif len(self.youngsList) == 0:
                u.inputError("Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot3phase(self.myks, self.shearList, self.phaseA, self.phaseB, self.phaseC, 'Shear Modulus (GPa)',type = "shearMod")

        if plot_name == "Youngs Modulus":
            if len(self.listWidget) != 3:
                u.inputError("Exact 3 phases needed.. nothing to do")
            elif len(self.youngsList) == 0:
                u.inputError("Nothing calculated. Please run calculation first: (Range) Run")
            else:
                plot3phase(self.myks, self.youngsList, self.phaseA, self.phaseB, self.phaseC, 'Youngs Modulus (GPa)',type = "youngsMod")

#    def plot3phaseAction(self):
#        #TODO: import 3phase plot here, move calculation to crystal Stack
#        if len(self.listWidget) != 3:
#            u.inputError("Exact 3 phases needed.. nothing to do")
#        elif len(self.youngsList) == 0:
#            u.inputError("Nothing calculated. Please run calculation first: (Range) Run")
#        else:
#            plot3phase(self.myks, self.youngsList, self.phaseA, self.phaseB, self.phaseC)

    def indexChanged(self):
        crys_name = self.comboBox.currentText()
        print("box changed to", crys_name,"")

        if crys_name == "cubic":
            print("name is cubic")

            self.label_C11.setText("C11")
            self.label_C12.setText("C12")
            self.label_C13.setText("C13")
            self.label_C22.setText("C22")
            self.label_C23.setText("C23")
            self.label_C33.setText("C33")
            self.label_C44.setText("C44")
            self.label_C55.setText("C55")
            self.label_C66.setText("C66")

            self.groupBox_C11.setVisible(True)
            self.groupBox_C12.setVisible(True)
            self.groupBox_C13.setVisible(True)
            self.groupBox_C22.setVisible(True)
            self.groupBox_C23.setVisible(True)
            self.groupBox_C33.setVisible(True)
            self.groupBox_C44.setVisible(True)
            self.groupBox_C55.setVisible(True)
            self.groupBox_C66.setVisible(True)

            self.label_C13.setText("C44")

            self.groupBox_C22.setVisible(False)
            self.groupBox_C23.setVisible(False)
            self.groupBox_C33.setVisible(False)
            self.groupBox_C44.setVisible(False)
            self.groupBox_C55.setVisible(False)
            self.groupBox_C66.setVisible(False)

            self.lineEdit_C11.setText("106.75")
            self.lineEdit_C12.setText("60.41")
            self.lineEdit_C13.setText("28.34")

        elif crys_name == "tetragonal":
            print("name is tetragonal")

            self.label_C11.setText("C11")
            self.label_C12.setText("C12")
            self.label_C13.setText("C13")
            self.label_C22.setText("C22")
            self.label_C23.setText("C23")
            self.label_C33.setText("C33")
            self.label_C44.setText("C44")
            self.label_C55.setText("C55")
            self.label_C66.setText("C66")

            self.groupBox_C11.setVisible(True)
            self.groupBox_C12.setVisible(True)
            self.groupBox_C13.setVisible(True)
            self.groupBox_C22.setVisible(True)
            self.groupBox_C23.setVisible(True)
            self.groupBox_C33.setVisible(True)
            self.groupBox_C44.setVisible(True)
            self.groupBox_C55.setVisible(True)
            self.groupBox_C66.setVisible(True)

            self.label_C22.setText("C33")
            self.label_C23.setText("C44")
            self.label_C33.setText("C66")

            self.groupBox_C44.setVisible(False)
            self.groupBox_C55.setVisible(False)
            self.groupBox_C66.setVisible(False)

            self.lineEdit_C11.setText("275.12")
            self.lineEdit_C12.setText("178.97")
            self.lineEdit_C13.setText("151.56")
            self.lineEdit_C22.setText("164.86")
            self.lineEdit_C23.setText("54.35")
            self.lineEdit_C33.setText("113.12")

        elif crys_name == "trigonal":
            print("name is trigonal")

            self.label_C11.setText("C11")
            self.label_C12.setText("C12")
            self.label_C13.setText("C13")
            self.label_C22.setText("C22")
            self.label_C23.setText("C23")
            self.label_C33.setText("C33")
            self.label_C44.setText("C44")
            self.label_C55.setText("C55")
            self.label_C66.setText("C66")

            self.groupBox_C11.setVisible(True)
            self.groupBox_C12.setVisible(True)
            self.groupBox_C13.setVisible(True)
            self.groupBox_C22.setVisible(True)
            self.groupBox_C23.setVisible(True)
            self.groupBox_C33.setVisible(True)
            self.groupBox_C44.setVisible(True)
            self.groupBox_C55.setVisible(True)
            self.groupBox_C66.setVisible(True)

            self.label_C22.setText("C14")
            self.label_C23.setText("C33")
            self.label_C33.setText("C44")

            self.groupBox_C44.setVisible(False)
            self.groupBox_C55.setVisible(False)
            self.groupBox_C66.setVisible(False)

            self.lineEdit_C11.setText("101.30")
            self.lineEdit_C12.setText("34.50")
            self.lineEdit_C13.setText("29.20")
            self.lineEdit_C22.setText("20.90")
            self.lineEdit_C23.setText("45.")
            self.lineEdit_C33.setText("39.30")

        elif crys_name == "hexagonal":
            print("name is hexagonal")

            self.label_C11.setText("C11")
            self.label_C12.setText("C12")
            self.label_C13.setText("C13")
            self.label_C22.setText("C22")
            self.label_C23.setText("C23")
            self.label_C33.setText("C33")
            self.label_C44.setText("C44")
            self.label_C55.setText("C55")
            self.label_C66.setText("C66")

            self.groupBox_C11.setVisible(True)
            self.groupBox_C12.setVisible(True)
            self.groupBox_C13.setVisible(True)
            self.groupBox_C22.setVisible(True)
            self.groupBox_C23.setVisible(True)
            self.groupBox_C33.setVisible(True)
            self.groupBox_C44.setVisible(True)
            self.groupBox_C55.setVisible(True)
            self.groupBox_C66.setVisible(True)

            self.label_C22.setText("C33")
            self.label_C23.setText("C55")

            self.groupBox_C33.setVisible(False)
            self.groupBox_C44.setVisible(False)
            self.groupBox_C55.setVisible(False)
            self.groupBox_C66.setVisible(False)

            self.lineEdit_C11.setText("166.7")
            self.lineEdit_C12.setText("13.1")
            self.lineEdit_C13.setText("66.3")
            self.lineEdit_C22.setText("65.5")
            self.lineEdit_C23.setText("139.6")

        elif  crys_name == "orthorombic":
            print("name is orthorombic")

            self.label_C11.setText("C11")
            self.label_C12.setText("C12")
            self.label_C13.setText("C13")
            self.label_C22.setText("C22")
            self.label_C23.setText("C23")
            self.label_C33.setText("C33")
            self.label_C44.setText("C44")
            self.label_C55.setText("C55")
            self.label_C66.setText("C66")

            self.groupBox_C11.setVisible(True)
            self.groupBox_C12.setVisible(True)
            self.groupBox_C13.setVisible(True)
            self.groupBox_C22.setVisible(True)
            self.groupBox_C23.setVisible(True)
            self.groupBox_C33.setVisible(True)
            self.groupBox_C44.setVisible(True)
            self.groupBox_C55.setVisible(True)
            self.groupBox_C66.setVisible(True)

            self.lineEdit_C11.setText("71.6")
            self.lineEdit_C12.setText("26.1")
            self.lineEdit_C13.setText("29.7")
            self.lineEdit_C22.setText("63.2")
            self.lineEdit_C23.setText("29.7")
            self.lineEdit_C33.setText("137.8")
            self.lineEdit_C44.setText("19.6")
            self.lineEdit_C55.setText("24.8")
            self.lineEdit_C66.setText("42.3")

    def sumOfConc(self):
        sum = 0.0
        for k in self.concentrationList:
            sum += float(k)
        return sum

    def calculateConcentration(self):
        if len(self.concentrationList) == 0:
            return 1.0
        else:
            return (1.0 - self.sumOfConc())

    def refreshListWidget(self):
        self.listWidget.clear()
        self.listWidget_2.clear()
        self.listWidget_3.clear()

        for k in self.crystalList:
            self.listWidget.addItem(k)

        for k in self.concentrationList:
            self.listWidget_2.addItem(k)

        for k in self.labelList:
            self.listWidget_3.addItem(k)

    def removeAction(self):
        if self.listWidget.currentItem(): # if item is selected
            val = self.listWidget.currentItem()
            line = self.listWidget.row(val)

            self.crystalList.remove(val.text())

            val_2 = self.listWidget_2.item(line)
            self.concentrationList.remove(val_2.text())

            val_3 = self.listWidget_3.item(line)
            self.labelList.remove(val_3.text())

            self.refreshListWidget()

    def cancelAction(self):
        sys.exit()


    def runRangeAction(self):
        self.kList[:] = []
        self.valList[:] = []
        self.youngsList[:] = []
        self.bulkList[:] = []
        self.shearList[:] = []
        self.myks[:] = []
        if len(self.crystalList) == 0:
            u.inputError("crystal list empty, nothing to do")
        elif self.sumOfConc() != 1:
            u.inputError("Check the concentration!\nSum must be 1")
            self.setConcentration()
        else:
            if len(self.crystalList) == 1:# if just one phase is set, do simple run!
                crys_1 = eval(self.crystalList[0])
                crys_1.setConc(float(self.concentrationList[0]))
                self.crystal = crys_1
            elif len(self.crystalList) == 2:# if 2 phase is set, do range run
                self.phaseA = self.labelList[0]
                self.phaseB = self.labelList[1]
                myvals = []
                myks = []
                crys_1 =  eval(str(self.listWidget.item(0).text()))
                crys_2 =  eval(str(self.listWidget.item(1).text()))
                steps = 10

                for k in range(steps+1):
                    conc = 1 - k*1 / steps    # Thanks to Python 3, we do not need to explicitly state "float"
                    self.kList.append([conc, 1-conc])
                    #print conc
                    self.printTxt( "concentration phase A: %s" % (conc) )

                    crys_1.setConc(1 - k*1 / steps)    # Thanks to Python 3, we do not need to explicitly state "float"
                    crys_2.setConc(k*1 / steps)    # Thanks to Python 3, we do not need to explicitly state "float"
                    mynew = crys_1 + crys_2
#                    if k == 0:
#                        mynew = crys_1
#                        self.pt1Dbulk = (mynew.voigt_bulk)
#                        self.pt1Ubulk = (mynew.reuss_bulk)
#                        self.pt1Dshear = mynew.voigt_shear
#                        self.pt1Ushear = mynew.reuss_shear
#                    if k == 10:
#                        mynew = crys_2
#                        self.pt2Dbulk = (mynew.voigt_bulk)
#                        self.pt2Ubulk = (mynew.reuss_bulk)
#                        self.pt2Dshear = mynew.voigt_shear
#                        self.pt2Ushear = mynew.reuss_shear
                    mynew.setK0andMue0()
                    mynew.setYoungsMod()
                    myvals.append(mynew.E * 100.) # TODO: it's not redundant, same as youngsList, remove one of them
                    self.myks.append(k*1 / steps)    # Thanks to Python 3, we do not need to explicitly state "float"
                    self.progressBar.setValue(k*100 // steps)
                    self.youngsList.append(mynew.E * 100.)
                    self.bulkList.append(mynew.K0 * 100.)
                    self.shearList.append(mynew.mue0 * 100.)
                    if k == 0:
                        mynew = crys_1
                        self.pt1Dbulk = float(mynew.reuss_bulk) * 100.
                        self.pt1Ubulk = float(mynew.voigt_bulk) * 100.
                        self.pt1Dshear = float(mynew.reuss_shear) * 100.
                        self.pt1Ushear = float(mynew.voigt_shear) * 100.
                    if k == 10:
                        mynew = crys_2
                        self.pt2Dbulk = float(mynew.reuss_bulk) * 100.
                        self.pt2Ubulk = float(mynew.voigt_bulk) * 100.
                        self.pt2Dshear = float(mynew.reuss_shear) * 100.
                        self.pt2Ushear = float(mynew.voigt_shear) * 100.
                self.valList = myvals

            elif len(self.crystalList) == 3:
                self.phaseA = self.labelList[0]
                self.phaseB = self.labelList[1]
                self.phaseC = self.labelList[2]

                steps = 10
                crys_1 =  eval(str(self.listWidget.item(0).text()))
                crys_2 =  eval(str(self.listWidget.item(1).text()))
                crys_3 =  eval(str(self.listWidget.item(2).text()))
                myks = []
                myvals = []

                for k in range(steps+1):
                    conc1 = np.round(1 - k*1./steps, 8)
                    for l in range(k+1):
                        conc2 = np.round(1 - conc1 - l*1./steps, 8)
                        conc3 = np.round(1 - conc1 - conc2, 8)
                        self.kList.append([conc1, conc2, conc3])
                        self.printTxt("concentrations of elements A,B,C: %s %s %s" % (conc1, conc2, conc3))

                        myks.append([conc1, conc2, conc3])
                        crys_1.setConc( conc1 )
                        crys_2.setConc( conc2 )
                        crys_3.setConc( conc3 )
                        mynew = crys_1 + crys_2 + crys_3
                        mynew.setK0andMue0()
                        mynew.setYoungsMod()
                        val = float(mynew.E * 100.)
                        self.valList.append(val)
                        self.youngsList.append(float(mynew.E * 100.))
                        self.bulkList.append(float(mynew.K0 * 100.))
                        self.shearList.append(float(mynew.mue0 * 100.))

                        self.progressBar.setValue(k*100 // steps)
                self.myks = myks

    def runAction(self):
        # check concentration = 1
        if len(self.crystalList) == 0:
            u.inputError("crystal lists empty, nothing to do")
        elif np.round(self.sumOfConc(),1) != 1.0:
            u.inputError("Check the concentration!\nSum must be 1.0, it is %s" % (self.sumOfConc()))
            self.setConcentration()
        else:
            if len(self.crystalList) == 1:
                crys_1 = eval(self.crystalList[0])
                crys_1.setConc(float(self.concentrationList[0]))
                self.crystal = crys_1
            elif len(self.crystalList) == 2:
                crys_1 = eval(self.crystalList[0])
                crys_1.setConc(float(self.concentrationList[0]))
                crys_2 = eval(self.crystalList[1])
                crys_2.setConc(float(self.concentrationList[1]))
                self.crystal = crys_1 + crys_2
            elif len(self.crystalList) == 3:
                crys_1 = eval(self.crystalList[0])
                crys_1.setConc(float(self.concentrationList[0]))
                crys_2 = eval(self.crystalList[1])
                crys_2.setConc(float(self.concentrationList[1]))
                crys_3 = eval(self.crystalList[2])
                crys_3.setConc(float(self.concentrationList[2]))
                self.crystal = crys_1 + crys_2 + crys_3
            elif len(self.crystalList) == 4:
                crys_1 = eval(self.crystalList[0])
                crys_1.setConc(float(self.concentrationList[0]))
                crys_2 = eval(self.crystalList[1])
                crys_2.setConc(float(self.concentrationList[1]))
                crys_3 = eval(self.crystalList[2])
                crys_3.setConc(float(self.concentrationList[2]))
                crys_4 = eval(self.crystalList[3])
                crys_4.setConc(float(self.concentrationList[3]))
                self.crystal = crys_1 + crys_2 + crys_3 + crys_4

            self.crystal.setK0andMue0()
            self.crystal.setYoungsMod()
            self.crystal.setpoissonratio()
            self.crystal.areamoduli()

            self.printTxt("Calculated values: Bulk modulus (K0), Shear modulus (mue), Young's modulus (E), Poission ratio (v)\n")
            self.printTxt("Bulk modulus = " + str('%.2f' %(self.crystal.K0 * 100.))+ " (GPa)\nShear modulus = " + str('%.2f' %(self.crystal.mue0 * 100.)) +
                          " (GPa)\nYoung's modulus = " + str('%.2f' %(self.crystal.E * 100.)) + " (GPa)\nPoission ratio = " + str('%.4f' %(self.crystal.v)))
#            self.printTxt("Area moduli (A1) = " + str(self.crystal.A1 * 100.) + "\nArea moduli (A2) = " + str(self.crystal.A2 * 100.) + "\nArea moduli (A3) = " + str(self.crystal.A3 * 100))

            if len(self.crystalList) == 1:
                if self.crystalname == 'cubic' or self.crystalname == 'orthorombic' or self.crystalname == 'tetragonal':
                    self.printTxt("Area moduli [001] = " + str('%.2f' %(self.crystal.A1 * 100.)) + "\nArea moduli [110] = "
                                  + str('%.2f' %(self.crystal.A2 * 100.)) + "\nArea moduli [111] = " + str('%.2f' %(self.crystal.A3 * 100)))

                self.printTxt("Reuss Bulk  = "  + str('%.3f' %(self.crystal.reuss_bulk * 100)))
                self.printTxt("Voigt Bulk  = "  + str('%.3f' %(self.crystal.voigt_bulk * 100)))
                self.printTxt("Reuss Shear = "  + str('%.3f' %(self.crystal.reuss_shear * 100)))
                self.printTxt("Voigt Shear = "  + str('%.3f' %(self.crystal.voigt_shear * 100)))


    def addAction(self):
        if self.comboBox.currentText()=="cubic":
            try:# check if input are just numbers
                eval(str(self.lineEdit_C11.text()))
                eval(str(self.lineEdit_C12.text()))
                eval(str(self.lineEdit_C13.text()))
                self.crystalname = "cubic"
            except:
                u.inputError("Just numbers allowed")
                self.crystalname = None

        if self.comboBox.currentText()=="tetragonal":
            try:# check if input are just numbers
                eval(str(self.lineEdit_C11.text()))
                eval(str(self.lineEdit_C12.text()))
                eval(str(self.lineEdit_C13.text()))
                eval(str(self.lineEdit_C22.text()))
                eval(str(self.lineEdit_C23.text()))
                eval(str(self.lineEdit_C33.text()))
                self.crystalname = "tetragonal"
            except:
                u.inputError("Just numbers allowed")
                self.crystalname = None

        if self.comboBox.currentText()=="trigonal":
            try:# check if input are just numbers
                eval(str(self.lineEdit_C11.text()))
                eval(str(self.lineEdit_C12.text()))
                eval(str(self.lineEdit_C13.text()))
                eval(str(self.lineEdit_C22.text()))
                eval(str(self.lineEdit_C23.text()))
                eval(str(self.lineEdit_C33.text()))
                self.crystalname = "trigonal"
            except:
                u.inputError("Just numbers allowed")
                self.crystalname = None

        if self.comboBox.currentText()=="hexagonal":
            try:# check if input are just numbers
                eval(str(self.lineEdit_C11.text()))
                eval(str(self.lineEdit_C12.text()))
                eval(str(self.lineEdit_C13.text()))
                eval(str(self.lineEdit_C22.text()))
                eval(str(self.lineEdit_C23.text()))
                self.crystalname = "hexagonal"
            except:
                u.inputError("Just numbers allowed")
                self.crystalname = None

        if self.comboBox.currentText()=="orthorombic":
            try:# check if input are just numbers
                eval(str(self.lineEdit_C11.text()))
                eval(str(self.lineEdit_C12.text()))
                eval(str(self.lineEdit_C13.text()))
                eval(str(self.lineEdit_C22.text()))
                eval(str(self.lineEdit_C23.text()))
                eval(str(self.lineEdit_C33.text()))
                eval(str(self.lineEdit_C44.text()))
                eval(str(self.lineEdit_C55.text()))
                eval(str(self.lineEdit_C66.text()))
                self.crystalname = "orthorombic"
            except:
                u.inputError("Just numbers allowed")
                self.crystalname = None

        if self.crystalname == "cubic":
            self.crystal = "cubic(%s, %s, %s)" % (float(self.lineEdit_C11.text()),
                                                  float(self.lineEdit_C12.text()),
                                                  float(self.lineEdit_C13.text()))
        elif self.crystalname == "tetragonal":
            self.crystal = "tetragonal(%s, %s, %s, %s, %s, %s)" % (float(self.lineEdit_C11.text()),
                                                                   float(self.lineEdit_C12.text()),
                                                                   float(self.lineEdit_C13.text()),
                                                                   float(self.lineEdit_C22.text()),
                                                                   float(self.lineEdit_C23.text()),
                                                                   float(self.lineEdit_C33.text()))
        elif self.crystalname == "trigonal":
            self.crystal = "trigonal(%s, %s, %s, %s, %s, %s)" % (float(self.lineEdit_C11.text()),
                                                                 float(self.lineEdit_C12.text()),
                                                                 float(self.lineEdit_C13.text()),
                                                                 float(self.lineEdit_C22.text()),
                                                                 float(self.lineEdit_C23.text()),
                                                                 float(self.lineEdit_C33.text()))
        elif self.crystalname == "hexagonal":
            self.crystal = "hexagonal(%s, %s, %s, %s, %s)" %(float(self.lineEdit_C11.text()),
                                                             float(self.lineEdit_C12.text()),
                                                             float(self.lineEdit_C13.text()),
                                                             float(self.lineEdit_C22.text()),
                                                             float(self.lineEdit_C23.text()))
        elif self.crystalname == "orthorombic":
            self.crystal = "orthorombic(%s, %s, %s, %s, %s, %s, %s, %s, %s)" % (float(self.lineEdit_C11.text()),
                                                                                float(self.lineEdit_C12.text()),
                                                                                float(self.lineEdit_C13.text()),
                                                                                float(self.lineEdit_C22.text()),
                                                                                float(self.lineEdit_C23.text()),
                                                                                float(self.lineEdit_C33.text()),
                                                                                float(self.lineEdit_C44.text()),
                                                                                float(self.lineEdit_C55.text()),
                                                                                float(self.lineEdit_C66.text()))
        if self.crystalname != None and len(self.crystalList) < 4:
            self.crystalList.append(str(self.crystal))
            self.concentrationList.append(str(self.calculateConcentration()))

#            if len(self.labelList) == 0 and "Phase A" not in self.labelList:
#                self.labelList.append("Phase A")
#            elif len(self.labelList) == 1 and "Phase B" not in self.labelList:
#                self.labelList.append("Phase B")
#            elif len(self.labelList) == 2 and "Phase C" not in self.labelList:
#                self.labelList.append("Phase C")
#            elif len(self.labelList) == 3 and "Phase D" not in self.labelList:
#                self.labelList.append("Phase D")
#            else:
            self.labelList.append(self.insertFreeListLabel())

            self.refreshListWidget()


    def insertFreeListLabel(self):
        labels = ["Phase A", "Phase B", "Phase C", "Phase D"]
        #compare current self.labelList with labels
        for elem in labels:
            if elem not in self.labelList:
                return elem

#
#            self.listWidget.clear()
#            for ind in self.crystalList:
#                self.listWidget.insertItem(0, ind) #.insertPlainText(QtCore.QString(crystalname + "("+self.lineEdit_1.text() +")"))
##        self.updateListWidget()
#        #self.textBrowser.insertPlainText(QtCore.QString("\n"))


if __name__ == "__main__":
    import sys, os
    app = QtGui.QApplication(sys.argv)

    #MainWindow = QtGui.QMainWindow()
    #ui = MainWindowGui()
    #ui.setupUi(MainWindow)
    mywindow = MainWindowGui()
    command = "wmic desktopmonitor get screenheight, screenwidth"
    screensize = os.system(command)

    #mywindow.setGeometry(400, 400)
    #print mywindow.geometry()
    mywindow.show()
    #MainWindow.show()
    sys.exit(app.exec_())


#sys.stdout = old_stdout
#
#print(mystdout.getvalue())
