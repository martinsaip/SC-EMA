__author__ = ['titrian','aydin']

try:
    from PyQt4 import QtGui, QtCore

    def inputError(text):
        abort_txt = "Warning: " + text + "\n"
        qSelectDialog = QtGui.QMessageBox()
        qSelectDialog.setWindowTitle("Error")
        qSelectDialog.setText(abort_txt)
        qSelectDialog.setStandardButtons(qSelectDialog.Ok)

        qSelectDialog.exec_()
except:
    import warnings
    def inputError(text):
        warnings.warn ("Warning: " + text)

