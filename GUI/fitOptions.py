# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fitOptions.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_FormFunctions(object):
    def setupUi(self, FormFunctions):
        FormFunctions.setObjectName(_fromUtf8("FormFunctions"))
        FormFunctions.resize(423, 234)
        FormFunctions.setMinimumSize(QtCore.QSize(423, 234))
        FormFunctions.setMaximumSize(QtCore.QSize(423, 234))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Bitstream Vera Sans"))
        FormFunctions.setFont(font)
        self.labelScipyFunc = QtGui.QLabel(FormFunctions)
        self.labelScipyFunc.setGeometry(QtCore.QRect(40, 59, 61, 41))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Bitstream Vera Sans"))
        self.labelScipyFunc.setFont(font)
        self.labelScipyFunc.setAlignment(QtCore.Qt.AlignCenter)
        self.labelScipyFunc.setObjectName(_fromUtf8("labelScipyFunc"))
        self.pushButtonChooseFunc = QtGui.QPushButton(FormFunctions)
        self.pushButtonChooseFunc.setGeometry(QtCore.QRect(20, 100, 110, 51))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Bitstream Vera Sans"))
        self.pushButtonChooseFunc.setFont(font)
        self.pushButtonChooseFunc.setObjectName(_fromUtf8("pushButtonChooseFunc"))
        self.labelFunction = QtGui.QLabel(FormFunctions)
        self.labelFunction.setGeometry(QtCore.QRect(240, 40, 141, 51))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Bitstream Vera Sans"))
        self.labelFunction.setFont(font)
        self.labelFunction.setAlignment(QtCore.Qt.AlignCenter)
        self.labelFunction.setObjectName(_fromUtf8("labelFunction"))
        self.lineEditConstants = QtGui.QLineEdit(FormFunctions)
        self.lineEditConstants.setGeometry(QtCore.QRect(290, 200, 113, 21))
        self.lineEditConstants.setText(_fromUtf8(""))
        self.lineEditConstants.setObjectName(_fromUtf8("lineEditConstants"))
        self.labelConstants = QtGui.QLabel(FormFunctions)
        self.labelConstants.setGeometry(QtCore.QRect(220, 200, 61, 16))
        self.labelConstants.setObjectName(_fromUtf8("labelConstants"))
        self.lineEditFunc = QtGui.QLineEdit(FormFunctions)
        self.lineEditFunc.setGeometry(QtCore.QRect(220, 110, 191, 51))
        self.lineEditFunc.setObjectName(_fromUtf8("lineEditFunc"))
        self.fitpushButton = QtGui.QPushButton(FormFunctions)
        self.fitpushButton.setGeometry(QtCore.QRect(20, 190, 110, 32))
        self.fitpushButton.setObjectName(_fromUtf8("fitpushButton"))

        self.retranslateUi(FormFunctions)
        QtCore.QMetaObject.connectSlotsByName(FormFunctions)

    def retranslateUi(self, FormFunctions):
        FormFunctions.setWindowTitle(_translate("FormFunctions", "Form", None))
        self.labelScipyFunc.setText(_translate("FormFunctions", "Scipy\n"
"Functions", None))
        self.pushButtonChooseFunc.setText(_translate("FormFunctions", "Choose\n"
"Function", None))
        self.labelFunction.setText(_translate("FormFunctions", "Write the Fit function\n"
"in Python format", None))
        self.lineEditConstants.setPlaceholderText(_translate("FormFunctions", "a, b, c", None))
        self.labelConstants.setText(_translate("FormFunctions", "Constants", None))
        self.lineEditFunc.setPlaceholderText(_translate("FormFunctions", "a*np.exp((b*x**2)/c**2)", None))
        self.fitpushButton.setText(_translate("FormFunctions", "Fit", None))

