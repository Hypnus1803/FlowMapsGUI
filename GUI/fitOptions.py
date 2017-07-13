# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fitOptions.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_FormFunctions(object):
    def setupUi(self, FormFunctions):
        FormFunctions.setObjectName("FormFunctions")
        FormFunctions.resize(423, 234)
        FormFunctions.setMinimumSize(QtCore.QSize(423, 234))
        FormFunctions.setMaximumSize(QtCore.QSize(423, 234))
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans")
        FormFunctions.setFont(font)
        self.labelScipyFunc = QtWidgets.QLabel(FormFunctions)
        self.labelScipyFunc.setGeometry(QtCore.QRect(40, 59, 61, 41))
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans")
        self.labelScipyFunc.setFont(font)
        self.labelScipyFunc.setAlignment(QtCore.Qt.AlignCenter)
        self.labelScipyFunc.setObjectName("labelScipyFunc")
        self.pushButtonChooseFunc = QtWidgets.QPushButton(FormFunctions)
        self.pushButtonChooseFunc.setGeometry(QtCore.QRect(20, 100, 110, 51))
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans")
        self.pushButtonChooseFunc.setFont(font)
        self.pushButtonChooseFunc.setObjectName("pushButtonChooseFunc")
        self.labelFunction = QtWidgets.QLabel(FormFunctions)
        self.labelFunction.setGeometry(QtCore.QRect(240, 40, 141, 51))
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans")
        self.labelFunction.setFont(font)
        self.labelFunction.setAlignment(QtCore.Qt.AlignCenter)
        self.labelFunction.setObjectName("labelFunction")
        self.lineEditConstants = QtWidgets.QLineEdit(FormFunctions)
        self.lineEditConstants.setGeometry(QtCore.QRect(290, 200, 113, 21))
        self.lineEditConstants.setText("")
        self.lineEditConstants.setObjectName("lineEditConstants")
        self.labelConstants = QtWidgets.QLabel(FormFunctions)
        self.labelConstants.setGeometry(QtCore.QRect(220, 200, 61, 16))
        self.labelConstants.setObjectName("labelConstants")
        self.lineEditFunc = QtWidgets.QLineEdit(FormFunctions)
        self.lineEditFunc.setGeometry(QtCore.QRect(220, 110, 191, 51))
        self.lineEditFunc.setObjectName("lineEditFunc")
        self.fitpushButton = QtWidgets.QPushButton(FormFunctions)
        self.fitpushButton.setGeometry(QtCore.QRect(20, 190, 110, 32))
        self.fitpushButton.setObjectName("fitpushButton")

        self.retranslateUi(FormFunctions)
        QtCore.QMetaObject.connectSlotsByName(FormFunctions)

    def retranslateUi(self, FormFunctions):
        _translate = QtCore.QCoreApplication.translate
        FormFunctions.setWindowTitle(_translate("FormFunctions", "Form"))
        self.labelScipyFunc.setText(_translate("FormFunctions", "Scipy\n"
"Functions"))
        self.pushButtonChooseFunc.setText(_translate("FormFunctions", "Choose\n"
"Function"))
        self.labelFunction.setText(_translate("FormFunctions", "Write the Fit function\n"
"in Python format"))
        self.lineEditConstants.setPlaceholderText(_translate("FormFunctions", "a, b, c"))
        self.labelConstants.setText(_translate("FormFunctions", "Constants"))
        self.lineEditFunc.setPlaceholderText(_translate("FormFunctions", "a*np.exp((b*x**2)/c**2)"))
        self.fitpushButton.setText(_translate("FormFunctions", "Fit"))

