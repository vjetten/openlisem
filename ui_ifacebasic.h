/********************************************************************************
** Form generated from reading UI file 'ifacebasic.ui'
**
** Created: Thu Mar 4 22:07:01 2010
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IFACEBASIC_H
#define UI_IFACEBASIC_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QToolButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ifacebasicClass
{
public:
    QPushButton *runButton;
    QCheckBox *checkChannel;
    QCheckBox *checkErosion;
    QWidget *layoutWidget;
    QGridLayout *gridLayout_4;
    QLabel *label_8;
    QLabel *label_29;
    QLabel *label_9;
    QLabel *label_30;
    QWidget *layoutWidget1;
    QHBoxLayout *horizontalLayout_2;
    QGridLayout *gridLayout;
    QLabel *label_14;
    QLabel *label_11;
    QLabel *label_2;
    QLabel *label_15;
    QLabel *label_3;
    QLabel *label_10;
    QLabel *label_4;
    QLabel *label_16;
    QLabel *label_5;
    QLabel *label_7;
    QLabel *label_12;
    QLabel *label_6;
    QLabel *label;
    QGridLayout *gridLayout_2;
    QLabel *label_17;
    QLabel *label_13;
    QLabel *label_19;
    QLabel *label_18;
    QLabel *label_21;
    QLabel *label_20;
    QLabel *label_23;
    QLabel *label_22;
    QLabel *label_25;
    QLabel *label_24;
    QLabel *label_27;
    QLabel *label_26;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_31;
    QLineEdit *E_runfilename;
    QToolButton *toolButton_runfilename;
    QProgressBar *progressBar;
    QLabel *label_28;

    void setupUi(QWidget *ifacebasicClass)
    {
        if (ifacebasicClass->objectName().isEmpty())
            ifacebasicClass->setObjectName(QString::fromUtf8("ifacebasicClass"));
        ifacebasicClass->resize(625, 272);
        runButton = new QPushButton(ifacebasicClass);
        runButton->setObjectName(QString::fromUtf8("runButton"));
        runButton->setGeometry(QRect(20, 10, 61, 23));
        runButton->setFlat(false);
        checkChannel = new QCheckBox(ifacebasicClass);
        checkChannel->setObjectName(QString::fromUtf8("checkChannel"));
        checkChannel->setGeometry(QRect(20, 60, 70, 17));
        checkErosion = new QCheckBox(ifacebasicClass);
        checkErosion->setObjectName(QString::fromUtf8("checkErosion"));
        checkErosion->setGeometry(QRect(20, 80, 70, 17));
        layoutWidget = new QWidget(ifacebasicClass);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(120, 130, 161, 41));
        gridLayout_4 = new QGridLayout(layoutWidget);
        gridLayout_4->setSpacing(6);
        gridLayout_4->setContentsMargins(11, 11, 11, 11);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_4->setContentsMargins(0, 0, 0, 0);
        label_8 = new QLabel(layoutWidget);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_4->addWidget(label_8, 0, 0, 1, 1);

        label_29 = new QLabel(layoutWidget);
        label_29->setObjectName(QString::fromUtf8("label_29"));

        gridLayout_4->addWidget(label_29, 0, 1, 1, 1);

        label_9 = new QLabel(layoutWidget);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout_4->addWidget(label_9, 1, 0, 1, 1);

        label_30 = new QLabel(layoutWidget);
        label_30->setObjectName(QString::fromUtf8("label_30"));

        gridLayout_4->addWidget(label_30, 1, 1, 1, 1);

        layoutWidget1 = new QWidget(ifacebasicClass);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(118, 10, 491, 112));
        horizontalLayout_2 = new QHBoxLayout(layoutWidget1);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label_14 = new QLabel(layoutWidget1);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        gridLayout->addWidget(label_14, 0, 0, 1, 1);

        label_11 = new QLabel(layoutWidget1);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        gridLayout->addWidget(label_11, 1, 0, 1, 1);

        label_2 = new QLabel(layoutWidget1);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 1, 1, 2);

        label_15 = new QLabel(layoutWidget1);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        gridLayout->addWidget(label_15, 2, 0, 1, 1);

        label_3 = new QLabel(layoutWidget1);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 2, 1, 1, 2);

        label_10 = new QLabel(layoutWidget1);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout->addWidget(label_10, 3, 0, 1, 1);

        label_4 = new QLabel(layoutWidget1);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 3, 1, 1, 2);

        label_16 = new QLabel(layoutWidget1);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        gridLayout->addWidget(label_16, 4, 0, 1, 1);

        label_5 = new QLabel(layoutWidget1);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout->addWidget(label_5, 4, 1, 1, 1);

        label_7 = new QLabel(layoutWidget1);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout->addWidget(label_7, 4, 2, 1, 1);

        label_12 = new QLabel(layoutWidget1);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        gridLayout->addWidget(label_12, 5, 0, 1, 1);

        label_6 = new QLabel(layoutWidget1);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout->addWidget(label_6, 5, 1, 1, 2);

        label = new QLabel(layoutWidget1);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 1, 1, 2);


        horizontalLayout_2->addLayout(gridLayout);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setSpacing(6);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        label_17 = new QLabel(layoutWidget1);
        label_17->setObjectName(QString::fromUtf8("label_17"));

        gridLayout_2->addWidget(label_17, 0, 0, 1, 1);

        label_13 = new QLabel(layoutWidget1);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        gridLayout_2->addWidget(label_13, 0, 1, 1, 1);

        label_19 = new QLabel(layoutWidget1);
        label_19->setObjectName(QString::fromUtf8("label_19"));

        gridLayout_2->addWidget(label_19, 1, 0, 1, 1);

        label_18 = new QLabel(layoutWidget1);
        label_18->setObjectName(QString::fromUtf8("label_18"));

        gridLayout_2->addWidget(label_18, 1, 1, 1, 1);

        label_21 = new QLabel(layoutWidget1);
        label_21->setObjectName(QString::fromUtf8("label_21"));

        gridLayout_2->addWidget(label_21, 2, 0, 1, 1);

        label_20 = new QLabel(layoutWidget1);
        label_20->setObjectName(QString::fromUtf8("label_20"));

        gridLayout_2->addWidget(label_20, 2, 1, 1, 1);

        label_23 = new QLabel(layoutWidget1);
        label_23->setObjectName(QString::fromUtf8("label_23"));

        gridLayout_2->addWidget(label_23, 3, 0, 1, 1);

        label_22 = new QLabel(layoutWidget1);
        label_22->setObjectName(QString::fromUtf8("label_22"));

        gridLayout_2->addWidget(label_22, 3, 1, 1, 1);

        label_25 = new QLabel(layoutWidget1);
        label_25->setObjectName(QString::fromUtf8("label_25"));

        gridLayout_2->addWidget(label_25, 4, 0, 1, 1);

        label_24 = new QLabel(layoutWidget1);
        label_24->setObjectName(QString::fromUtf8("label_24"));

        gridLayout_2->addWidget(label_24, 4, 1, 1, 1);

        label_27 = new QLabel(layoutWidget1);
        label_27->setObjectName(QString::fromUtf8("label_27"));

        gridLayout_2->addWidget(label_27, 5, 0, 1, 1);

        label_26 = new QLabel(layoutWidget1);
        label_26->setObjectName(QString::fromUtf8("label_26"));

        gridLayout_2->addWidget(label_26, 5, 1, 1, 1);


        horizontalLayout_2->addLayout(gridLayout_2);

        widget = new QWidget(ifacebasicClass);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(10, 180, 589, 24));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_31 = new QLabel(widget);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        horizontalLayout->addWidget(label_31);

        E_runfilename = new QLineEdit(widget);
        E_runfilename->setObjectName(QString::fromUtf8("E_runfilename"));

        horizontalLayout->addWidget(E_runfilename);

        toolButton_runfilename = new QToolButton(widget);
        toolButton_runfilename->setObjectName(QString::fromUtf8("toolButton_runfilename"));
        toolButton_runfilename->setEnabled(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resources/ofolderopen.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton_runfilename->setIcon(icon);

        horizontalLayout->addWidget(toolButton_runfilename);

        progressBar = new QProgressBar(ifacebasicClass);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 240, 601, 23));
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);
        label_28 = new QLabel(ifacebasicClass);
        label_28->setObjectName(QString::fromUtf8("label_28"));
        label_28->setGeometry(QRect(10, 210, 46, 16));

        retranslateUi(ifacebasicClass);

        QMetaObject::connectSlotsByName(ifacebasicClass);
    } // setupUi

    void retranslateUi(QWidget *ifacebasicClass)
    {
        ifacebasicClass->setWindowTitle(QApplication::translate("ifacebasicClass", "ifacebasic", 0, QApplication::UnicodeUTF8));
        runButton->setText(QApplication::translate("ifacebasicClass", "Run", 0, QApplication::UnicodeUTF8));
        checkChannel->setText(QApplication::translate("ifacebasicClass", "Channel", 0, QApplication::UnicodeUTF8));
        checkErosion->setText(QApplication::translate("ifacebasicClass", "Erosion", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("ifacebasicClass", "TextLabel", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("ifacebasicClass", "(min)", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("ifacebasicClass", "TextLabel", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("ifacebasicClass", "sec model runtime", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("ifacebasicClass", "MB (%)", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("ifacebasicClass", "Rain tot", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("ifacebasicClass", "Water Vol", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("ifacebasicClass", "Qtot", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("ifacebasicClass", "Infil (kin wave)", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("ifacebasicClass", "Interception", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("ifacebasicClass", "MBs (%)", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("ifacebasicClass", "MBs (%)", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("ifacebasicClass", "Splash tot", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_21->setText(QApplication::translate("ifacebasicClass", "Flow tot", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_23->setText(QApplication::translate("ifacebasicClass", "Soil loss", 0, QApplication::UnicodeUTF8));
        label_22->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_25->setText(QApplication::translate("ifacebasicClass", "Sed vol", 0, QApplication::UnicodeUTF8));
        label_24->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("ifacebasicClass", "Deposition", 0, QApplication::UnicodeUTF8));
        label_26->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("ifacebasicClass", "Runfile: ", 0, QApplication::UnicodeUTF8));
        toolButton_runfilename->setText(QApplication::translate("ifacebasicClass", "...", 0, QApplication::UnicodeUTF8));
        progressBar->setFormat(QApplication::translate("ifacebasicClass", "%p% (%v-%m)", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("ifacebasicClass", "TextLabel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ifacebasicClass: public Ui_ifacebasicClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IFACEBASIC_H
