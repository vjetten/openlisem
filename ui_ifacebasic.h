/********************************************************************************
** Form generated from reading UI file 'ifacebasic.ui'
**
** Created: Thu May 13 13:32:44 2010
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
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_31;
    QLineEdit *E_runfilename;
    QToolButton *toolButton_runfilename;
    QToolButton *toolButton_ShowRunfile;
    QProgressBar *progressBar;
    QLabel *label_28;
    QWidget *layoutWidget1;
    QGridLayout *gridLayout_3;
    QLabel *label_time;
    QLabel *label_29;
    QLabel *label_runtime;
    QLabel *label_30;
    QLabel *label_endruntime;
    QLabel *label_33;
    QLabel *label_endtime;
    QLabel *label_38;
    QLabel *label_dx;
    QLabel *label_46;
    QLabel *label_area;
    QLabel *label_48;
    QWidget *layoutWidget2;
    QGridLayout *gridLayout_2;
    QLabel *label_19;
    QLabel *label_splashdet;
    QLabel *label_flowdet;
    QLabel *label_27;
    QLabel *label_dep;
    QLabel *label_25;
    QLabel *label_sedvol;
    QLabel *label_41;
    QLabel *label_detch;
    QLabel *label_43;
    QLabel *label_depch;
    QLabel *label_39;
    QLabel *label_sedvolch;
    QLabel *label_23;
    QLabel *label_soilloss;
    QLabel *label_17;
    QLabel *label_MBs;
    QLabel *label;
    QLabel *label_soillosskgha;
    QLabel *label_21;
    QWidget *widget;
    QGridLayout *gridLayout;
    QLabel *label_raintot;
    QLabel *label_15;
    QLabel *label_watervoltot;
    QLabel *label_16;
    QLabel *label_infiltot;
    QLabel *label_12;
    QLabel *label_interctot;
    QLabel *label_10;
    QLabel *label_qtot;
    QLabel *label_7;
    QLabel *label_qtotm3;
    QLabel *label_35;
    QLabel *label_qpeak;
    QLabel *label_14;
    QLabel *label_MB;
    QLabel *label_11;

    void setupUi(QWidget *ifacebasicClass)
    {
        if (ifacebasicClass->objectName().isEmpty())
            ifacebasicClass->setObjectName(QString::fromUtf8("ifacebasicClass"));
        ifacebasicClass->resize(601, 409);
        ifacebasicClass->setStyleSheet(QString::fromUtf8("Cleanlooks"));
        runButton = new QPushButton(ifacebasicClass);
        runButton->setObjectName(QString::fromUtf8("runButton"));
        runButton->setGeometry(QRect(20, 10, 61, 23));
        runButton->setFlat(false);
        layoutWidget = new QWidget(ifacebasicClass);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(90, 10, 501, 24));
        horizontalLayout = new QHBoxLayout(layoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_31 = new QLabel(layoutWidget);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        horizontalLayout->addWidget(label_31);

        E_runfilename = new QLineEdit(layoutWidget);
        E_runfilename->setObjectName(QString::fromUtf8("E_runfilename"));

        horizontalLayout->addWidget(E_runfilename);

        toolButton_runfilename = new QToolButton(layoutWidget);
        toolButton_runfilename->setObjectName(QString::fromUtf8("toolButton_runfilename"));
        toolButton_runfilename->setEnabled(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/fileopen.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton_runfilename->setIcon(icon);

        horizontalLayout->addWidget(toolButton_runfilename);

        toolButton_ShowRunfile = new QToolButton(layoutWidget);
        toolButton_ShowRunfile->setObjectName(QString::fromUtf8("toolButton_ShowRunfile"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/new.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton_ShowRunfile->setIcon(icon1);

        horizontalLayout->addWidget(toolButton_ShowRunfile);

        progressBar = new QProgressBar(ifacebasicClass);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 380, 581, 23));
        progressBar->setAutoFillBackground(false);
        progressBar->setMaximum(99);
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);
        label_28 = new QLabel(ifacebasicClass);
        label_28->setObjectName(QString::fromUtf8("label_28"));
        label_28->setGeometry(QRect(10, 330, 581, 41));
        label_28->setFrameShape(QFrame::Panel);
        label_28->setFrameShadow(QFrame::Sunken);
        label_28->setWordWrap(true);
        layoutWidget1 = new QWidget(ifacebasicClass);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(90, 40, 341, 53));
        gridLayout_3 = new QGridLayout(layoutWidget1);
        gridLayout_3->setSpacing(6);
        gridLayout_3->setContentsMargins(11, 11, 11, 11);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setContentsMargins(0, 0, 0, 0);
        label_time = new QLabel(layoutWidget1);
        label_time->setObjectName(QString::fromUtf8("label_time"));
        label_time->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_time, 0, 0, 1, 1);

        label_29 = new QLabel(layoutWidget1);
        label_29->setObjectName(QString::fromUtf8("label_29"));

        gridLayout_3->addWidget(label_29, 0, 1, 1, 1);

        label_runtime = new QLabel(layoutWidget1);
        label_runtime->setObjectName(QString::fromUtf8("label_runtime"));
        label_runtime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_runtime, 1, 0, 1, 1);

        label_30 = new QLabel(layoutWidget1);
        label_30->setObjectName(QString::fromUtf8("label_30"));

        gridLayout_3->addWidget(label_30, 1, 1, 1, 1);

        label_endruntime = new QLabel(layoutWidget1);
        label_endruntime->setObjectName(QString::fromUtf8("label_endruntime"));
        label_endruntime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_endruntime, 1, 2, 1, 1);

        label_33 = new QLabel(layoutWidget1);
        label_33->setObjectName(QString::fromUtf8("label_33"));

        gridLayout_3->addWidget(label_33, 1, 3, 1, 1);

        label_endtime = new QLabel(layoutWidget1);
        label_endtime->setObjectName(QString::fromUtf8("label_endtime"));
        label_endtime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_endtime, 0, 2, 1, 1);

        label_38 = new QLabel(layoutWidget1);
        label_38->setObjectName(QString::fromUtf8("label_38"));

        gridLayout_3->addWidget(label_38, 0, 3, 1, 1);

        label_dx = new QLabel(layoutWidget1);
        label_dx->setObjectName(QString::fromUtf8("label_dx"));
        label_dx->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_dx, 2, 0, 1, 1);

        label_46 = new QLabel(layoutWidget1);
        label_46->setObjectName(QString::fromUtf8("label_46"));

        gridLayout_3->addWidget(label_46, 2, 1, 1, 1);

        label_area = new QLabel(layoutWidget1);
        label_area->setObjectName(QString::fromUtf8("label_area"));
        label_area->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_area, 2, 2, 1, 1);

        label_48 = new QLabel(layoutWidget1);
        label_48->setObjectName(QString::fromUtf8("label_48"));

        gridLayout_3->addWidget(label_48, 2, 3, 1, 1);

        layoutWidget2 = new QWidget(ifacebasicClass);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(340, 100, 251, 221));
        gridLayout_2 = new QGridLayout(layoutWidget2);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        label_19 = new QLabel(layoutWidget2);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_19, 0, 0, 1, 1);

        label_splashdet = new QLabel(layoutWidget2);
        label_splashdet->setObjectName(QString::fromUtf8("label_splashdet"));
        label_splashdet->setFrameShape(QFrame::StyledPanel);
        label_splashdet->setFrameShadow(QFrame::Sunken);
        label_splashdet->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_splashdet, 0, 1, 1, 1);

        label_flowdet = new QLabel(layoutWidget2);
        label_flowdet->setObjectName(QString::fromUtf8("label_flowdet"));
        label_flowdet->setFrameShape(QFrame::StyledPanel);
        label_flowdet->setFrameShadow(QFrame::Sunken);
        label_flowdet->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_flowdet, 1, 1, 1, 1);

        label_27 = new QLabel(layoutWidget2);
        label_27->setObjectName(QString::fromUtf8("label_27"));
        label_27->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_27, 2, 0, 1, 1);

        label_dep = new QLabel(layoutWidget2);
        label_dep->setObjectName(QString::fromUtf8("label_dep"));
        label_dep->setFrameShape(QFrame::StyledPanel);
        label_dep->setFrameShadow(QFrame::Sunken);
        label_dep->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_dep, 2, 1, 1, 1);

        label_25 = new QLabel(layoutWidget2);
        label_25->setObjectName(QString::fromUtf8("label_25"));
        label_25->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_25, 3, 0, 1, 1);

        label_sedvol = new QLabel(layoutWidget2);
        label_sedvol->setObjectName(QString::fromUtf8("label_sedvol"));
        label_sedvol->setFrameShape(QFrame::StyledPanel);
        label_sedvol->setFrameShadow(QFrame::Sunken);
        label_sedvol->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_sedvol, 3, 1, 1, 1);

        label_41 = new QLabel(layoutWidget2);
        label_41->setObjectName(QString::fromUtf8("label_41"));
        label_41->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_41, 4, 0, 1, 1);

        label_detch = new QLabel(layoutWidget2);
        label_detch->setObjectName(QString::fromUtf8("label_detch"));
        label_detch->setFrameShape(QFrame::StyledPanel);
        label_detch->setFrameShadow(QFrame::Sunken);
        label_detch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_detch, 4, 1, 1, 1);

        label_43 = new QLabel(layoutWidget2);
        label_43->setObjectName(QString::fromUtf8("label_43"));
        label_43->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_43, 5, 0, 1, 1);

        label_depch = new QLabel(layoutWidget2);
        label_depch->setObjectName(QString::fromUtf8("label_depch"));
        label_depch->setFrameShape(QFrame::StyledPanel);
        label_depch->setFrameShadow(QFrame::Sunken);
        label_depch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_depch, 5, 1, 1, 1);

        label_39 = new QLabel(layoutWidget2);
        label_39->setObjectName(QString::fromUtf8("label_39"));
        label_39->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_39, 6, 0, 1, 1);

        label_sedvolch = new QLabel(layoutWidget2);
        label_sedvolch->setObjectName(QString::fromUtf8("label_sedvolch"));
        label_sedvolch->setFrameShape(QFrame::StyledPanel);
        label_sedvolch->setFrameShadow(QFrame::Sunken);
        label_sedvolch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_sedvolch, 6, 1, 1, 1);

        label_23 = new QLabel(layoutWidget2);
        label_23->setObjectName(QString::fromUtf8("label_23"));
        label_23->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_23, 7, 0, 1, 1);

        label_soilloss = new QLabel(layoutWidget2);
        label_soilloss->setObjectName(QString::fromUtf8("label_soilloss"));
        label_soilloss->setFrameShape(QFrame::StyledPanel);
        label_soilloss->setFrameShadow(QFrame::Sunken);
        label_soilloss->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_soilloss, 7, 1, 1, 1);

        label_17 = new QLabel(layoutWidget2);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        QPalette palette;
        QBrush brush(QColor(128, 128, 128, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        QBrush brush1(QColor(122, 122, 122, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_17->setPalette(palette);
        label_17->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_17, 9, 0, 1, 1);

        label_MBs = new QLabel(layoutWidget2);
        label_MBs->setObjectName(QString::fromUtf8("label_MBs"));
        QPalette palette1;
        palette1.setBrush(QPalette::Active, QPalette::WindowText, brush);
        palette1.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        palette1.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_MBs->setPalette(palette1);
        label_MBs->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_MBs, 9, 1, 1, 1);

        label = new QLabel(layoutWidget2);
        label->setObjectName(QString::fromUtf8("label"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label, 8, 0, 1, 1);

        label_soillosskgha = new QLabel(layoutWidget2);
        label_soillosskgha->setObjectName(QString::fromUtf8("label_soillosskgha"));
        label_soillosskgha->setFrameShape(QFrame::StyledPanel);
        label_soillosskgha->setFrameShadow(QFrame::Sunken);
        label_soillosskgha->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_soillosskgha, 8, 1, 1, 1);

        label_21 = new QLabel(layoutWidget2);
        label_21->setObjectName(QString::fromUtf8("label_21"));
        label_21->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_21, 1, 0, 1, 1);

        widget = new QWidget(ifacebasicClass);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(90, 100, 221, 181));
        gridLayout = new QGridLayout(widget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_raintot = new QLabel(widget);
        label_raintot->setObjectName(QString::fromUtf8("label_raintot"));
        label_raintot->setAutoFillBackground(false);
        label_raintot->setFrameShape(QFrame::StyledPanel);
        label_raintot->setFrameShadow(QFrame::Sunken);
        label_raintot->setLineWidth(1);
        label_raintot->setMidLineWidth(0);
        label_raintot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_raintot, 0, 1, 1, 1);

        label_15 = new QLabel(widget);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_15, 1, 0, 1, 1);

        label_watervoltot = new QLabel(widget);
        label_watervoltot->setObjectName(QString::fromUtf8("label_watervoltot"));
        label_watervoltot->setFrameShape(QFrame::StyledPanel);
        label_watervoltot->setFrameShadow(QFrame::Sunken);
        label_watervoltot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_watervoltot, 1, 1, 1, 1);

        label_16 = new QLabel(widget);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_16, 2, 0, 1, 1);

        label_infiltot = new QLabel(widget);
        label_infiltot->setObjectName(QString::fromUtf8("label_infiltot"));
        label_infiltot->setFrameShape(QFrame::StyledPanel);
        label_infiltot->setFrameShadow(QFrame::Sunken);
        label_infiltot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_infiltot, 2, 1, 1, 1);

        label_12 = new QLabel(widget);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_12, 3, 0, 1, 1);

        label_interctot = new QLabel(widget);
        label_interctot->setObjectName(QString::fromUtf8("label_interctot"));
        label_interctot->setFrameShape(QFrame::StyledPanel);
        label_interctot->setFrameShadow(QFrame::Sunken);
        label_interctot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_interctot, 3, 1, 1, 1);

        label_10 = new QLabel(widget);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_10, 4, 0, 1, 1);

        label_qtot = new QLabel(widget);
        label_qtot->setObjectName(QString::fromUtf8("label_qtot"));
        label_qtot->setFrameShape(QFrame::StyledPanel);
        label_qtot->setFrameShadow(QFrame::Sunken);
        label_qtot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_qtot, 4, 1, 1, 1);

        label_7 = new QLabel(widget);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_7, 5, 0, 1, 1);

        label_qtotm3 = new QLabel(widget);
        label_qtotm3->setObjectName(QString::fromUtf8("label_qtotm3"));
        label_qtotm3->setFrameShape(QFrame::StyledPanel);
        label_qtotm3->setFrameShadow(QFrame::Sunken);
        label_qtotm3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_qtotm3, 5, 1, 1, 1);

        label_35 = new QLabel(widget);
        label_35->setObjectName(QString::fromUtf8("label_35"));
        label_35->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_35, 6, 0, 1, 1);

        label_qpeak = new QLabel(widget);
        label_qpeak->setObjectName(QString::fromUtf8("label_qpeak"));
        label_qpeak->setFrameShape(QFrame::StyledPanel);
        label_qpeak->setFrameShadow(QFrame::Sunken);
        label_qpeak->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_qpeak, 6, 1, 1, 1);

        label_14 = new QLabel(widget);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        QPalette palette2;
        QBrush brush2(QColor(129, 129, 129, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette2.setBrush(QPalette::Active, QPalette::WindowText, brush2);
        palette2.setBrush(QPalette::Inactive, QPalette::WindowText, brush2);
        palette2.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_14->setPalette(palette2);
        label_14->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_14, 7, 0, 1, 1);

        label_MB = new QLabel(widget);
        label_MB->setObjectName(QString::fromUtf8("label_MB"));
        QPalette palette3;
        palette3.setBrush(QPalette::Active, QPalette::WindowText, brush2);
        palette3.setBrush(QPalette::Inactive, QPalette::WindowText, brush2);
        palette3.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_MB->setPalette(palette3);
        label_MB->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_MB, 7, 1, 1, 1);

        label_11 = new QLabel(widget);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_11, 0, 0, 1, 1);


        retranslateUi(ifacebasicClass);

        QMetaObject::connectSlotsByName(ifacebasicClass);
    } // setupUi

    void retranslateUi(QWidget *ifacebasicClass)
    {
        ifacebasicClass->setWindowTitle(QApplication::translate("ifacebasicClass", "openLisem (1.0 beta)", 0, QApplication::UnicodeUTF8));
        runButton->setText(QApplication::translate("ifacebasicClass", "Run", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("ifacebasicClass", "Runfile: ", 0, QApplication::UnicodeUTF8));
        toolButton_runfilename->setText(QApplication::translate("ifacebasicClass", "...", 0, QApplication::UnicodeUTF8));
        toolButton_ShowRunfile->setText(QApplication::translate("ifacebasicClass", "...", 0, QApplication::UnicodeUTF8));
        progressBar->setFormat(QApplication::translate("ifacebasicClass", "%p% (%v-%m)", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("ifacebasicClass", "!!!", 0, QApplication::UnicodeUTF8));
        label_time->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("ifacebasicClass", "sim time (min)", 0, QApplication::UnicodeUTF8));
        label_runtime->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("ifacebasicClass", "runtime  (min)", 0, QApplication::UnicodeUTF8));
        label_endruntime->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_33->setText(QApplication::translate("ifacebasicClass", "Est. total run time (min)", 0, QApplication::UnicodeUTF8));
        label_endtime->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_38->setText(QApplication::translate("ifacebasicClass", "End time (min)", 0, QApplication::UnicodeUTF8));
        label_dx->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_46->setText(QApplication::translate("ifacebasicClass", "Cellsize (m)", 0, QApplication::UnicodeUTF8));
        label_area->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_48->setText(QApplication::translate("ifacebasicClass", "ha", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("ifacebasicClass", "Splash detachment (ton)", 0, QApplication::UnicodeUTF8));
        label_splashdet->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_flowdet->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("ifacebasicClass", "Deposition (ton)", 0, QApplication::UnicodeUTF8));
        label_dep->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_25->setText(QApplication::translate("ifacebasicClass", "Sed volume in flow (ton)", 0, QApplication::UnicodeUTF8));
        label_sedvol->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_41->setText(QApplication::translate("ifacebasicClass", "Detachment channel (ton)", 0, QApplication::UnicodeUTF8));
        label_detch->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_43->setText(QApplication::translate("ifacebasicClass", "Deposition channel (ton)", 0, QApplication::UnicodeUTF8));
        label_depch->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_39->setText(QApplication::translate("ifacebasicClass", "Sed vol in channel (ton)", 0, QApplication::UnicodeUTF8));
        label_sedvolch->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_23->setText(QApplication::translate("ifacebasicClass", "Soil loss (ton)", 0, QApplication::UnicodeUTF8));
        label_soilloss->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("ifacebasicClass", "MBs (%)", 0, QApplication::UnicodeUTF8));
        label_MBs->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("ifacebasicClass", "Soil loss (kg/ha)", 0, QApplication::UnicodeUTF8));
        label_soillosskgha->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_21->setText(QApplication::translate("ifacebasicClass", "Flow detachment (ton)", 0, QApplication::UnicodeUTF8));
        label_raintot->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("ifacebasicClass", "Water Vol (mm)", 0, QApplication::UnicodeUTF8));
        label_watervoltot->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("ifacebasicClass", "Infil (mm)", 0, QApplication::UnicodeUTF8));
        label_infiltot->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("ifacebasicClass", "Interception (mm)", 0, QApplication::UnicodeUTF8));
        label_interctot->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("ifacebasicClass", "Qtot (mm)", 0, QApplication::UnicodeUTF8));
        label_qtot->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("ifacebasicClass", "Qtot (m3)", 0, QApplication::UnicodeUTF8));
        label_qtotm3->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_35->setText(QApplication::translate("ifacebasicClass", "Qpeak (l/s)", 0, QApplication::UnicodeUTF8));
        label_qpeak->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("ifacebasicClass", "MB (%)", 0, QApplication::UnicodeUTF8));
        label_MB->setText(QApplication::translate("ifacebasicClass", "0", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("ifacebasicClass", "Rain tot (mm)", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ifacebasicClass: public Ui_ifacebasicClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IFACEBASIC_H
