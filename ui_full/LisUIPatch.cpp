/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include "lisemqt.h"

#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QSslSocket>
#include <QSysInfo>


// NOTE: on windows for a working version, openssl must be installed in msys. Then the file
// "qopensslbackend.dll" must be installed in the subdirectory tls where the exe is

//-------------------------------------------------------------------------------------

void lisemqt::downloadPatch(QString latestVersion)
{
    QString patchName = QString("openLisemSetup_%1.exe").arg(latestVersion);
    QString urltxt = QString("https://github.com/vjetten/openlisem/releases/download/lisem_bin/%1?raw=true").arg(patchName);

    // Ask the user if they want to download and install the patch
    QMessageBox::StandardButton replyButton;
    replyButton = QMessageBox::question(nullptr, "LISEM version check",
                                        QString("A newer LISEM version is available (%1).\nDo you want to download and install the new version?").arg(latestVersion),
                                        QMessageBox::Yes | QMessageBox::No);

    if (replyButton == QMessageBox::Yes) {
        QString downloadsDir = QStandardPaths::writableLocation(QStandardPaths::DownloadLocation);
        QString filePath = QDir(downloadsDir).filePath(patchName);

        // Create a file to save the patch
        QFile *file = new QFile(filePath);
        if (!file->open(QIODevice::WriteOnly)) {
            qDebug() << "Error: Unable to open patch file.";
            delete file;
            return;
        }

        // Show download progress
        QProgressDialog *progressDialog = new QProgressDialog("Downloading patch...", "Cancel", 0, 100, nullptr);
        progressDialog->setWindowModality(Qt::WindowModal);
        progressDialog->setMinimumDuration(0);

        // Proceed with the download
        QNetworkReply *reply = manager->get(QNetworkRequest(QUrl(urltxt)));

        QObject::connect(reply, &QNetworkReply::downloadProgress, [=](qint64 bytesReceived, qint64 bytesTotal) mutable {
            if (bytesTotal > 0) {
                progressDialog->setValue(static_cast<int>((bytesReceived * 100) / bytesTotal));
                if (bytesReceived == bytesTotal) {
                    progressDialog->setLabelText("Finalizing...");
                }
            }
        });

        QObject::connect(progressDialog, &QProgressDialog::canceled, [=]() mutable {
            reply->abort();
            file->close();
            file->remove();
            QMessageBox::information(nullptr, "Download Cancelled", "The download has been cancelled.");
            progressDialog->deleteLater();
        });

        QObject::connect(reply, &QNetworkReply::readyRead, [=]() mutable {
            file->write(reply->readAll());
        });

        QObject::connect(reply, &QNetworkReply::finished, [=]() mutable {
            int err = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
            qDebug() << err;
            // Handle error scenario
            if (err != 200) {
                file->close();
                file->remove();
                QMessageBox::critical(nullptr, "Error downloading", QString("Cannot download file, return code GitHub server = %1").arg(err));
                progressDialog->deleteLater();
                return;
            }

            // Close the new file
            file->close();
            reply->deleteLater();
            file->deleteLater();

            if (reply->error() == QNetworkReply::NoError) {
                progressDialog->setValue(100);
                progressDialog->deleteLater();

                QMessageBox::StandardButton installButton;
                installButton = QMessageBox::question(nullptr, "Install new version",
                                                      "The new version has been downloaded successfully in your Download folder. Do you want to close Lisem and install the new version? (the old version will be uninstalled, your list of runfiles is preserved)",
                                                      QMessageBox::Yes | QMessageBox::No);

                if (installButton == QMessageBox::Yes) {
                    // Execute the downloaded file
                    QProcess::startDetached(filePath);

                    // Close the current application
                    QApplication::quit();
                }
            } else {
                // Handle download failure
                QMessageBox::critical(nullptr, "Download Failed", "Failed to download the new version.");
                qDebug() << "Error:" << reply->errorString();
                progressDialog->deleteLater();
            }
        });

        // Show the progress dialog
        progressDialog->show();
    } else {
        replyButton = QMessageBox::question(nullptr, "New LISEM version",
                                            QString("Do you want to continue checking for new versions?\nYou can activate this again in the Advanced Options."),
                                            QMessageBox::Yes | QMessageBox::No);
        checkforpatch = replyButton == QMessageBox::Yes;
        toolButton_version->setChecked(checkforpatch  );
    }
}

//-------------------------------------------------------------------------------------
bool lisemqt::isNewVersionAvailable(QString &GitHubVersion)
{
    // Assuming version strings are in the format "major.minor.patch"
    QStringList currentParts = QString(VERSIONNR).split(".");//currentVersion.split(".");
    QStringList githubParts = GitHubVersion.split(".");
    QString revision;
    QString revisionGIT;
    bool beta = false;
    bool betaGIT = false;
    int size = currentParts.size();
    int sizeGIT = githubParts.size();

    // if this is a beta version do not check format "beta.R1"
    if (currentParts[size-2].toUpper().contains("BETA")) {
        beta = true;
        revision = currentParts[size-1];
    }

    if (githubParts[sizeGIT-2].toUpper().contains("BETA")) {
        betaGIT = true;
        revisionGIT = githubParts[sizeGIT-1];
    }

 //   if (beta && !betaGIT)
 //       return false;
    // do not update a beta version, do nothing with revision numbers for now

    // case current 7.4.8 and online 7.4.9 or 7.4.9 and online 7.5
    // use full numbvers, so 7.5.0 and not 7.5
    for (int i = 0; i < qMin(size, sizeGIT); ++i) {
        int currentPart = currentParts.at(i).toInt();
        int githubPart = githubParts.at(i).toInt();
        if (currentPart < githubPart)
            return true;
        else if (currentPart > githubPart)
            return false;
    }
    return size > sizeGIT;
}
//-------------------------------------------------------------------------------------
QString lisemqt::getLatestVersionFromGitHub()
{
    QEventLoop loop;
    manager = new QNetworkAccessManager();
    QNetworkReply *reply = manager->get(QNetworkRequest(QUrl("https://raw.githubusercontent.com/vjetten/openlisem/main_C/ReleaseVersion.txt")));
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();

    QString latestVersion;
    if (reply->error() == QNetworkReply::NoError) {
        QByteArray response = reply->readAll();
        QString content(response);
        QRegularExpression re(R"#(VERSIONNR "([^"]+)")#");
        QRegularExpressionMatch match = re.match(content);
        if (match.hasMatch()) {
            latestVersion = match.captured(1);
        }
    }  else {
        // Handle the network error silently
        qDebug() << "Network error: " << reply->errorString();
    }
    reply->deleteLater();
//qDebug() << "latest version git" << latestVersion;
    return latestVersion;
}
//-------------------------------------------------------------------------------------

void lisemqt::CheckVersion()
{
    QString latestVersion = getLatestVersionFromGitHub();
    if (!latestVersion.isEmpty() && isNewVersionAvailable(latestVersion)) {
qDebug() << "download" << latestVersion;
#ifdef Q_OS_WIN
        downloadPatch(latestVersion);
#elif defined(Q_OS_LINUX)
        QMessageBox::information(nullptr, "Update Available", "A new version is available. "
                                                              "Please download it manually from the GitHub repository.");
#endif

    } else {
        if (latestVersion.isEmpty()) {
            qDebug() << "Cannot check updates online.";
            int ret = QMessageBox::warning(this, "openLISEM","Cannot check updates online.");
        } else {
            QMessageBox::StandardButton replyButton;
            replyButton = QMessageBox::question(nullptr, "openLISEM version check",
                                                QString("No new version is available.\nDo you want to continue checking for new versions?\nYou can activate this again in the Advanced Options."),
                                                QMessageBox::Yes | QMessageBox::No);
            checkforpatch = replyButton == QMessageBox::Yes;
            toolButton_version->setChecked(checkforpatch  );
        }
    }
}




