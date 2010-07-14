// This file is part of Agros2D.
//
// Agros2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Agros2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Agros2D.  If not, see <http://www.gnu.org/licenses/>.
//
// hp-FEM group (http://hpfem.org/)
// University of Nevada, Reno (UNR) and University of West Bohemia, Pilsen
// Email: agros2d@googlegroups.com, home page: http://hpfem.org/agros2d/

#include "logdialog.h"

LogDialog::LogDialog(QWidget *parent) : QDialog(parent)
{
    setWindowIcon(icon("log"));
    setWindowTitle(tr("Progress log"));
    setWindowFlags(Qt::Window);

    createControls();

    QSettings settings;
    setMinimumSize(sizeHint());
    restoreGeometry(settings.value("LogDialog/Geometry", saveGeometry()).toByteArray());
}

LogDialog::~LogDialog()
{
    QSettings settings;
    settings.setValue("LogDialog/Geometry", saveGeometry());

    delete btnSaveLog;
    delete btnClose;
}

void LogDialog::createControls()
{
    lstMessages = new QTextEdit(this);
    lstMessages->setReadOnly(true);

    btnSaveLog = new QPushButton(tr("Save log"));
    connect(btnSaveLog, SIGNAL(clicked()), this, SLOT(doSaveLog()));

    btnClose = new QPushButton(tr("Close"));
    connect(btnClose, SIGNAL(clicked()), this, SLOT(doClose()));

    QHBoxLayout *layoutButtons = new QHBoxLayout();
    layoutButtons->addStretch();
    layoutButtons->addWidget(btnSaveLog);
    layoutButtons->addWidget(btnClose);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->addWidget(lstMessages);
    layout->addLayout(layoutButtons);

    setLayout(layout);
}

void LogDialog::showDialog()
{
    loadProgressLog();

    show();
    activateWindow();
    raise();
}

void LogDialog::loadProgressLog()
{
    lstMessages->clear();

    QFile file(tempProblemDir() + "/messages.log");
    if (file.exists())
    {
        if (file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            QString message;
            QTextStream messages(&file);
            while (!messages.atEnd())
            {
                message = messages.readLine();
                if (message[1].isNumber())
                    lstMessages->setTextColor(QColor(Qt::gray));
                else
                    lstMessages->setTextColor(QColor(Qt::black));

                lstMessages->insertPlainText(message + "\n");
            }
        }
        btnSaveLog->setEnabled(true);
    }
    else
    {
        lstMessages->setTextColor(QColor(Qt::gray));
        lstMessages->insertPlainText(tr("No messages..."));
        btnSaveLog->setEnabled(false);
    }
}


void LogDialog::doClose()
{
    hide();
}

void LogDialog::doSaveLog()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save log"), "data", tr("Log files (*.log)"));

    QFile file(fileName);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QTextStream messages(&file);
        messages << lstMessages->toPlainText();
    }
}
