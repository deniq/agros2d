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

#ifndef OPTIONSDIALOG_H
#define OPTIONSDIALOG_H

#include "util.h"

struct SceneViewSettings;
class ScriptEditor;

class LineEditDouble;
class ValueLineEdit;

class ConfigComputerDialog : public QDialog
{
    Q_OBJECT
public:
    ConfigComputerDialog(QWidget *parent);

private slots:
    void doCurrentItemChanged(QListWidgetItem *current, QListWidgetItem *previous);

    void doAccept();
    void doReject();

private:
    QListWidget *lstView;
    QStackedWidget *pages;
    QWidget *panMain;
    QWidget *panSolver;

    // main
    QComboBox *cmbGUIStyle;
    QComboBox *cmbLanguage;
    // LineEditDouble *txtCollaborationServerURL;

    // show result in line edit value widget
    QCheckBox *chkLineEditValueShowResult;

    // log std out
    QCheckBox *chkLogStdOut;

    // development
    QCheckBox *chkDiscreteSaveMatrixRHS;
    QComboBox *cmbDumpFormat;

    // cache
    QSpinBox *txtCacheSize;

    // threads
    QSpinBox *txtNumOfThreads;

    // grid
    QCheckBox *chkShowGrid;

    // scene font
    QComboBox *cmbRulersFont;
    QSpinBox *txtRulersFontSizes;
    QComboBox *cmbPostFont;
    QSpinBox *txtPostFontSizes;

    // workspace other
    QCheckBox *chkShowAxes;
    QCheckBox *chkShowRulers;

    void load();
    void save();

    void createControls();
    QWidget *createMainWidget();
    QWidget *createSolverWidget();

    void fillComboBoxPhysicField(QComboBox *cmbPhysicField);
};


#endif // OPTIONSDIALOG_H
