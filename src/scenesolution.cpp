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

#include "scenesolution.h"

#include "scene.h"
#include "scenemarkerdialog.h"
#include "sceneview_post2d.h"
#include "sceneview_post3d.h"
#include "progressdialog.h"
#include "hermes2d/solver.h"
#include "hermes2d/module.h"
#include "hermes2d/module_agros.h"
#include "hermes2d/problem.h"

template <typename Scalar>
SceneSolution<Scalar>::SceneSolution(FieldInfo* fieldInfo) : m_fieldInfo(fieldInfo)
{
    logMessage("SceneSolution::SceneSolution()");


    //m_timeStep = -1;
    //m_isSolving = false;

    m_meshInitial = NULL;
    m_slnContourView = NULL;
    m_slnScalarView = NULL;
    m_slnVectorXView = NULL;
    m_slnVectorYView = NULL;

    connect(this, SIGNAL(processedRangeContour()), sceneView(), SLOT(processedRangeContour()));
    connect(this, SIGNAL(processedRangeScalar()), sceneView(), SLOT(processedRangeScalar()));
    connect(this, SIGNAL(processedRangeVector()), sceneView(), SLOT(processedRangeVector()));

}

template <typename Scalar>
SceneSolution<Scalar>::~SceneSolution()
{
//    delete m_progressDialog;
//    delete m_progressItemMesh;
//    delete m_progressItemSolve;
//    delete m_progressItemSolveAdaptiveStep;
//    delete m_progressItemProcessView;
}

template <typename Scalar>
void SceneSolution<Scalar>::clear(bool all)
{
    logMessage("SceneSolution::clear()");

//     m_linInitialMeshView.free();
//     m_linSolutionMeshView.free();
//     m_linContourView.free();
//     m_linScalarView.free();
//     m_vecVectorView.free();

    // solution array
    if (all)
    {
        //m_timeStep = -1;

//        for (int i = 0; i < m_solutionArrayList.size(); i++)
//            delete m_solutionArrayList.at(i);
        m_solutionArrayList.clear();

        // mesh
        if (m_meshInitial)
        {
            delete m_meshInitial;
            m_meshInitial = NULL;
        }
    }

    // countour
    if (m_slnContourView)
    {
        delete m_slnContourView;
        m_slnContourView = NULL;
    }
    
    // scalar
    if (m_slnScalarView)
    {
        delete m_slnScalarView;
        m_slnScalarView = NULL;
    }
    
    // vector
    if (m_slnVectorXView)
    {
        delete m_slnVectorXView;
        m_slnVectorXView = NULL;
    }
    if (m_slnVectorYView)
    {
        delete m_slnVectorYView;
        m_slnVectorYView = NULL;
    }

    //m_progressDialog->clear();
}


template <typename Scalar>
void SceneSolution<Scalar>::loadMeshInitial(QDomElement element)
{
    logMessage("SceneSolution::loadMeshInitial()");

    QDomText text = element.childNodes().at(0).toText();

    // write content (saved mesh)
    QString fileName = tempProblemFileName() + ".mesh";
    QByteArray content;
    content.append(text.nodeValue());
    writeStringContentByteArray(fileName, QByteArray::fromBase64(content));

    Hermes::Hermes2D::Mesh *mesh = readMeshFromFile(tempProblemFileName() + ".xml");
    // refineMesh(mesh, true, true);

    setMeshInitial(mesh);
}

template <typename Scalar>
void SceneSolution<Scalar>::saveMeshInitial(QDomDocument *doc, QDomElement element)
{
    logMessage("SceneSolution::saveMeshInitial()");

    assert(0);
//    if (isMeshed())
//    {
//        QString fileName = tempProblemFileName() + ".mesh";

//        writeMeshFromFile(fileName, m_meshInitial);

//        // read content
//        QDomText text = doc->createTextNode(readFileContentByteArray(fileName).toBase64());
//        element.appendChild(text);
//    }
}

template <typename Scalar>
void SceneSolution<Scalar>::loadSolution(QDomElement element)
{
    assert(0); //TODO
//    logMessage("SceneSolution::loadSolution()");

//    Hermes::vector<SolutionArray<Scalar> *> solutionArrayList;

//    // constant solution cannot be saved
//    if (Util::scene()->problemInfo()->analysisType() == AnalysisType_Transient)
//    {
//        Util::scene()->problemInfo()->initialCondition.evaluate(true);

//        SolutionArray<Scalar> *solutionArray = new SolutionArray<Scalar>();
//        // solutionArray->order = new Hermes::Hermes2D::Views::Orderizer();
//        solutionArray->sln = new Hermes::Hermes2D::Solution<Scalar>();
//        // FIXME
//        // solutionArray->sln->set_const(Util::scene()->sceneSolution()->meshInitial(), Util::scene()->problemInfo()->initialCondition.number());
//        solutionArray->adaptiveError = 0.0;
//        solutionArray->adaptiveSteps = 0.0;
//        solutionArray->time = 0.0;

//        solutionArrayList.push_back(solutionArray);
//    }

//    QDomNode n = element.firstChild();
//    while(!n.isNull())
//    {
//        SolutionArray<Scalar> *solutionArray = new SolutionArray<Scalar>();
//        solutionArray->load(n.toElement());

//        // add to the array
//        solutionArrayList.push_back(solutionArray);

//        n = n.nextSibling();
//    }

//    if (solutionArrayList.size() > 0)
//        setSolutionArrayList(solutionArrayList);
}

template <typename Scalar>
void SceneSolution<Scalar>::saveSolution(QDomDocument *doc, QDomElement element)
{
    assert(0); //TODO
//    logMessage("SceneSolution::saveSolution()");

//    if (isSolved())
//    {
//        // constant solution cannot be saved
//        int start = (Util::scene()->problemInfo()->analysisType() != AnalysisType_Transient) ? 0 : 1;

//        for (int i = start; i < timeStepCount(); i++)
//        {
//            QDomNode eleSolution = doc->createElement("solution");
//            m_solutionArrayList.at(i)->save(doc, eleSolution.toElement());
//            element.appendChild(eleSolution);
//        }
//    }
}

template <typename Scalar>
SolutionArray<Scalar> SceneSolution<Scalar>::solutionArray(int i)
{
    logMessage("SceneSolution::solutionArray()");

//    assert(i == 0);
    return m_solutionArrayList.at(i);

//    int currentTimeStep = i;
//    if (isSolved() && currentTimeStep < timeStepCount() * Util::scene()->fieldInfo("TODO")->module()->number_of_solution())
//    {
//        // default
//        if (currentTimeStep == -1)
//            currentTimeStep = m_timeStep * Util::scene()->fieldInfo("TODO")->module()->number_of_solution();

//        if (m_solutionArrayList.at(currentTimeStep))
//            return m_solutionArrayList.at(currentTimeStep);
//    }
//    return NULL;
}

template <typename Scalar>
Hermes::Hermes2D::Solution<Scalar> *SceneSolution<Scalar>::sln(int i)
{
    logMessage("SceneSolution::sln()");

    return solutionArray(i).sln.get();
}

template <typename Scalar>
Hermes::Hermes2D::Space<Scalar> *SceneSolution<Scalar>::space(int i)
{
    logMessage("SceneSolution::space()");

    return solutionArray(i).space.get();
}

//template <typename Scalar>
//double SceneSolution<Scalar>::adaptiveError()
//{
//    assert(0); //TODO
////    logMessage("SceneSolution::adaptiveError()");

////    return (isSolved()) ? m_solutionArrayList.at(m_timeStep * Util::scene()->problemInfo()->module()->number_of_solution())->adaptiveError : 100.0;
//}

//template <typename Scalar>
//int SceneSolution<Scalar>::adaptiveSteps()
//{
//    assert(0); //TODO
////    logMessage("SceneSolution::adaptiveSteps()");

////    return (isSolved()) ? m_solutionArrayList.at(m_timeStep * Util::scene()->problemInfo()->module()->number_of_solution())->adaptiveSteps : 0.0;
//}

template <typename Scalar>
int SceneSolution<Scalar>::findElementInVectorizer(Hermes::Hermes2D::Views::Vectorizer &vec, const Point &point) const
{
    logMessage("SceneSolution::findTriangleInVectorizer()");

    double4* vecVert = vec.get_vertices();
    int3* vecTris = vec.get_triangles();

    for (int i = 0; i < vec.get_num_triangles(); i++)
    {
        bool inElement = true;

        // triangle element
        for (int l = 0; l < 3; l++)
        {
            int k = l + 1;
            if (k == 3)
                k = 0;

            double z = (vecVert[vecTris[i][k]][0] - vecVert[vecTris[i][l]][0]) * (point.y - vecVert[vecTris[i][l]][1]) -
                    (vecVert[vecTris[i][k]][1] - vecVert[vecTris[i][l]][1]) * (point.x - vecVert[vecTris[i][l]][0]);

            if (z < 0)
            {
                inElement = false;
                break;
            }
        }

        if (inElement)
            return i;
    }

    return -1;
}

#define sign(x) (( x > 0 ) - ( x < 0 ))

template <typename Scalar>
int SceneSolution<Scalar>::findElementInMesh(Hermes::Hermes2D::Mesh *mesh, const Point &point) const
{
    logMessage("SceneSolution::findTriangleInMesh()");

    for (int i = 0; i < mesh->get_num_active_elements(); i++)
    {
        Hermes::Hermes2D::Element *element = mesh->get_element_fast(i);

        bool inElement = false;
        int j;
        int npol = (element->is_triangle()) ? 3 : 4;

        for (int i = 0, j = npol-1; i < npol; j = i++) {
            if ((((element->vn[i]->y <= point.y) && (point.y < element->vn[j]->y)) ||
                 ((element->vn[j]->y <= point.y) && (point.y < element->vn[i]->y))) &&
                    (point.x < (element->vn[j]->x - element->vn[i]->x) * (point.y - element->vn[i]->y)
                     / (element->vn[j]->y - element->vn[i]->y) + element->vn[i]->x))
                inElement = !inElement;
        }

        if (inElement)
            return i;
    }
    
    return -1;
}

template <typename Scalar>
void SceneSolution<Scalar>::setMeshInitial(Hermes::Hermes2D::Mesh *meshInitial)
{
    if (m_meshInitial)
    {
        delete m_meshInitial;
    }

    m_meshInitial = meshInitial;
}

template <typename Scalar>
void SceneSolution<Scalar>::setSolutionArrayList(Hermes::vector<SolutionArray<Scalar> > solutionArrayList)
{
    logMessage("SceneSolution::setSolutionArrayList()");

//    for (int i = 0; i < m_solutionArrayList.size(); i++)
//        delete m_solutionArrayList.at(i);
    m_solutionArrayList.clear();

    m_solutionArrayList = solutionArrayList;

    // if (!isSolving())
    //setTimeStep(timeStepCount() - 1);
}

//template <typename Scalar>
//void SceneSolution<Scalar>::setTimeStep(int timeStep, bool showViewProgress)
//{
//    logMessage("SceneSolution::setTimeStep()");

//    m_timeStep = timeStep;
//    if (!isSolved()) return;

//    emit timeStepChanged(showViewProgress);
//}

//template <typename Scalar>
//int SceneSolution<Scalar>::timeStepCount() const
//{
//    logMessage("SceneSolution::timeStepCount()");

//    return (m_solutionArrayList.size() > 0) ? m_solutionArrayList.size() / Util::scene()->fieldInfo("TODO")->module()->number_of_solution() : 0;
//}

//template <typename Scalar>
//double SceneSolution<Scalar>::time() const
//{
//    assert(0); //TODO
////    logMessage("SceneSolution::time()");

////    if (isSolved())
////    {
////        if (m_solutionArrayList.at(m_timeStep * Util::scene()->problemInfo()->module()->number_of_solution())->sln)
////            return m_solutionArrayList.at(m_timeStep * Util::scene()->problemInfo()->module()->number_of_solution())->time;
////    }
////    return 0.0;
//}

template <typename Scalar>
void SceneSolution<Scalar>::setSlnContourView(ViewScalarFilter<Scalar> *slnScalarView)
{
    logMessage("SceneSolution::setSlnContourView()");

    if (m_slnContourView)
    {
        delete m_slnContourView;
        m_slnContourView = NULL;
    }
    
    m_slnContourView = slnScalarView;
    m_linContourView.process_solution(m_slnContourView,
                                      Hermes::Hermes2D::H2D_FN_VAL_0,
                                      Util::config()->linearizerQuality);

    // deformed shape
    if (Util::config()->deformContour)
        m_fieldInfo->module()->deform_shape(m_linContourView.get_vertices(), m_linContourView.get_num_vertices());
}

template <typename Scalar>
void SceneSolution<Scalar>::setSlnScalarView(ViewScalarFilter<Scalar> *slnScalarView)
{
    logMessage("SceneSolution::setSlnScalarView()");

    if (m_slnScalarView)
    {
        delete m_slnScalarView;
        m_slnScalarView = NULL;
    }
    
    m_slnScalarView = slnScalarView;
    // QTime time;
    // time.start();
    m_linScalarView.process_solution(m_slnScalarView,
                                     Hermes::Hermes2D::H2D_FN_VAL_0,
                                     Util::config()->linearizerQuality);

    // deformed shape
    if (Util::config()->deformScalar)
        m_fieldInfo->module()->deform_shape(m_linScalarView.get_vertices(),
                                                             m_linScalarView.get_num_vertices());
}

template <typename Scalar>
void SceneSolution<Scalar>::setSlnVectorView(ViewScalarFilter<Scalar> *slnVectorXView, ViewScalarFilter<Scalar> *slnVectorYView)
{
    logMessage("SceneSolution::setSlnVectorView()");

    if (m_slnVectorXView)
    {
        delete m_slnVectorXView;
        m_slnVectorXView = NULL;
    }
    if (m_slnVectorYView)
    {
        delete m_slnVectorYView;
        m_slnVectorYView = NULL;
    }
    
    m_slnVectorXView = slnVectorXView;
    m_slnVectorYView = slnVectorYView;

    m_vecVectorView.process_solution(m_slnVectorXView, m_slnVectorYView,
                                     Hermes::Hermes2D::H2D_FN_VAL_0, Hermes::Hermes2D::H2D_FN_VAL_0,
                                     Hermes::Hermes2D::Views::HERMES_EPS_LOW);

    // deformed shape
    if (Util::config()->deformVector)
        m_fieldInfo->module()->deform_shape(m_vecVectorView.get_vertices(),
                                                             m_vecVectorView.get_num_vertices());
}

template <typename Scalar>
void SceneSolution<Scalar>::setOrderView(shared_ptr<Hermes::Hermes2D::Space<Scalar> > space)
{
    logMessage("SceneSolution::setOrderView()");

    m_orderView.process_space(space.get());
}

template <typename Scalar>
void SceneSolution<Scalar>::processView(bool showViewProgress)
{
    int step = 0;

    // process order
    Util::scene()->activeSceneSolution()->processOrder();

    if (sceneView()->sceneViewSettings().showSolutionMesh)
    {
        step++;
        //emit message(tr("Processing solution mesh cache"), false, step);
        Util::scene()->activeSceneSolution()->processSolutionMesh();
    }
    if (sceneView()->sceneViewSettings().showContours)
    {
        step++;
        //emit message(tr("Processing contour view cache"), false, step);
        Util::scene()->activeSceneSolution()->processRangeContour();
    }
    if (sceneView()->sceneViewSettings().postprocessorShow == SceneViewPostprocessorShow_ScalarView ||
            sceneView()->sceneViewSettings().postprocessorShow == SceneViewPostprocessorShow_ScalarView3D ||
            sceneView()->sceneViewSettings().postprocessorShow == SceneViewPostprocessorShow_ScalarView3DSolid)
    {
        step++;
        //emit message(tr("Processing scalar view cache"), false, step);
        cout << "process Range Scalar" << endl;
        Util::scene()->activeSceneSolution()->processRangeScalar();
    }
    if (sceneView()->sceneViewSettings().showVectors)
    {
        step++;
        //emit message(tr("Processing vector view cache"), false, step);
        Util::scene()->activeSceneSolution()->processRangeVector();
    }


//    if (showViewProgress)
//    {
//        m_progressDialog->clear();
//        m_progressDialog->appendProgressItem(m_progressItemProcessView);
//        m_progressDialog->run(showViewProgress);
//    }
//    else
//    {
//        m_progressItemProcessView->setSteps();
//        m_progressItemProcessView->run(true);
//    }
}

template <typename Scalar>
void SceneSolution<Scalar>::processSolutionMesh()
{
    logMessage("SceneSolution::processSolutionMesh()");

    if (Util::problem()->isSolved())
    {
        InitialCondition<double> initial(sln()->get_mesh(), 0.0);
        m_linSolutionMeshView.process_solution(&initial);

        emit processedSolutionMesh();
    }
}

template <typename Scalar>
void SceneSolution<Scalar>::processRangeContour()
{
    logMessage("SceneSolution::processRangeContour()");

    if (Util::problem()->isSolved() && sceneView()->sceneViewSettings().contourPhysicFieldVariable != "")
    {
        ViewScalarFilter<Scalar> *viewScalarFilter;
        if (m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().contourPhysicFieldVariable)->is_scalar)
            viewScalarFilter = m_fieldInfo->module()->view_scalar_filter(m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().contourPhysicFieldVariable),
                                                                                          PhysicFieldVariableComp_Scalar);
        else
            viewScalarFilter = m_fieldInfo->module()->view_scalar_filter(m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().contourPhysicFieldVariable),
                                                                                          PhysicFieldVariableComp_Magnitude);

        setSlnContourView(viewScalarFilter);
        emit processedRangeContour();
    }
}

template <typename Scalar>
void SceneSolution<Scalar>::processRangeScalar()
{
    logMessage("SceneSolution::processRangeScalar()");

    if (Util::problem()->isSolved() && sceneView()->sceneViewSettings().scalarPhysicFieldVariable != "")
    {
        ViewScalarFilter<Scalar> *viewScalarFilter = m_fieldInfo->module()->view_scalar_filter(m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().scalarPhysicFieldVariable),
                                                                                                                sceneView()->sceneViewSettings().scalarPhysicFieldVariableComp);

        setSlnScalarView(viewScalarFilter);
        emit processedRangeScalar();
    }
}

template <typename Scalar>
void SceneSolution<Scalar>::processRangeVector()
{
    logMessage("SceneSolution::processRangeVector()");

    if (Util::problem()->isSolved() && sceneView()->sceneViewSettings().vectorPhysicFieldVariable != "")
    {
        ViewScalarFilter<Scalar> *viewVectorXFilter = m_fieldInfo->module()->view_scalar_filter(m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().vectorPhysicFieldVariable),
                                                                                                                 PhysicFieldVariableComp_X);

        ViewScalarFilter<Scalar> *viewVectorYFilter = m_fieldInfo->module()->view_scalar_filter(m_fieldInfo->module()->get_variable(sceneView()->sceneViewSettings().vectorPhysicFieldVariable),
                                                                                                                 PhysicFieldVariableComp_Y);

        setSlnVectorView(viewVectorXFilter, viewVectorYFilter);
        emit processedRangeVector();
    }
}

template <typename Scalar>
void SceneSolution<Scalar>::processOrder()
{
    logMessage("SceneSolution::processOrder()");

    if (Util::problem()->isSolved())
        setOrderView(m_solutionArrayList.at(0).space);
                //TODO timedependence rpoblemsm_timeStep * Util::scene()->problemInfo()->module()->number_of_solution())->space);
}

template <typename Scalar>
void SceneSolution<Scalar>::setSolutionArray(QList<SolutionArray<Scalar> > solutionArrays)
{
    assert(solutionArrays.size() == fieldInfo()->module()->number_of_solution());

    for(int i = 0; i < solutionArrays.size(); i++)
        m_solutionArrayList.push_back(solutionArrays[i]);
}

template class SceneSolution<double>;
