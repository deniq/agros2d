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

#include "meshgenerator.h"

#include "gui.h"
#include "scene.h"

#include "scenebasic.h"
#include "scenenode.h"
#include "sceneedge.h"
#include "scenelabel.h"

#include "sceneview_common.h"
#include "scenemarker.h"
#include "scenemarkerdialog.h"
#include "logview.h"

#include "hermes2d/module.h"
#include "hermes2d/module_agros.h"
#include "hermes2d/problem.h"

MeshGeneratorTriangle::MeshGeneratorTriangle() : QObject()
{    
}

bool MeshGeneratorTriangle::mesh()
{
    m_isError = false;

    QFile::remove(tempProblemFileName() + ".mesh");

    // create triangle files
    if (writeToTriangle())
    {
        Util::log()->printDebug(tr("Mesh generator"), tr("Poly file was created"));

        // exec triangle
        QProcess processTriangle;
        processTriangle.setStandardOutputFile(tempProblemFileName() + ".triangle.out");
        processTriangle.setStandardErrorFile(tempProblemFileName() + ".triangle.err");
        connect(&processTriangle, SIGNAL(finished(int)), this, SLOT(meshTriangleCreated(int)));

        QString triangleBinary = "triangle";
        if (QFile::exists(QApplication::applicationDirPath() + QDir::separator() + "triangle.exe"))
            triangleBinary = "\"" + QApplication::applicationDirPath() + QDir::separator() + "triangle.exe\"";
        if (QFile::exists(QApplication::applicationDirPath() + QDir::separator() + "triangle"))
            triangleBinary = QApplication::applicationDirPath() + QDir::separator() + "triangle";

        processTriangle.start(QString(Util::config()->commandTriangle).
                              arg(triangleBinary).
                              arg(tempProblemFileName()));

        if (!processTriangle.waitForStarted(100000))
        {
            Util::log()->printError(tr("Mesh generator"), tr("could not start Triangle"));
            processTriangle.kill();

            m_isError = true;
        }
        else
        {
            // copy triangle files
            if ((!Util::config()->deleteTriangleMeshFiles) && (!Util::scene()->problemInfo()->fileName.isEmpty()))
            {
                QFileInfo fileInfoOrig(Util::scene()->problemInfo()->fileName);

                QFile::copy(tempProblemFileName() + ".poly", fileInfoOrig.absolutePath() + "/" + fileInfoOrig.baseName() + ".poly");
                QFile::copy(tempProblemFileName() + ".node", fileInfoOrig.absolutePath() + "/" + fileInfoOrig.baseName() + ".node");
                QFile::copy(tempProblemFileName() + ".edge", fileInfoOrig.absolutePath() + "/" + fileInfoOrig.baseName() + ".edge");
                QFile::copy(tempProblemFileName() + ".ele", fileInfoOrig.absolutePath() + "/" + fileInfoOrig.baseName() + ".ele");
            }

            while (!processTriangle.waitForFinished()) {}
        }
    }
    else
    {
        m_isError = true;
    }

    return !m_isError;
}

void MeshGeneratorTriangle::meshTriangleCreated(int exitCode)
{
    if (exitCode == 0)
    {
        Util::log()->printMessage(tr("Mesh generator"), tr("mesh files were created"));

        // convert triangle mesh to hermes mesh
        if (triangleToHermes2D())
        {
            Util::log()->printMessage(tr("Mesh generator"), tr("mesh was converted to Hermes2D mesh file"));

            // copy triangle files
            if ((!Util::config()->deleteHermes2DMeshFile) && (!Util::scene()->problemInfo()->fileName.isEmpty()))
            {
                QFileInfo fileInfoOrig(Util::scene()->problemInfo()->fileName);

                QFile::copy(tempProblemFileName() + ".mesh", fileInfoOrig.absolutePath() + "/" + fileInfoOrig.baseName() + ".mesh");
            }

            //  remove triangle temp files
            QFile::remove(tempProblemFileName() + ".poly");
            QFile::remove(tempProblemFileName() + ".node");
            QFile::remove(tempProblemFileName() + ".edge");
            QFile::remove(tempProblemFileName() + ".ele");
            QFile::remove(tempProblemFileName() + ".neigh");
            QFile::remove(tempProblemFileName() + ".triangle.out");
            QFile::remove(tempProblemFileName() + ".triangle.err");
            Util::log()->printMessage(tr("Mesh generator"), tr("mesh files were deleted"));

            // load mesh
            QMap<FieldInfo*, Hermes::Hermes2D::Mesh*> meshes = readMeshesFromFile(tempProblemFileName() + ".xml");

            // FIXME: jinak
            Util::problem()->setMeshesInitial(meshes);
        }
        else
        {
            m_isError = true;
            QFile::remove(Util::scene()->problemInfo()->fileName + ".mesh");
        }
    }
    else
    {
        m_isError = true;
        QString errorMessage = readFileContent(Util::scene()->problemInfo()->fileName + ".triangle.out");
        Util::log()->printError(tr("Mesh generator"), errorMessage);
    }
}

bool MeshGeneratorTriangle::writeToTriangle()
{
    // basic check
    if (Util::scene()->nodes->length() < 3)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of nodes (%1 < 3)").arg(Util::scene()->nodes->length()));
        return false;
    }
    if (Util::scene()->edges->length() < 3)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of edges (%1 < 3)").arg(Util::scene()->edges->length()));
        return false;
    }

    // save current locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    QDir dir;
    dir.mkdir(QDir::temp().absolutePath() + "/agros2d");
    QFile file(tempProblemFileName() + ".poly");

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not create Triangle poly mesh file (%1)").arg(file.errorString()));
        return false;
    }
    QTextStream out(&file);


    // nodes
    QString outNodes;
    int nodesCount = 0;
    for (int i = 0; i<Util::scene()->nodes->length(); i++)
    {
        outNodes += QString("%1  %2  %3  %4\n").
                arg(i).
                arg(Util::scene()->nodes->at(i)->point.x, 0, 'f', 10).
                arg(Util::scene()->nodes->at(i)->point.y, 0, 'f', 10).
                arg(0);
        nodesCount++;
    }

    // edges
    QString outEdges;
    int edgesCount = 0;
    for (int i = 0; i<Util::scene()->edges->length(); i++)
    {
        if (Util::scene()->edges->at(i)->angle == 0)
        {
            // line
            outEdges += QString("%1  %2  %3  %4\n").
                    arg(edgesCount).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeStart)).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeEnd)).
                    arg(i+1);
            edgesCount++;
        }
        else
        {
            // arc
            // add pseudo nodes
            Point center = Util::scene()->edges->at(i)->center();
            double radius = Util::scene()->edges->at(i)->radius();
            double startAngle = atan2(center.y - Util::scene()->edges->at(i)->nodeStart->point.y,
                                      center.x - Util::scene()->edges->at(i)->nodeStart->point.x) - M_PI;

            int segments = Util::scene()->edges->at(i)->segments();
            double theta = deg2rad(Util::scene()->edges->at(i)->angle) / double(segments);

            int nodeStartIndex = 0;
            int nodeEndIndex = 0;
            for (int j = 0; j < segments; j++)
            {
                double arc = startAngle + j*theta;

                double x = radius * cos(arc);
                double y = radius * sin(arc);

                nodeEndIndex = nodesCount+1;
                if (j == 0)
                {
                    nodeStartIndex = Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeStart);
                    nodeEndIndex = nodesCount;
                }
                if (j == segments - 1)
                {
                    nodeEndIndex = Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeEnd);
                }
                if ((j > 0) && (j < segments))
                {
                    outNodes += QString("%1  %2  %3  %4\n").
                            arg(nodesCount).
                            arg(center.x + x, 0, 'f', 10).
                            arg(center.y + y, 0, 'f', 10).
                            arg(0);
                    nodesCount++;
                }
                outEdges += QString("%1  %2  %3  %4\n").
                        arg(edgesCount).
                        arg(nodeStartIndex).
                        arg(nodeEndIndex).
                        arg(i+1);
                edgesCount++;
                nodeStartIndex = nodeEndIndex;
            }
        }
    }

    // holes
    int holesCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
        if (label->markersCount() == 0)
            holesCount++;

    QString outHoles = QString("%1\n").arg(holesCount);
    holesCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
    {
        if (label->markersCount() == 0)
        {
            outHoles += QString("%1  %2  %3\n").
                    arg(holesCount).
                    // arg(Util::scene()->labels->items().indexOf(label) + 1).
                    arg(label->point.x, 0, 'f', 10).
                    arg(label->point.y, 0, 'f', 10);

            holesCount++;
        }
    }

    // labels
    QString outLabels;
    int labelsCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
    {
        if (label->markersCount() > 0)
        {
            outLabels += QString("%1  %2  %3  %4  %5\n").
                    arg(labelsCount).
                    arg(label->point.x, 0, 'f', 10).
                    arg(label->point.y, 0, 'f', 10).
                    // arg(labelsCount + 1). // triangle returns zero region number for areas without marker, markers must start from 1
                    arg(Util::scene()->labels->items().indexOf(label) + 1).
                    arg(label->area);
            labelsCount++;
        }
    }

    outNodes.insert(0, QString("%1 2 0 1\n").
                    arg(nodesCount)); // + additional Util::scene()->nodes
    out << outNodes;
    outEdges.insert(0, QString("%1 1\n").
                    arg(edgesCount)); // + additional edges
    out << outEdges;
    out << outHoles;
    outLabels.insert(0, QString("%1 1\n").
                     arg(labelsCount)); // - holes
    out << outLabels;

    file.waitForBytesWritten(0);
    file.close();

    // set system locale
    setlocale(LC_NUMERIC, plocale);

    return true;
}

bool MeshGeneratorTriangle::triangleToHermes2D()
{
    int k;

    // save current locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    QDomDocument doc;
    QDomProcessingInstruction instr = doc.createProcessingInstruction("xml", "version='1.0' encoding='UTF-8' standalone='no'");
    doc.appendChild(instr);

    // main document
    QDomElement eleMesh = doc.createElement("domain:domain");
    eleMesh.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    eleMesh.setAttribute("xmlns:domain", "XMLSubdomains");
    eleMesh.setAttribute("xsi:schemaLocation", QString("XMLSubdomains %1/subdomains_h2d_xml.xsd").arg(datadir() + "/resources/xsd"));
    doc.appendChild(eleMesh);

    QDomElement eleVertices = doc.createElement("vertices");
    eleMesh.appendChild(eleVertices);
    QDomElement eleElements = doc.createElement("elements");
    eleMesh.appendChild(eleElements);
    QDomElement eleEdges = doc.createElement("edges");
    eleMesh.appendChild(eleEdges);
    QDomElement eleCurves = doc.createElement("curves");
    eleMesh.appendChild(eleCurves);
    QDomElement eleSubdomains = doc.createElement("subdomains");
    eleMesh.appendChild(eleSubdomains);

    QFile fileNode(tempProblemFileName() + ".node");
    if (!fileNode.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not read Triangle node file"));
        return false;
    }
    QTextStream inNode(&fileNode);

    QFile fileEdge(tempProblemFileName() + ".edge");
    if (!fileEdge.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not read Triangle edge file"));
        return false;
    }
    QTextStream inEdge(&fileEdge);

    QFile fileEle(tempProblemFileName() + ".ele");
    if (!fileEle.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not read Triangle ele file"));
        return false;
    }
    QTextStream inEle(&fileEle);

    QFile fileNeigh(tempProblemFileName() + ".neigh");
    if (!fileNeigh.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not read Triangle neigh file"));
        return false;
    }
    QTextStream inNeigh(&fileNeigh);

    // triangle nodes
    sscanf(inNode.readLine().toStdString().c_str(), "%i", &k);
    QList<Point> nodeList;
    for (int i = 0; i<k; i++)
    {
        int marker, n;
        double x, y;

        sscanf(inNode.readLine().toStdString().c_str(), "%i   %lf %lf %i", &n, &x, &y, &marker);
        nodeList.append(Point(x, y));
    }

    // triangle edges
    sscanf(inEdge.readLine().toStdString().c_str(), "%i", &k);
    QList<MeshEdge> edgeList;
    QSet<int> edgeMarkersCheck;
    for (int i = 0; i<k; i++)
    {
        int node[2];
        int marker, n;

        sscanf(inEdge.readLine().toStdString().c_str(), "%i	%i	%i	%i", &n, &node[0], &node[1], &marker);

        edgeList.append(MeshEdge(node[0], node[1], marker - 1)); // marker conversion from triangle, where it starts from 1
    }
    int edgeCountLinear = edgeList.count();

    FieldInfo *fieldInfoTMP;
    foreach (FieldInfo *fieldInfo, Util::scene()->fieldInfos())
    {
        fieldInfoTMP = fieldInfo;
    }

    // check for multiple edge markers
    foreach (SceneEdge *edge, Util::scene()->edges->items())
    {
        // if (Util::scene()->boundaries.indexOf(edge->getMarker(fieldInfoTMP) > 0)
        {
            // if (!edgeMarkersCheck.contains(i))
            {
                // emit message(tr("Edge marker '%1' is not present in mesh file.").
                //              arg(i), true, 0);
                // return false;
            }
        }
    }
    // edgeMarkersCheck.clear();

    // TODO: move
    // no edge marker
    /*
    if (edgeCountLinear < 1)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of edge markers"));
        return false;
    }
    */

    // triangle elements
    sscanf(inEle.readLine().toStdString().c_str(), "%i", &k);
    QList<MeshElement> elementList;
    QSet<int> labelMarkersCheck;
    for (int i = 0; i < k; i++)
    {
        int node[6];
        int marker, n;

        sscanf(inEle.readLine().toStdString().c_str(), "%i	%i	%i	%i	%i	%i	%i	%i",
               &n, &node[0], &node[1], &node[2], &node[3], &node[4], &node[5], &marker);
        assert(marker > 0);

        if (Util::scene()->problemInfo()->meshType == MeshType_Triangle ||
                Util::scene()->problemInfo()->meshType == MeshType_QuadJoin ||
                Util::scene()->problemInfo()->meshType == MeshType_QuadRoughDivision)
        {
            elementList.append(MeshElement(node[0], node[1], node[2], marker - 1)); // marker conversion from triangle, where it starts from 1
        }

        if (Util::scene()->problemInfo()->meshType == MeshType_QuadFineDivision)
        {
            // add additional node
            nodeList.append(Point((nodeList[node[0]].x + nodeList[node[1]].x + nodeList[node[2]].x) / 3.0,
                                  (nodeList[node[0]].y + nodeList[node[1]].y + nodeList[node[2]].y) / 3.0));
            // add three quad elements
            elementList.append(MeshElement(node[4], node[0], node[5], nodeList.count() - 1, marker - 1)); // marker conversion from triangle, where it starts from 1
            elementList.append(MeshElement(node[5], node[1], node[3], nodeList.count() - 1, marker - 1)); // marker conversion from triangle, where it starts from 1
            elementList.append(MeshElement(node[3], node[2], node[4], nodeList.count() - 1, marker - 1)); // marker conversion from triangle, where it starts from 1
        }

        if (marker == 0)
        {
            Util::log()->printError(tr("Mesh generator"), tr("some areas have no label marker"));
            return false;
        }

        labelMarkersCheck.insert(marker - 1);
    }
    int elementCountLinear = elementList.count();

    // check for multiple label markers
    foreach (SceneLabel *label, Util::scene()->labels->items())
    {
        // if (Util::scene()->materials.indexOf(Util::scene()->labels[i]->material) > 0)
        {
            //            if (!labelMarkersCheck.contains(i + 1))
            //            {
            //                emit message(tr("Label marker '%1' is not present in mesh file (multiple label markers in one area).").
            //                             arg(i), true, 0);
            //                return false;
            //            }
        }
    }
    labelMarkersCheck.clear();

    // TODO: move
    /*
    if (elementList.count() < 1)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of label markers"));
        return false;
    }
    */

    // triangle neigh
    sscanf(inNeigh.readLine().toStdString().c_str(), "%i", &k);
    for (int i = 0; i < k; i++)
    {
        int n;
        int ele_1, ele_2, ele_3;

        sscanf(inNeigh.readLine().toStdString().c_str(), "%i	%i	%i	%i", &n, &ele_1, &ele_2, &ele_3);
        elementList[i].neigh[0] = ele_1;
        elementList[i].neigh[1] = ele_2;
        elementList[i].neigh[2] = ele_3;
    }

    // heterogeneous mesh
    // element division
    if (Util::scene()->problemInfo()->meshType == MeshType_QuadFineDivision)
    {
        for (int i = 0; i < edgeCountLinear; i++)
        {
            if (edgeList[i].marker != -1)
            {
                for (int j = 0; j < elementList.count() / 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        if (edgeList[i].node[0] == elementList[3*j + k].node[1] && edgeList[i].node[1] == elementList[3*j + (k + 1) % 3].node[1])
                        {
                            edgeList.append(MeshEdge(edgeList[i].node[0], elementList[3*j + (k + 1) % 3].node[0], edgeList[i].marker));
                            edgeList[i].node[0] = elementList[3*j + (k + 1) % 3].node[0];
                        }
                    }
                }
            }
        }
    }

    if (Util::scene()->problemInfo()->meshType == MeshType_QuadRoughDivision)
    {
        for (int i = 0; i < elementCountLinear; i++)
        {
            // check same material
            if (elementList[i].isActive &&
                    (elementList[i].neigh[0] != -1) &&
                    (elementList[i].neigh[1] != -1) &&
                    (elementList[i].neigh[2] != -1) &&
                    elementList[elementList[i].neigh[0]].isActive &&
                    elementList[elementList[i].neigh[1]].isActive &&
                    elementList[elementList[i].neigh[2]].isActive &&
                    (elementList[i].marker == elementList[elementList[i].neigh[0]].marker) &&
                    (elementList[i].marker == elementList[elementList[i].neigh[1]].marker) &&
                    (elementList[i].marker == elementList[elementList[i].neigh[2]].marker))
            {
                // add additional node
                nodeList.append(Point((nodeList[elementList[i].node[0]].x + nodeList[elementList[i].node[1]].x + nodeList[elementList[i].node[2]].x) / 3.0,
                                      (nodeList[elementList[i].node[0]].y + nodeList[elementList[i].node[1]].y + nodeList[elementList[i].node[2]].y) / 3.0));

                // add three quad elements
                for (int nd = 0; nd < 3; nd++)
                    for (int neigh = 0; neigh < 3; neigh++)
                        for (int neigh_nd = 0; neigh_nd < 3; neigh_nd++)
                            if ((elementList[i].node[(nd + 0) % 3] == elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 1) % 3]) &&
                                    (elementList[i].node[(nd + 1) % 3] == elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 0) % 3]))
                                elementList.append(MeshElement(elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 1) % 3],
                                                               elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 2) % 3],
                                                               elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 0) % 3],
                                                               nodeList.count() - 1, elementList[i].marker));

                elementList[i].isUsed = false;
                elementList[i].isActive = false;
                for (int k = 0; k < 3; k++)
                {
                    elementList[elementList[i].neigh[k]].isUsed = false;
                    elementList[elementList[i].neigh[k]].isActive = false;
                }
            }
        }
    }

    if (Util::scene()->problemInfo()->meshType == MeshType_QuadJoin)
    {
        for (int i = 0; i < elementCountLinear; i++)
        {
            // check same material
            if (elementList[i].isActive)
            {
                // add quad elements
                for (int nd = 0; nd < 3; nd++)
                    for (int neigh = 0; neigh < 3; neigh++)
                        for (int neigh_nd = 0; neigh_nd < 3; neigh_nd++)
                            if (elementList[i].isActive &&
                                    elementList[i].neigh[neigh] != -1 &&
                                    elementList[elementList[i].neigh[neigh]].isActive &&
                                    elementList[i].marker == elementList[elementList[i].neigh[neigh]].marker)
                                if ((elementList[i].node[(nd + 0) % 3] == elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 1) % 3]) &&
                                        (elementList[i].node[(nd + 1) % 3] == elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 0) % 3]))
                                {
                                    int tmp_node[3];
                                    for (int k = 0; k < 3; k++)
                                        tmp_node[k] = elementList[i].node[k];

                                    Point quad_check[4];
                                    quad_check[0] = nodeList[tmp_node[(nd + 1) % 3]];
                                    quad_check[1] = nodeList[tmp_node[(nd + 2) % 3]];
                                    quad_check[2] = nodeList[tmp_node[(nd + 0) % 3]];
                                    quad_check[3] = nodeList[elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 2) % 3]];

                                    if ((!Hermes::Hermes2D::Mesh::same_line(quad_check[0].x, quad_check[0].y, quad_check[1].x, quad_check[1].y, quad_check[2].x, quad_check[2].y)) &&
                                            (!Hermes::Hermes2D::Mesh::same_line(quad_check[0].x, quad_check[0].y, quad_check[1].x, quad_check[1].y, quad_check[3].x, quad_check[3].y)) &&
                                            (!Hermes::Hermes2D::Mesh::same_line(quad_check[0].x, quad_check[0].y, quad_check[2].x, quad_check[2].y, quad_check[3].x, quad_check[3].y)) &&
                                            (!Hermes::Hermes2D::Mesh::same_line(quad_check[1].x, quad_check[1].y, quad_check[2].x, quad_check[2].y, quad_check[3].x, quad_check[3].y)) &&
                                            Hermes::Hermes2D::Mesh::is_convex(quad_check[1].x - quad_check[0].x, quad_check[1].y - quad_check[0].y, quad_check[2].x - quad_check[0].x, quad_check[2].y - quad_check[0].y) &&
                                            Hermes::Hermes2D::Mesh::is_convex(quad_check[2].x - quad_check[0].x, quad_check[2].y - quad_check[0].y, quad_check[3].x - quad_check[0].x, quad_check[3].y - quad_check[0].y) &&
                                            Hermes::Hermes2D::Mesh::is_convex(quad_check[2].x - quad_check[1].x, quad_check[2].y - quad_check[1].y, quad_check[3].x - quad_check[1].x, quad_check[3].y - quad_check[1].y) &&
                                            Hermes::Hermes2D::Mesh::is_convex(quad_check[3].x - quad_check[1].x, quad_check[3].y - quad_check[1].y, quad_check[0].x - quad_check[1].x, quad_check[0].y - quad_check[1].y))
                                    {
                                        // regularity check
                                        bool regular = true;
                                        for (int k = 0; k < 4; k++)
                                        {
                                            double length_1 = (quad_check[k] - quad_check[(k + 1) % 4]).magnitude();
                                            double length_2 = (quad_check[(k + 1) % 4] - quad_check[(k + 2) % 4]).magnitude();
                                            double length_together = (quad_check[k] - quad_check[(k + 2) % 4]).magnitude();

                                            if ((length_1 + length_2) / length_together < 1.03)
                                                regular = false;
                                        }

                                        if (!regular)
                                            break;

                                        elementList[i].node[0] = tmp_node[(nd + 1) % 3];
                                        elementList[i].node[1] = tmp_node[(nd + 2) % 3];
                                        elementList[i].node[2] = tmp_node[(nd + 0) % 3];
                                        elementList[i].node[3] = elementList[elementList[i].neigh[neigh]].node[(neigh_nd + 2) % 3];

                                        elementList[i].isActive = false;

                                        elementList[elementList[i].neigh[neigh]].isUsed = false;
                                        elementList[elementList[i].neigh[neigh]].isActive = false;

                                        break;
                                    }
                                }
            }
        }
    }

    // edges
    int countEdges = 0;
    for (int i = 0; i < edgeList.count(); i++)
    {
        if (edgeList[i].isUsed && edgeList[i].marker != -1)
        {
            //TODO - no "inner edges" in new xml mesh file format - remove?
            // inner edge marker (minus markers are ignored)
//            int marker = - (edgeList[i].marker-1);
//            foreach (FieldInfo *fieldInfo, Util::scene()->fieldInfos())
//                if (Util::scene()->edges->at(edgeList[i].marker-1)->getMarker(fieldInfoTMP)
//                        != SceneBoundaryContainer::getNone(fieldInfoTMP))
//                {
//                    // boundary marker
//                    marker = edgeList[i].marker;

//                    break;
//                }

//            countEdges++;

            int marker = edgeList[i].marker;


            //assert(countEdges == i+1);
            QDomElement eleEdge = doc.createElement("edge");
            eleEdge.setAttribute("v1", edgeList[i].node[0]);
            eleEdge.setAttribute("v2", edgeList[i].node[1]);
            eleEdge.setAttribute("i", i);
            eleEdge.setAttribute("marker", marker);

            eleEdges.appendChild(eleEdge);
        }
    }

    // curves
    int countCurves = 0;
    if (Util::config()->curvilinearElements)
    {
        for (int i = 0; i<edgeList.count(); i++)
        {
            if (edgeList[i].marker != -1)
            {
                // curve
                if (Util::scene()->edges->at(edgeList[i].marker)->angle > 0.0)
                {
                    countCurves++;
                    int segments = Util::scene()->edges->at(edgeList[i].marker)->segments();

                    // subdivision angle and chord
                    double theta = deg2rad(Util::scene()->edges->at(edgeList[i].marker)->angle) / double(segments);
                    double chord = 2 * Util::scene()->edges->at(edgeList[i].marker)->radius() * sin(theta / 2.0);

                    // length of short chord
                    double chordShort = (nodeList[edgeList[i].node[1]] - nodeList[edgeList[i].node[0]]).magnitude();

                    // direction
                    Point center = Util::scene()->edges->at(edgeList[i].marker)->center();
                    int direction = (((nodeList[edgeList[i].node[0]].x-center.x)*(nodeList[edgeList[i].node[1]].y-center.y) -
                                      (nodeList[edgeList[i].node[0]].y-center.y)*(nodeList[edgeList[i].node[1]].x-center.x)) > 0) ? 1 : -1;

                    double angle = direction * theta * chordShort / chord;

                    QDomElement eleArc = doc.createElement("arc");
                    eleArc.setAttribute("v1", edgeList[i].node[0]);
                    eleArc.setAttribute("v2", edgeList[i].node[1]);
                    eleArc.setAttribute("angle", QString("%1").arg(rad2deg(angle)));

                    eleCurves.appendChild(eleArc);
                }
            }
        }

        // move nodes (arcs)
        for (int i = 0; i<edgeList.count(); i++)
        {
           // assert(edgeList[i].marker >= 0); // markers changed to marker - 1, check...
            if (edgeList[i].marker != -1)
            {
                // curve
                if (Util::scene()->edges->at(edgeList[i].marker)->angle > 0.0)
                {
                    // angle
                    Point center = Util::scene()->edges->at(edgeList[i].marker)->center();
                    double pointAngle1 = atan2(center.y - nodeList[edgeList[i].node[0]].y,
                                               center.x - nodeList[edgeList[i].node[0]].x) - M_PI;

                    double pointAngle2 = atan2(center.y - nodeList[edgeList[i].node[1]].y,
                                               center.x - nodeList[edgeList[i].node[1]].x) - M_PI;

                    nodeList[edgeList[i].node[0]].x = center.x + Util::scene()->edges->at(edgeList[i].marker)->radius() * cos(pointAngle1);
                    nodeList[edgeList[i].node[0]].y = center.y + Util::scene()->edges->at(edgeList[i].marker)->radius() * sin(pointAngle1);

                    nodeList[edgeList[i].node[1]].x = center.x + Util::scene()->edges->at(edgeList[i].marker)->radius() * cos(pointAngle2);
                    nodeList[edgeList[i].node[1]].y = center.y + Util::scene()->edges->at(edgeList[i].marker)->radius() * sin(pointAngle2);
                }
            }
        }
    }

    // nodes
    for (int i = 0; i<nodeList.count(); i++)
    {
        QDomElement eleVertex = doc.createElement("vertex");
        eleVertex.setAttribute("i", i);
        eleVertex.setAttribute("x", QString("%1").arg(nodeList[i].x));
        eleVertex.setAttribute("y", QString("%1").arg(nodeList[i].y));

        eleVertices.appendChild(eleVertex);
    }

    // elements
    for (int i = 0; i<elementList.count(); i++)
    {
        if (elementList[i].isUsed)
        {
            QDomElement eleElement = doc.createElement(QString("domain:%1").arg(elementList[i].isTriangle() ? "triangle" : "quad"));
            eleElement.setAttribute("v1", elementList[i].node[0]);
            eleElement.setAttribute("v2", elementList[i].node[1]);
            eleElement.setAttribute("v3", elementList[i].node[2]);
            if (!elementList[i].isTriangle())
                eleElement.setAttribute("v4", elementList[i].node[3]);
            eleElement.setAttribute("i", i);
            eleElement.setAttribute("marker", QString("%1").arg(abs(elementList[i].marker)));

            eleElements.appendChild(eleElement);
        }
    }

    // find edge neighbours

    // for each vertex list elements that it belogns to
    QList<QSet<int> > vertexElements;
    vertexElements.reserve(nodeList.count());
    for (int i = 0; i < nodeList.count(); i++)
        vertexElements.push_back(QSet<int>());

    for (int i = 0; i < elementList.count(); i++)
    {
        if (elementList[i].isUsed)
        {
            for(int elemNode = 0; elemNode < (elementList[i].isTriangle() ? 3 : 4); elemNode++)
            {
                vertexElements[elementList[i].node[elemNode]].insert(i);
            }
        }
    }

    for (int i = 0; i < edgeList.count(); i++)
    {
        if (edgeList[i].isUsed && edgeList[i].marker != -1)
        {
            QSet<int> neighbours = vertexElements[edgeList[i].node[0]];
            neighbours.intersect(vertexElements[edgeList[i].node[1]]);
            assert((neighbours.size() > 0) && (neighbours.size() <= 2));
            edgeList[i].neighElem[0] = neighbours.values()[0];
            if(neighbours.size() == 2)
                edgeList[i].neighElem[1] = neighbours.values()[1];
        }
    }


    // subdomains
    foreach(FieldInfo* fieldInfo, Util::scene()->fieldInfos())
    {
        QDomElement eleSubdomain = doc.createElement("subdomain");
        eleSubdomains.appendChild(eleSubdomain);
        eleSubdomain.setAttribute("name", fieldInfo->fieldId());

        QDomElement eleSubElements = doc.createElement("elements");
        eleSubdomain.appendChild(eleSubElements);

        for (int i = 0; i<elementList.count(); i++)
        {
            if (elementList[i].isUsed && (Util::scene()->labels->at(elementList[i].marker)->getMarker(fieldInfo) != SceneMaterialContainer::getNone(fieldInfo)))
            {
                QDomElement eleSubElement = doc.createElement("i");
                eleSubElements.appendChild(eleSubElement);
                QDomText number = doc.createTextNode(QString::number(i));
                eleSubElement.appendChild(number);
            }
        }

        QDomElement eleBoundaryEdges = doc.createElement("boundary_edges");
        eleSubdomain.appendChild(eleBoundaryEdges);
        QDomElement eleInnerEdges = doc.createElement("inner_edges");
        eleSubdomain.appendChild(eleInnerEdges);

        for (int i = 0; i < edgeList.count(); i++)
        {
            QDomElement eleEdge = doc.createElement("i");
            QDomText number = doc.createTextNode(QString::number(i));
            eleEdge.appendChild(number);

            //assert(edgeList[i].marker >= 0);
            if (edgeList[i].isUsed && edgeList[i].marker != -1)
            {
                int numNeighWithField = 0;
                for(int neigh_i = 0; neigh_i < 2; neigh_i++)
                {
                    int neigh = edgeList[i].neighElem[neigh_i];
                    if(neigh != -1)
                    {
                        if(Util::scene()->labels->at(elementList[neigh].marker)->getMarker(fieldInfo)
                                != SceneMaterialContainer::getNone(fieldInfo))
                            numNeighWithField++;
                    }
                }

                // edge has boundary condition prescribed for this field
                bool hasFieldBoundaryCondition = (Util::scene()->edges->at(edgeList[i].marker)->getMarker(fieldInfo)
                                                  != SceneBoundaryContainer::getNone(fieldInfo));

                if(numNeighWithField == 1)
                {
                    // edge is on "boundary" of the field, should have boundary condition prescribed

                    // todo: instead of assert put some warning and stop calculation, as if boundary condition is not assigned at all
                    assert(hasFieldBoundaryCondition);

                    eleBoundaryEdges.appendChild(eleEdge);
                }
                else if(numNeighWithField == 2)
                {
                    // todo: we could enforce not to have boundary conditions prescribed inside:
                    // assert(!hasFieldBoundaryCondition);

                    eleInnerEdges.appendChild(eleEdge);
                }
            }
        }
    }


    nodeList.clear();
    edgeList.clear();
    elementList.clear();

    fileNode.close();
    fileEdge.close();
    fileEle.close();
    fileNeigh.close();

    // save to file
    QFile file(tempProblemFileName() + ".xml");
    if (!file.open(QIODevice::WriteOnly))
        return false;

    QTextStream out(&file);
    doc.save(out, 4);

    file.waitForBytesWritten(0);
    file.close();

    // set system locale
    setlocale(LC_NUMERIC, plocale);

    return true;
}