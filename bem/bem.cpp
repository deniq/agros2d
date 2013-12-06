#include <QTextStream>
#include <QTime>
#include <QtTest/QTest>
#include "util.h"
#include "util/global.h"

#include "../agros2d-library/scene.h"
#include "../agros2d-library/scenemarker.h"
#include "../agros2d-library/scenebasic.h"
#include "../agros2d-library/scenenode.h"
#include "../agros2d-library/sceneedge.h"
#include "../agros2d-library/scenelabel.h"
#include "../agros2d-library/hermes2d/module.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"
#include "../hermes2d/include/function/exact_solution.h"

#include <cblas.h>
#include <lapacke.h>

#include "bem.h"
#include "mesh.h"
#include "quad_tables.h"

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{    
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
    m_polyOrder = 0;
}


void Bem::addPhysics()
{ 
    QList<Segment> edgeComponents;
    Hermes::Hermes2D::Element *e;
    int j = 0;
    for(int j = 0; j < Agros2D::scene()->edges->items().count(); j++)
    {
        SceneBoundary *boundary = Agros2D::scene()->edges->items().at(j)->marker(m_fieldInfo);
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            double value= boundary->value(boundaryType.id()).number();
            bool isEssential = (boundaryType.essential().count() > 0);

            int nElement = 0;
            for_all_active_elements(e, m_mesh)
            {
                nElement++;
                QList<Node *> points;
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    Node * point = new Node(e->vn[i]->x, e->vn[i]->y);
                    points.append(point);

                    if(e->vn[i]->bnd == 1)
                    {
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Node * firstPoint = new Node(e->vn[i]->x, e->vn[i]->y);
                            Node * secondPoint = new Node(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y);

                            Node * fP;
                            int index = -1;
                            for(int l = 0; l < mesh.m_nodes.count(); l++)
                            {
                                if(*mesh.m_nodes[l] == *firstPoint)
                                {
                                    index = l;
                                    break;
                                }
                            }

                            if(index == -1)
                            {

                                firstPoint->globalIndex = mesh.m_nodes.count();
                                mesh.m_nodes.append(firstPoint);
                                fP = mesh.m_nodes[mesh.m_nodes.count() - 1];
                            }
                            else
                            {
                                fP = mesh.m_nodes[index];
                            }

                            Node * sP;
                            index = -1;
                            for(int l = 0; l < mesh.m_nodes.count(); l++)
                            {
                                if(*mesh.m_nodes[l] == *secondPoint)
                                {
                                    index = l;
                                    break;
                                }
                            }

                            if(index == -1)
                            {
                                secondPoint->globalIndex = mesh.m_nodes.count();
                                mesh.m_nodes.append(secondPoint);
                                sP =  mesh.m_nodes[mesh.m_nodes.count() - 1];
                            }
                            else
                            {
                                sP = mesh.m_nodes[index];
                            }

                            Segment segment(fP, sP);
                            segment.m_mesh = & mesh;
                            segment.setElement(e);
                            segment.setValue(value);
                            segment.setEssential(isEssential);
                            segment.setEdgeId(j);
                            edgeComponents.append(segment);
                        }
                    }
                }
                Element element(points);
                element.setArea(e->get_area());

                mesh.m_elements.append(element);
                // ToDo: read material properties

            }
            mesh.m_nElement = nElement;
        }
    }

    for(int i = 0; i < edgeComponents.count(); i++)
    {
        mesh.m_segments.append((edgeComponents.at(i)));
    }


    if(m_polyOrder == 1)
    {
        for(int i = 0; i < mesh.m_segments.count(); i++)
        {

            Segment & segment = mesh.m_segments[i];
            if (segment.isEssential())
            {
                segment.firstNode().isEssential = true;
                segment.firstNode().value = segment.value();
                segment.lastNode().isEssential = true;
                segment.lastNode().value = segment.value();
            }
            else
            {
                segment.firstNode().normalDerivation = segment.derivation();
                segment.lastNode().normalDerivation = segment.derivation();
            }

            mesh.m_nodes.append(&segment.firstNode());
            mesh.m_nodes.append(&segment.lastNode());
            segment.m_points.append(&segment.firstNode());
            segment.m_points.append(&segment.lastNode());
        }
    }


    if(m_polyOrder == 0)
    {
        mesh.m_nodes.clear();

        for(int i = 0; i < mesh.m_segments.count(); i++)
        {
            Segment & segment = mesh.m_segments[i];
            Node * node = new Node(segment.gravity());
            if (segment.isEssential())
            {
                node->isEssential = true;
                node->value = segment.value();
            }
            else
            {
                node->normalDerivation = segment.derivation();
            }

            node->globalIndex = i;

            mesh.m_nodes.append(node);
            segment.m_points.append(node);
        }
    }
}

QString Bem::toString()
{
    QString output = "";

    foreach(Segment component, mesh.m_segments)
    {
        output += "Edge: ";
        output += QString::number(component.id());
        output +=  "\n";
        output += "Boundary conditiom type: ";
        output += QString::number(component.isEssential());
        output +=  "\n";
        output += "Boundary conditiom value: ";
        output += QString::number(component.value());
        output +=  "\n";
        output += "Edge components: \n";
        output += "First Point:";
        output += " ";
        output += QString::number(component.firstNode().x);
        output += " ";
        output += QString::number(component.firstNode().y);
        output += "\n";
        output += "Second Point:";
        output += " ";
        output += QString::number(component.lastNode().x);
        output += " ";
        output += QString::number(component.lastNode().y);
        output += "\n";

    }
    return output;
}


Node Bem::globalCoordinates(double xi, Segment segment)
{
    Node v;
    int n = segment.geometricOrder();
    BemVector Sf(n);


    Sf = shapeFunction(n, xi);

    QList<Node> points;
    if(n == 1);
    points.append(segment.firstNode());
    points.append(segment.lastNode());


    for(int i = 0; i <= n; i++)
    {
        v = v + points[i] * Sf(i);
    }
    return v;
}


BemVector Bem::shapeFunction(int polyOrder, double xi)
{     
    BemVector Ni = BemVector(polyOrder + 1);
    if (polyOrder == 0)
    {
        Ni(0) = 1;
        return Ni;
    }

    Ni(0) = 0.5 * (1 - xi);
    Ni(1) = 0.5 * (1 + xi);

    if (polyOrder == 1)
        return Ni;

    Ni(2) = 1.0 - xi * xi;
    Ni(0) = Ni(0) - 0.5 * Ni(2);
    Ni(2) = Ni(1) - 0.5 * Ni(2);
    return Ni;
}

BemVector Bem::shapeFunctionDerivative(int polyOrder, double xi)
{    
    BemVector Dn = BemVector(polyOrder);

    Dn(0) = -0.5;
    Dn(1) = 0.5;

    if((polyOrder == 1) || (polyOrder == 0))
        return Dn;

    Dn(2) = -2.0 * xi;
    Dn(0) = Dn(0) - 0.5 * Dn(2);
    Dn(1) = Dn(1) - 0.5 * Dn(2);
    return Dn;
}

double Bem::jacobian(int polyOrder, double xi, Segment segment)
{
    double dGamma = 0;

    if((polyOrder == 1) || (polyOrder == 0))
    {
        double x = (segment.lastNode().x - segment.firstNode().x);
        double y = (segment.lastNode().y - segment.firstNode().y);

        dGamma = 0.5 * sqrt(x*x + y*y);
    }

    return dGamma;
}

Node Bem::normalVector(double xi, Segment segment)
{
    int polyOrder = segment.geometricOrder();
    BemVector Dn = shapeFunctionDerivative(polyOrder, xi);
    Node dVec;
    Node n;

    QList<Node> points;
    if(polyOrder == 1);
    {
        points.append(segment.firstNode());
        points.append(segment.lastNode());
    }

    for(int i = 0; i < polyOrder + 1; i++)
    {
        dVec = dVec  + points[i] * Dn(i);
    }
    n.x = dVec.y;
    n.y = - dVec.x;
    n = n * 1 / jacobian(polyOrder, xi, segment);
    return n;
}


template<typename Scalar>
BemSolution<Scalar>::BemSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{    
    this->mesh = mesh;
}


template<typename Scalar>
Scalar BemSolution<Scalar>::value(double x, double y) const
{
    return m_bem->getValue(x, y);
}

template<typename Scalar>
void BemSolution<Scalar>::derivatives(double x, double y, Scalar &dx, Scalar &dy) const
{
    dx = 0;
    dy = 0;
}

template<typename Scalar>
Hermes::Hermes2D::MeshFunction<Scalar>* BemSolution<Scalar>::clone() const
{
    if(this->sln_type == Hermes::Hermes2D::HERMES_SLN)
        return Hermes::Hermes2D::Solution<Scalar>::clone();
    BemSolution<Scalar>* sln = new BemSolution<Scalar>(this->mesh);
    sln->setSolver(m_bem);
    return sln;
}

double Bem::kernel_length(Node refNode, Segment segment, double xi)
{
    return 1;
}


double Bem::kernel_laplace2D(Node refNode, Segment segment, double xi)
{
    double r = globalCoordinates(xi, segment).distanceOf(refNode);
    return 1.0 / (2 * M_PI) * log(1/r);
}


double Bem::kernel_laplace2D_derivation(Node refNode, Segment segment, double xi)
{
    Node r = globalCoordinates(xi, segment) - refNode;
    Node n = normalVector(xi, segment);
    double rNorm = globalCoordinates(xi, segment).distanceOf(refNode);
    return - 1.0 / (2 * M_PI) * (n.x *  r.x  + n.y * r.y) / (rNorm * rNorm);
}

void Bem::solve()
{            
    int order = 7;
    int n = mesh.m_segments.count();

    BemMatrix dU(n, n);
    BemMatrix dT(n, n);
    dU.clear();
    dT.clear();

    // Loop over all nodes
    Node node;
    for (int i = 0; i < n; i++)
    {
        node = * mesh.m_nodes[i];

        // Loop over all elements
        for (int j = 0; j < n; j++)
        {
            Segment segment = mesh.m_segments[j];
            // Loop over element nodes
            for(int k = 0; k <= m_polyOrder; k++)
            {
                int index = segment.m_points[k]->globalIndex;
                for(int l = 0; l <= order; l++)
                {
                    double  jac;
                    jac =jacobian(m_polyOrder, gaussCoords[order][l], segment);
                    if(i == index)
                    {
                        double xi;
                        double dxdb = 2;
                        if(m_polyOrder == 0)
                        {
                            xi = 1.0;
                        }
                        else
                        {
                            if(k == 0)
                            {
                                xi = 2 * gaussLeguerreCoords[order][l] - 1.0;
                                dxdb = 2.0;
                            }
                            else
                            {
                                xi = 1.0 - 2.0 * gaussLeguerreCoords[order][l];
                                dxdb = 2.0;
                            }
                        }
                        dU(i, index) +=  - 1 / ( 2 * M_PI) * log(segment.length() / 2) * shapeFunction(m_polyOrder, gaussCoords[order][l])(k) *  jac * gaussWeights[order][l] + 1 / (2 * M_PI) * dxdb * shapeFunction(m_polyOrder, xi)(k) *  jac * gaussLeguerreWeights[order][l];
                        dT(i, index) = 0.5;

                    } else
                    {
                        dU(i, index) += shapeFunction(m_polyOrder, gaussCoords[order][l])(k) * kernel_laplace2D(node, segment, gaussCoords[order][l]) * jac * gaussWeights[order][l];
                        dT(i, index) += shapeFunction(m_polyOrder, gaussCoords[order][l])(k) * kernel_laplace2D_derivation(node, segment, gaussCoords[order][l]) * jac * gaussWeights[order][l];
                    }
                }
            }
        }
    }

    BemVector bp(n);
    bp.clear();


    // Rearranging matrices
    BemMatrix A(n, n);
    BemMatrix C(n,n);
    BemVector rsv(n);

    for (int j = 0; j < n; j++)
    {
        if(mesh.m_nodes[j]->isEssential)
        {
            for(int i = 0; i < n; i++)
            {
                A(i, j) =  - dU(i, j);
                C(i, j) =  - dT(i, j);
            }
            rsv(j) = mesh.m_nodes[j]->value;
        }

        else
        {
            for(int i = 0; i < n; i++)
            {
                A(i, j) =     dT(i, j);
                C(i, j) =     dU(i, j);
            }
            rsv(j) = mesh.m_nodes[j]->normalDerivation;
        }
    }

    // ToDo: Poisson - right side vector

    qDebug() << dU.toString();
    qDebug() << dT.toString();

    BemVector b(n);
    b = C * rsv + bp;

    BemVector results(n);
    results = A.solve(b);
    qDebug() << results.toString();

    for (int i = 0; i < n; i++)
    {
        if(mesh.m_nodes[i]->isEssential)
        {
            mesh.m_nodes[i]->normalDerivation = results(i);

        }
        else
        {
            mesh.m_nodes[i]->value = results(i);
        }
    }

    for (int i = 0; i < n; i++)
    {
        qDebug() << mesh.m_nodes[i]->toString();
        qDebug() << "Value" <<  mesh.m_nodes[i]->value;
        qDebug() << mesh.m_nodes[i]->normalDerivation;
    }

    domainSolution();
}

void Bem::domainSolution()
{

    for(int i = 0; i < mesh.m_elements.count(); i++)
    {
        Element & element = mesh.m_elements[i];
        for(int i = 0; i < 3; i++)
        {
            element.nodeValues[i] = potential(element.m_nodes.at(i)->x, element.m_nodes.at(i)->y);
        }
    }
}

double Bem::getValue(double x, double y)
{
    double result = 0;
    for(int i = 0; i < mesh.m_elements.count(); i++)
    {
        Node p(x,y);
        if(mesh.m_elements[i].containsPoint(p))
        {
            result = mesh.m_elements[i].value(p);
        }
    }

    return result;
}

double Bem::potential(double x, double y)
{
    int order = 7;
    int n = mesh.m_segments.count();

    double u = 0;
    Node p(x, y);
    n = mesh.m_segments.count();
    for(int i = 0; i < n; i++)
    {
        Segment edge = mesh.m_segments[i];
        if(p.distanceOf(edge.firstNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return edge.m_points[0]->value;

            return edge.firstNode().value;
        }

        if(p.distanceOf(edge.lastNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return edge.m_points[0]->value;

            return edge.lastNode().value;
        }

        if(p.distanceOf(edge.gravity()) < EPS_ZERO)
        {
            return edge.m_points[0]->value;
        }

        if(edge.distanceOf(p) < EPS_ZERO)
        {
            double xi  = edge.parametricCoordinate(p);
            double result = 0;

            for(int k = 0; k <= m_polyOrder; k++)
            {
                result +=  shapeFunction(m_polyOrder, xi)(k) * edge.m_points[k]->value;
            }
            return result;
        }

        QList<Node> nodes;
        if(m_polyOrder == 0)
            nodes.append(* edge.m_points[0]);
        else
        {
            nodes.append(edge.firstNode());
            nodes.append(edge.lastNode());
        }
        //        double H = 0;
        //        double G = 0;
        //        Node integ = integral(p, a, b);

        double dU = 0;
        double dT = 0;

        double delta_u = 0;
        for(int k = 0; k <= m_polyOrder; k++)
        {
            for(int l = 0; l <= order; l++)
            {
                double  jac =jacobian(m_polyOrder, gaussCoords[order][l], edge);
                dT += shapeFunction(m_polyOrder, gaussCoords[order][l])(k) * kernel_laplace2D_derivation(p, edge, gaussCoords[order][l]) * jac * gaussWeights[order][l];
                dU += shapeFunction(m_polyOrder, gaussCoords[order][l])(k) * kernel_laplace2D(p, edge, gaussCoords[order][l]) * jac * gaussWeights[order][l];
            }
            delta_u =   - dT * edge.m_points[k]->value + dU * edge.m_points[k]->normalDerivation;
            u = u + delta_u;
            dT = 0;
            dU = 0;
        }

        //        G = integ.y; // Analytic solution
        //        H = integ.x; // Analytic solution
    }

    //    for(int i = 0; i < m; i++)
    //    {
    //        double R = sqrt((p.x - m_elements[i].gravity().x) * (p.x - m_elements[i].gravity().x) - (p.y - m_elements[i].gravity().y) * (p.y - m_elements[i].gravity().y));
    //        u = u - 1 / (2 * M_PI) * m_elements[i].f() * log(R) * m_elements[i].araea();
    //    }

    // qDebug("Time elapsed: %f ms", t.elapsed());
    return u;
}


Node Bem::integral(Node v, Node pa, Node pb)
{
    // center
    Node center = (pa + pb) / 2;

    // shift of center of the edge on the position [0, 0]
    Node at = pa - center;
    Node bt = pb - center;

    // shift of reference point
    Node vt = v - center;

    // rotation
    double phi = - atan2(bt.y, bt.x);
    vt.rotate(phi);
    at.rotate(phi);
    bt.rotate(phi);

    double H;
    double G;

    if(v.distanceOf(center) < EPS_ZERO)
    {
        // singular point
        H = 0;
        double m = abs(bt.x);
        G = -1 / ( M_PI)  * (m * log(m) - m);
        return Node(H, G);
    } else
    {

        double x = 0;
        double y = 0;
        double a = (vt.x - at.x);
        double b = (vt.x - bt.x);

        if(abs(vt.y) < EPS_ZERO)
        {

            if((abs(vt.x - at.x) < EPS_ZERO) && (abs(vt.x - bt.x) < EPS_ZERO))
            {
                assert(0);
            }

            if(abs(b) < EPS_ZERO)
            {
                x = (a > 0) ? M_PI_2 : - M_PI_2;
                y = M_PI_4;
                G =  - ( a * (-2 + log(a * a)))/(4 * M_PI);
            } else
                if(abs(a) < EPS_ZERO)
                {
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    x = M_PI_4;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI);
                }
                else
                {
                    x = (a > 0) ? M_PI_2 : - M_PI_2;
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI) - (a * (-2 + log(a * a)))/(4 * M_PI);
                }
        }
        else
        {
            x = atan( a / vt.y);
            y = atan( b / vt.y);
            G = (2 * vt.y * y + b * (-2 + log(vt.y * vt.y + b * b)))/(4 * M_PI) -
                    (2 * vt.y * x + a * (-2 + log(vt.y * vt.y + a * a)))/(4 * M_PI);
        }
        H = (y - x)  / (2 * M_PI);
    }

    return Node(H, G);
}

template class BemSolution<double>;
