#include <QTextStream>

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

Element::Element(Node a, Node b, Node c)
{
    m_nodes.append(a);
    m_nodes.append(b);
    m_nodes.append(c);
    m_gravity = (a + b + c) / 3;
}

Element::Element(QList<Node> nodes)
{
    m_nodes = nodes;
    int i = 0;
    foreach (Node node, m_nodes) {
        m_gravity = m_gravity + node;
        i++;
    }
    m_gravity = m_gravity / i;

    m_area = 0;
    m_f = 0;
}

EdgeComponent::EdgeComponent(Node firstNode, Node secondNode)
{
    m_firstNode = firstNode;
    m_secondNode = secondNode;
    m_gravity(0) = (firstNode(0) + secondNode(0)) / 2;
    m_gravity(1) = (firstNode(1) + secondNode(1)) / 2;
    m_length = sqrt(pow(firstNode(0) - secondNode(0), 2) + pow(firstNode(1) - secondNode(1), 2));
    m_derivation = 0;
    m_value = 0;
}

bool EdgeComponent::isLyingNode(Node node)
{

    double dx = m_secondNode(0) - m_firstNode(0);
    double dy = m_secondNode(1) - m_firstNode(1);

    Node sp = m_firstNode;

    double t = ((node(0) - sp(0))*dx + (node(1) - sp(1))*dy);

    if (t < 0.0)
        t = 0.0;
    else if (t > (dx * dx + dy * dy))
        t = 1.0;
    else
        t /= (dx * dx + dy * dy);

    Node p(sp(0) + t*dx, sp(1) + t*dy);

    Node dv = node - p;
    return (dv.length() * dv.length() < EPS_ZERO);
}


Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{    
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
}


void Bem::addPhysics()
{

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
                QList<Node> nodes;
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    Node node(e->vn[i]->x, e->vn[i]->y, e->vn[i]->id);
                    nodes.append(node);
                    if(e->vn[i]->bnd == 1)
                    {
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Node firstNode(e->vn[i]->x, e->vn[i]->y, e->vn[i]->id);
                            Node secondNode(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y, e->vn[e->next_vert(i)]->id);

                            EdgeComponent component(firstNode, secondNode);
                            component.m_element = e;
                            component.m_value = value;
                            component.m_isEssential = isEssential;
                            component.m_edgeID = j;
                            m_edgeComponents.append(component);
                        }
                    }
                }
                Element element(nodes);
                element.setArea(e->get_area());
                m_elements.append(element);
                // ToDo: read material properties

            }
            // ToDo: improve
            m_nElement = nElement;
        }
    }
}

QString Bem::toString()
{
    QString output = "";

    foreach(EdgeComponent component, m_edgeComponents)
    {
        output += "Edge: ";
        output += QString::number(component.m_edgeID);
        output +=  "\n";
        output += "Boundary conditiom type: ";
        output += QString::number(component.m_isEssential);
        output +=  "\n";
        output += "Boundary conditiom value: ";
        output += QString::number(component.m_value);
        output +=  "\n";
        output += "Edge components: \n";
        output += "First node:";
        output += QString::number(component.firstNode().id());
        output += " ";
        output += QString::number(component.firstNode()(0));
        output += " ";
        output += QString::number(component.firstNode()(1));
        output += "\n";
        output += "Second node:";
        output += QString::number(component.secondNode().id());
        output += " ";
        output += QString::number(component.secondNode()(0));
        output += " ";
        output += QString::number(component.secondNode()(1));
        output += "\n";

    }
    return output;
}


template<typename Scalar>
BemSolution<Scalar>::BemSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{    
    this->mesh = mesh;
}


template<typename Scalar>
Scalar BemSolution<Scalar>::value(double x, double y) const
{
    return m_bem->potential(x, y);
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

void Bem::solve()
{        
    int n = m_edgeComponents.count();
    BemMatrix H(n, n);
    BemMatrix G(n, n);
    for (int i = 0; i < n; i++)
    {
        Node v(m_edgeComponents.at(i).gravity());
        for (int j = 0; j < n; j++)
        {
            Node a = m_edgeComponents[j].firstNode();
            Node b = m_edgeComponents[j].secondNode();
            Node integ = integral(v, a, b);

            H(i, j) = integ(0);
            G(i, j) = integ(1);

            if (i == j)
                H(i, j) = H(i, j) + 0.5;
        }
    }

    BemVector bp(n);
    for(int i = 0; i < n; i++ )
    {
        for(int j = 0; j < m_elements.count(); j++)
        {
            double R = sqrt((m_edgeComponents[i].gravity()(0) - m_elements[j].gravity()(0)) * (m_edgeComponents[i].gravity()(0) - m_elements[j].gravity()(0)) +
                            (m_edgeComponents[i].gravity()(1) - m_elements[j].gravity()(1)) * (m_edgeComponents[i].gravity()(1) - m_elements[j].gravity()(1)));
            bp(i) += 1/(2 * M_PI) * m_elements[j].f() * log(R) * m_elements[j].araea();
        }
    }

    // Rearranging matrices
    BemMatrix A(n, n);
    BemMatrix C(n,n);
    BemVector rsv(n);

    for (int i = 0; i < n; i++)
    {
        if(m_edgeComponents[i].m_isEssential)
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) = - G(j, i);
                C(j, i) = H(j, i);
            }
            rsv(i) = m_edgeComponents[i].m_value;
        }

        else
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) = - H(j, i);
                C(j, i) =   G(j, i);
            }
            rsv(i) = m_edgeComponents[i].m_derivation;
        }
    }

    // ToDo: Poisson - right side vector

    BemVector b = C * rsv + bp;
    BemVector results = A.solve(b);
    // qDebug() << results.toString();

    for (int i = 0; i < n; i++)
    {
        if(m_edgeComponents[i].m_isEssential)
        {
            m_edgeComponents[i].m_derivation = results(i);
        }
        else
        {
            m_edgeComponents[i].m_value = results(i);
        }
    }
}

double Bem::potential(double x, double y)
{    
    double u = 0;
    Node p(x, y);
    int n = m_edgeComponents.count();
    int m = m_elements.count();


    for(int i = 0; i < n; i++)
    {
        if(m_edgeComponents[i].isLyingNode(p))
        {
            return m_edgeComponents[i].m_value;
        }
        Node a = m_edgeComponents[i].firstNode();
        Node b = m_edgeComponents[i].secondNode();
        Node integ = integral(p, a, b);

        double H = integ(0);
        double G = integ(1);

        u = u - H * m_edgeComponents[i].m_value - G * m_edgeComponents[i].m_derivation;
    }

    //    for(int i = 0; i < m; i++)
    //    {
    //        double R = sqrt((p(0) - m_elements[i].gravity()(0)) * (p(0) - m_elements[i].gravity()(0)) - (p(1) - m_elements[i].gravity()(1)) * (p(1) - m_elements[i].gravity()(1)));
    //        u = u - 1 / (2 * M_PI) * m_elements[i].f() * log(R) * m_elements[i].araea();
    //    }

    return u;
}

Node Bem::integral(Node v, Node a, Node b)
{    
    // center
    Node center = (a + b) / 2;

    // shift of center of the edge on the position [0, 0]
    Node at = a - center;
    Node bt = b - center;

    // shift of reference point
    Node vt = v - center;

    // rotation
    double phi = - atan2(bt(1), bt(0));
    vt.rotate(phi);
    at.rotate(phi);
    bt.rotate(phi);


    double H;
    double G;

    if(v.distanceOf(center) < EPS_ZERO)
    {
        // singular point
        H = 0;
        double m = abs(bt(0));
        G = -1 / (2 * M_PI) * 2 * (m * log(m) - m);
    } else
    {
        // regular point
        H = atan((vt(0) - bt(0)) / vt(1)) / (2 * M_PI) - atan((vt(0) - at(0))/vt(1)) / (2 * M_PI);
        G = (2 * vt(1) * atan((vt(0) - bt(0))/vt(1)) + (vt(0) - bt(0))*(-2 + log(vt(1) * vt(1) + (vt(0) - bt(0)) * (vt(0) - bt(0)))))/(4 * M_PI) -
                (2 * vt(1) * atan((vt(0) - at(0))/vt(1)) + (vt(0) - at(0))*(-2 + log(vt(1) * vt(1) + (vt(0) - at(0)) * (vt(0) - at(0)))))/(4 * M_PI);
    }

    return Node(H, G);
}

template class BemSolution<double>;
