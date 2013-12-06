#include <QList>
#include <QString>
#include <QObject>

#include "../agros2d-library/meshgenerator.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"
#include "../hermes2d/include/function/exact_solution.h"
#include "../agros2d-library/hermes2d/solver.h"
#include "bem_matrix.h"
#include "mesh.h"

// ToDo: Do it by pointers
Element::Element(Node *a, Node *b, Node *c)
{
    m_nodes.append(a);
    m_nodes.append(b);
    m_nodes.append(c);
    m_gravity = (*a + *b + *c) / 3;
}

Element::Element(QList<Node *> points)
{
    m_nodes = points;
    int i = 0;
    foreach (Node * point, m_nodes) {
        m_gravity = m_gravity + (* point);
        i++;
    }
    m_gravity = m_gravity / i;

    m_area = 0;
    m_f = 0;
}


bool Element::containsPoint(Node x)
{
    double cp1 = 0;
    double cp2 = 0;
    for(int i = 0; i < 3; i++)
    {
        int j = (i+1) % 3;
        int k = (i+2) % 3;
        cp1 = ((*m_nodes.at(j) - *m_nodes.at(i)) % (x - *m_nodes.at(i)));
        cp2 = ((*m_nodes.at(j) - *m_nodes.at(i)) % (*m_nodes.at(k) - *m_nodes.at(i)));
        if((cp1 < 0) || (cp2 < 0))
            return false;
    }

    return true;
}

double Element::value(Node p)
{
    double area = araea();
    double x1 = this->m_nodes.at(0)->x;
    double y1 = this->m_nodes.at(0)->y;
    double x2 = this->m_nodes.at(1)->x;
    double y2 = this->m_nodes.at(1)->y;
    double x3 = this->m_nodes.at(2)->x;
    double y3 = this->m_nodes.at(2)->y;
    double N1 = 0.5 * ((x2 * y3 - x3 * y2) + (y2 - y3) * p.x + (x3 - x2) * p.y) / area;
    double N2 = 0.5 * ((x3 * y1 - x1 * y3) + (y3 - y1) * p.x + (x1 - x3) * p.y) / area;
    double N3 = 0.5 * ((x1 * y2 - x2 * y1) + (y1 - y2) * p.x + (x2 - x1) * p.y) / area;
    double result = N1 * nodeValues[0] +  N2 * nodeValues[1]  + N3 * nodeValues[2];
    return result;
}

Segment::Segment(Node * firstNode, Node * secondNode)
{
    m_nodes.append(firstNode);
    m_nodes.append(secondNode);
    m_gravity.x = (firstNode->x + secondNode->x) / 2;
    m_gravity.y = (firstNode->y + secondNode->y) / 2;
    m_length = sqrt(pow(firstNode->x - secondNode->x, 2) + pow(firstNode->y - secondNode->y, 2));
    m_geometricOrder = 1;
    m_approximationOrder = 0;
    m_derivation = 0;
    m_value = 0;
}

bool Segment::isLyingPoint(Node point)
{

    double dx = lastNode().x - firstNode().x;
    double dy = lastNode().y - firstNode().y;

    Node sp = firstNode();

    double t = ((point.x - sp.x)*dx + (point.y - sp.y)*dy);

    if (t < 0.0)
        t = 0.0;
    else if (t > (dx * dx + dy * dy))
        t = 1.0;
    else
        t /= (dx * dx + dy * dy);

    Node p(sp.x + t*dx, sp.y + t*dy);

    Node dv = point - p;
    return (dv.magnitudeSquared() < EPS_ZERO);
}


double Segment::distanceOf(Node point)
{
    double dx = lastNode().x - firstNode().x;
    double dy = lastNode().y - firstNode().y;

    Node sp = firstNode();

    double t = ((point.x - sp.x)*dx + (point.y - sp.y)*dy);

    if (t < 0.0)
        t = 0.0;
    else if (t > (dx * dx + dy * dy))
        t = 1.0;
    else
        t /= (dx * dx + dy * dy);

    Node p(sp.x + t*dx, sp.y + t*dy);
    Node dv = point - p;

    return dv.magnitude();
}

double Segment::parametricCoordinate(Node node)
{

    double dx = lastNode().x - firstNode().x;
    double dy = lastNode().y - firstNode().y;
    double x = node.x - firstNode().x;
    double y = node.y - firstNode().y;

    if(abs(x) > EPS_ZERO)
        return 2 * (x / dx) - 1;
    else
        if(abs(y) > EPS_ZERO)
            return 2 * (y / dy) - 1;
        else
            return -1;
}

bool Node::operator !=(const Node &vec) const
{
    return (!almostEqualRelAndAbs(vec.x, x, POINT_ABS_ZERO, POINT_REL_ZERO)
            || !almostEqualRelAndAbs(vec.y, y, POINT_ABS_ZERO, POINT_REL_ZERO));
}

bool Node::operator ==(const Node &vec) const
{
    return (almostEqualRelAndAbs(vec.x, x, POINT_ABS_ZERO, POINT_REL_ZERO)
            && almostEqualRelAndAbs(vec.y, y, POINT_ABS_ZERO, POINT_REL_ZERO));
}

Node Node::rotate(double angle)
{   double cos_angle = cos(angle);
    double sin_angle = sin(angle);
    double tmp_x = x * cos_angle -  y * sin_angle;
    y = x * sin_angle +  y * cos_angle;
    x = tmp_x;
    return (*this);
}
